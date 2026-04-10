import numpy as np
import pandas as pd
import h5py
import time
from sys import argv

# Registra el tiempo de inicio
inicio = time.time()

ruta_excel = "circle.xlsx" 
datos_excel = pd.read_excel(ruta_excel)

# Datos geométricos
r = 0.019
l = 0.019
c = 0.000075

# Datos de condiciones de operación
rpm = 2400
w = rpm*2*np.pi/60 # angular velocity
mu = 0.032

# Datos de criterios de convergencia
tr = 1e-8
ib = 20
tb = 1e-4
tol_phi = 1e-3 # tolerance attitude angle
ip = 5 # iterations attitude angle

# Datos de tamaño de la cuadrícula
N_theta = 360
N_z = 100
Gt = 20

# Unidades correspondientes
unidades_r = "m"
unidades_l = "m"
unidades_c = "m"
unidades_rpm = "RPM"
unidades_mu = "Pa*s"
unidades_tr = "-"
unidades_ib = "-"
unidades_tb = "-"
unidades_N_theta = "-"
unidades_N_z = "-"

# Valores de Eccentricity ratio
valores_eccentricity = [
    0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 0.91,
    0.92, 0.93, 0.94, 0.95, 0.96,
    0.97, 0.98, 0.985
    ]

run = int(argv[0][:-3])

# Crear del archivo HDF5
with h5py.File(f"cir_{run}.h5", "w") as filename:

    grupo_runs = filename.create_group(f"Run_{run}")

    theta_1 = datos_excel.iloc[run-1, 1]
    theta_f = datos_excel.iloc[run-1, 2]
    Nt_theta = datos_excel.iloc[run-1, 3]
    Nt_z = datos_excel.iloc[run-1, 4]
    d = datos_excel.iloc[run-1, 5]
    rt = datos_excel.iloc[run-1, 6]

    grupo_runs.create_dataset("initial_position", data=theta_1)
    grupo_runs["initial_position"].attrs["unit"] = ["°"]
    grupo_runs.create_dataset("fraction_final_position", data=theta_f)
    grupo_runs["fraction_final_position"].attrs["unit"] = ["-"]
    grupo_runs.create_dataset("textures_theta", data=Nt_theta)
    grupo_runs["textures_theta"].attrs["unit"] = ["-"]
    grupo_runs.create_dataset("textures_z", data=Nt_z)
    grupo_runs["textures_z"].attrs["unit"] = ["-"]
    grupo_runs.create_dataset("depth", data=d)
    grupo_runs["depth"].attrs["unit"] = ["m"]
    grupo_runs.create_dataset("radious_texture", data=rt)
    grupo_runs["radious_texture"].attrs["unit"] = ["m"]

    grupo_eccentricity = grupo_runs.create_group("Eccentricity_ratios")

    # Crear los subgrupos para cada valor de "Eccentricity ratio"
    for valor_e in valores_eccentricity:
        nombre_subgrupo = f"Eccentricity_ratio_{valor_e}"
        grupo_eccentricity_valor = grupo_eccentricity.create_group(nombre_subgrupo)

theta = np.linspace(0, 2*np.pi, N_theta) # coordinates in theta [radians]
z = np.linspace(0, l, N_z) # coordinates in z [m]

theta_2 = (360-theta_1)*theta_f+theta_1

# center of circles in theta
C_theta = np.linspace(theta_1, theta_2, (Nt_theta + 2))
C_theta = C_theta[1:-1]
C_theta = np.deg2rad(C_theta)*r

# center of circles in z
C_z = np.linspace(0, l, (Nt_z + 2))
C_z = C_z[1:-1]

# define ranges where there is texture
r_theta = np.array([C_theta - rt, C_theta + rt])
r_z = np.array([C_z - rt, C_z + rt])

# define logical mask to omit points between ranges
theta = theta*r
mask_theta = np.ones_like(theta, dtype=bool)
mask_z = np.ones_like(z, dtype=bool)

for i in range(np.shape(r_theta)[1]):
    mask_theta = mask_theta & np.logical_or(np.less(theta, r_theta[0,i]), np.greater(theta, r_theta[1,i]))

for i in range(np.shape(r_z)[1]):
    mask_z = mask_z & np.logical_or(np.less(z, r_z[0,i]), np.greater(z, r_z[1,i]))

# apply logical mask to coodinates
theta = np.compress(mask_theta, theta)
theta = np.concatenate((theta, r_theta[0, :], r_theta[1, :]))
z = np.compress(mask_z, z)
z = np.concatenate((z, r_z[0, :], r_z[1, :]))

# create refinement
for i in range(np.shape(r_theta)[1]):
    gr_theta = np.linspace(r_theta[0,i], r_theta[1,i], Gt)
    theta = np.concatenate((theta, gr_theta))

for i in range(np.shape(r_z)[1]):
    gr_z = np.linspace(r_z[0,i], r_z[1,i], Gt)
    z = np.concatenate((z, gr_z))

# sort coordinates and eliminate duplicates
theta = np.unique(theta)
theta = np.reshape(theta, (len(theta), 1))
z = np.unique(z)
z = np.reshape(z, (1, len(z)))
theta = theta/r
Z = z/l # dimensionless coordiantes in z

with h5py.File(f"cir_{run}.h5", "a") as filename:
    
    filename[f"Run_{run}"].create_dataset("coordinate_theta", data=np.rad2deg(theta))
    filename[f"Run_{run}"]["coordinate_theta"].attrs["unit"] = "°"
    filename[f"Run_{run}"].create_dataset("coordinate_z", data=z)
    filename[f"Run_{run}"]["coordinate_z"].attrs["unit"] = "m"
    filename[f"Run_{run}"].create_dataset("coordinate_Z", data=Z)
    filename[f"Run_{run}"]["coordinate_Z"].attrs["unit"] = "-"
    filename[f"Run_{run}"].create_dataset("angular_velocity", data=w)
    filename[f"Run_{run}"]["angular_velocity"].attrs["unit"] = "radians/s"
    filename[f"Run_{run}"].create_dataset("final_position", data=theta_2)
    filename[f"Run_{run}"]["final_position"].attrs["unit"] = "radians"

def film_thickness_circle(r, c, E, theta, z, rt, C_theta, C_z, d, phi):
    
    h = c*(1+E*np.cos(theta-phi)) # film thickness
    h = np.tile(h, (1, np.shape(z)[1])) # matrix film thickness

    # convert units
    theta = theta*r

    # film thikness
    for i in range(np.shape(theta)[0]):
        for j in range(np.shape(z)[1]):
            for k in range(len(C_theta)):
                for l in range(len(C_z)):
                    if ((theta[i] - C_theta[k])**2 + (z[0, j] - C_z[l])**2) <= rt**2:
                 
                        h[i,j] = h[i,j] + d
    return h

def reynolds_steady_circle(r, c, l, w, mu, E, theta, Z, rt, C_theta, C_z, d, tr, phi):
    
    P = np.zeros((np.shape(theta)[0], np.shape(Z)[1])) # boundary condition
    
    # film thickness [m]
    h = film_thickness_circle(r, c, E, theta, z, rt, C_theta, C_z, d, phi)
    
    H = h/c

    # coefficients
    a_1 = (H[1:-1, 1:-1] + H[2:, 1:-1])**3 / (8*(theta[1:-1] - theta[2:])**2) \
        + (H[1:-1, 1:-1] + H[:-2, 1:-1])**3 / (8*(theta[1:-1] - theta[:-2])**2) \
        + (r**2/l**2)*((H[1:-1, 1:-1] + H[1:-1, 2:])**3) / (8*(Z[0, 1:-1] - Z[0, 2:])**2) \
        + (r**2/l**2)*((H[1:-1, 1:-1] + H[1:-1, :-2])**3) / (8*(Z[0, 1:-1] - Z[0, :-2])**2)
    
    a_2 = (H[1:-1, 1:-1] + H[2:, 1:-1])**3 / (8*(theta[1:-1] - theta[2:])**2)
    
    a_3 = (H[1:-1, 1:-1] + H[:-2, 1:-1])**3 / (8*(theta[1:-1] - theta[:-2])**2)
    
    a_4 = (r**2/l**2)*((H[1:-1, 1:-1] + H[1:-1, 2:])**3) / (8*(Z[0, 1:-1] - Z[0, 2:])**2)
    
    a_5 = (r**2/l**2)*((H[1:-1, 1:-1] + H[1:-1, :-2])**3) / (8*(Z[0, 1:-1] - Z[0, :-2])**2)
    
    a_6 = (H[2:, 1:-1] - H[:-2, 1:-1]) / (2*(theta[2:] - theta[:-2]))
    
    a_2_a_1 = a_2/a_1
    a_3_a_1 = a_3/a_1
    a_4_a_1 = a_4/a_1
    a_5_a_1 = a_5/a_1
    a_6_a_1 = a_6/a_1
    
    # pressure [Pa]
    while 1:

        previous_k = np.sum(P)

        P_1 = a_2_a_1*P[2:,1:-1]
        P_2 = a_3_a_1*P[:-2,1:-1]
        P_3 = a_4_a_1*P[1:-1,2:]
        P_4 = a_5_a_1*P[1:-1,:-2]
        P[1:-1,1:-1] = P_1 + P_2 + P_3 + P_4 - a_6_a_1

        P[P < 0] = 0
        
        actual_k = np.sum(P)

        er = np.abs(actual_k - previous_k)/np.abs(actual_k)
        
        if er < tr:
            break
    
    # Variables
    p = 6*P*mu*w*r**2/c**2 # pressure
    
    fr = np.trapz(np.trapz(p*r*np.sin(theta), x=theta, axis=0), x=Z*l) # radial force [N]
    
    ft = -np.trapz(np.trapz(p*r*np.cos(theta), x=theta, axis=0), x=Z*l) # tangencial force [N]
    
    lcc = np.sqrt(ft**2 + fr**2) # load carrying capacity [N]

    dp_theta, dp_z = np.gradient(p); # derivative of pressure
    
    h_1 = np.where(np.logical_and(dp_theta == 0, p == 0), 0, h) # only take values where there is fluid
    h_2 = np.where(np.logical_and(dp_theta == 0, p == 0), h, 0) # only take values where there is cavitation
    h_2[:,0] = 0
    h_2[:,-1] = 0

    ff_1 = np.trapz(np.trapz(h_1*dp_theta/2, x=theta[:,0], axis=0), x=Z*l) # frictional force due change in pressure
    ff_2 = np.trapz(np.trapz(np.where(h_1 != 0, mu*w*r**2/h_1, 0), x=theta[:,0], axis=0), x=Z*l) # frictional force in fluid zone

    h_c = np.zeros(h_2.shape) # values where cavitation starts
    
    for i in range(h_2.shape[0]):

        for j in range(h_2.shape[1]):

            if h_2[i,j] != 0:
                
                if h_2[i-1,j] ==0:
                    
                    h_c[i,j] = h_2[i,j]
                else:
                    h_c[i,j] = h_c[i-1,j]

    # firctional force in cavitation zone
    ff_3 = np.trapz(np.trapz(np.where(h_2 != 0, h_c*mu*w*r**2/(h_2)**2, 0), x=theta[:,0], axis=0), x=Z*l)
    
    ff = ff_1 + ff_2 + ff_3 # frictional force

    COF = ff/lcc # coefficient of friction
    
    phi = np.rad2deg(np.arctan(fr/ft)) # attitude angle [°]
    
    tau = ft*r # torque [Nm]

    q_0 = np.trapz(h[:,0]**3*dp_z[:,0]*r/(12*mu), x=theta[:,0], axis=0) # flow rate in z=0

    q_l = -np.trapz(h[:,-1]**3*dp_z[:,-1]*r/(12*mu), x=theta[:,0], axis=0) # flow rate in z=l

    q = abs(q_0) + abs(q_l) # leakage flow rate

    Po = ft*r*w # Power [W]

    Po_l = np.trapz(np.trapz(mu*(w*r)**2*r/h, x=Z*l, axis=1), x=theta[:,0], axis=0) + \
        np.trapz(np.trapz(h**3*np.linalg.norm([dp_z, dp_theta])**2*r/(12*mu), x=Z*l, axis=1), x=theta[:,0], axis=0) # power loss [W]

    return p, P, lcc, COF, phi, h, H, tau, q, Po, Po_l, fr, ft, ff_1, ff_2, ff_3, ff, er
    
# Nombre del archivo HDF5
file_name = f"cir_{run}.h5"

with h5py.File(file_name, "a") as filename:
    eccentricity_ratios_group = filename[f"Run_{run}"]["Eccentricity_ratios"]

    for E in valores_eccentricity:
                    
        subgroup_name = f"Eccentricity_ratio_{E}"
            
        subgroup = eccentricity_ratios_group[subgroup_name]

        p, P, lcc, COF, phi_1, h, H, tau, q, Po, Po_l, fr, ft, ff_1, ff_2, ff_3, ff, er = \
            reynolds_steady_circle(r, c, l, w, mu, E, theta, Z, rt, C_theta, C_z, d, tr, 0)
        p, P, lcc, COF, phi_2, h, H, tau, q, Po, Po_l, fr, ft, ff_1, ff_2, ff_3, ff, er = \
            reynolds_steady_circle(r, c, l, w, mu, E, theta, Z, rt, C_theta, C_z, d, tr, np.deg2rad(phi_1))
                             
        if abs(phi_2) > tol_phi:
    
            for i in range(ip):

                phi_3 = phi_1 + phi_2
                p, P, lcc, COF, phi_4, h, H, tau, q, Po, Po_l, fr, ft, ff_1, ff_2, ff_3, ff, er = \
                    reynolds_steady_circle(r, c, l, w, mu, E, theta, Z, rt, C_theta, C_z, d, tr, np.deg2rad(phi_3))

                if abs(phi_4) < tol_phi:
                    phi = phi_3
                    ep = phi_4
                    break
                else:
                    phi_5 = phi_1 - phi_2
                    p, P, lcc, COF, phi_6, h, H, tau, q, Po, Po_l, fr, ft, ff_1, ff_2, ff_3, ff, er = \
                        reynolds_steady_circle(r, c, l, w, mu, E, theta, Z, rt, C_theta, C_z, d, tr, np.deg2rad(phi_5))

                    if abs(phi_6) < tol_phi:
                        phi = phi_5
                        ep = phi_6
                        break
        
                    if abs(phi_4) < abs(phi_6):
                        phi_1 = phi_3
                        phi_2 = phi_4
                    else:
                        phi_1 = phi_5
                        phi_2 = phi_6
        else:
            phi = phi_1
            ep = phi_2

        if E == 0.1:
            print('{:<25} {:<25} {:<25}'.format('Eccentricity Ratio', 'Error Attitude Angle', 'Error Pressure'))
        print('{:<25} {:<25} {:<25}'.format(E, ep[0], er))
    
        subgroup.create_dataset("pressure", data=p)
        subgroup["pressure"].attrs["unit"] = "Pa"
        subgroup.create_dataset('maximum_pressure', data=np.max(p))
        subgroup['maximum_pressure'].attrs['unit'] = 'Pa'
        subgroup.create_dataset("dimensionless_pressure", data=P)
        subgroup["dimensionless_pressure"].attrs["unit"] = "-"
        subgroup.create_dataset("load_carrying_capacity", data=lcc[0])
        subgroup["load_carrying_capacity"].attrs["unit"] = "N"
        subgroup.create_dataset("coefficient_of_friction", data=COF[0])
        subgroup["coefficient_of_friction"].attrs["unit"] = "-"
        subgroup.create_dataset("attitude_angle", data=phi[0])
        subgroup["attitude_angle"].attrs["unit"] = "°"
        subgroup.create_dataset("film_thickness", data=h)
        subgroup["film_thickness"].attrs["unit"] = "m"
        subgroup.create_dataset('minimum_film_thickness', data=np.min(h))
        subgroup['minimum_film_thickness'].attrs['unit'] = 'm'
        subgroup.create_dataset("dimensionless_film_thickness", data=H)
        subgroup["dimensionless_film_thickness"].attrs["unit"] = "-"
        subgroup.create_dataset("torque", data=tau[0])
        subgroup["torque"].attrs["unit"] = "N*m"
        subgroup.create_dataset("leakage_flow_rate", data=q)
        subgroup["leakage_flow_rate"].attrs["unit"] = "m^3/s"
        subgroup.create_dataset("power", data=Po[0])
        subgroup["power"].attrs["unit"] = "W"
        subgroup.create_dataset("power_loss", data=Po_l)
        subgroup["power_loss"].attrs["unit"] = "W"
        subgroup.create_dataset("radial_force", data=fr[0])
        subgroup["radial_force"].attrs["unit"] = "N"
        subgroup.create_dataset("tangential_force", data=ft[0])
        subgroup["tangential_force"].attrs["unit"] = "N"
        subgroup.create_dataset("frictional_force_due_change_pressure", data=ff_1[0])
        subgroup["frictional_force_due_change_pressure"].attrs["unit"] = "N"
        subgroup.create_dataset("frictional_force_in_fluid_zone", data=ff_2[0])
        subgroup["frictional_force_in_fluid_zone"].attrs["unit"] = "N"
        subgroup.create_dataset("frictional_force_in_cavitation_zone", data=ff_3[0])
        subgroup["frictional_force_in_cavitation_zone"].attrs["unit"] = "N"
        subgroup.create_dataset("frictional_force", data=ff[0])
        subgroup["frictional_force"].attrs["unit"] = "N"
        subgroup.create_dataset("error_reynolds", data=er)
        subgroup["error_reynolds"].attrs["unit"] = "-"
        subgroup.create_dataset("eccentricity_ratio", data=E)
        subgroup["eccentricity_ratio"].attrs["unit"] = "-"

# Registra el tiempo al final
fin = time.time()

# Calcula la diferencia de tiempo
tiempo_total = fin - inicio

# Calcula los días, horas, minutos y segundos
dias = int(tiempo_total // 86400)
horas = int((tiempo_total % 86400) // 3600)
minutos = int((tiempo_total % 3600) // 60)
segundos = tiempo_total % 60

tiempo_formateado = f"{dias}d {horas}h {minutos}m {segundos:.2f}s"

print("Tiempo total:", tiempo_formateado)