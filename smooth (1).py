import numpy as np
import h5py
import time

inicio = time.time()

# Datos geométricos (m)
r = 0.019
l = 0.019
c = 0.000075

# Datos de condiciones de operación
rpm = 1800
w = rpm*2*np.pi/60 # angular velocity
mu = 0.032

E = [
    0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 0.91,
    0.92, 0.93, 0.94, 0.95, 0.96,
    0.97, 0.98, 0.985
    ]

# Datos de criterios de convergencia
tr = 1e-8 # tolerance reynolds
ib = 20 # iterations bisection
tb = 1e-4 # tolerance nisection
tol_phi = 1e-3 # tolerance attitude angle
ip = 5 # iterations attitude angle

# Datos de tamaño de la cuadrícula
N_theta = 360
N_z = 100

theta = np.linspace(0, 2*np.pi, N_theta).reshape((N_theta, 1)) # coordinates in theta [radians]

z = np.linspace(0, l, N_z).reshape((1, N_z)) # coordinates in z [m]

Z = z/l # dimensionless coordiantes in z

def film_thickness_smooth(c, E, theta, divisions, phi):
    
    h = c*(1+E*np.cos(theta-phi)) # film thickness
    
    h = np.repeat(h, divisions, axis=1) # matrix film thickness
    
    return h

def reynolds_steady_smooth(r, c, l, w, mu, E, theta, Z, tr, phi):
    
    P = np.zeros((np.shape(theta)[0], np.shape(Z)[1])) # boundary condition
    
    # film thickness [m]
    h = film_thickness_smooth(c, E, theta, np.shape(Z)[1], phi)
    
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

    h_1 = np.where(np.logical_and(dp_theta == 0 , p == 0), 0, h) # only take values where there is fluid
    h_2 = np.where(np.logical_and(dp_theta == 0 , p == 0), h, 0) # only take values where there is cavitation
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

    ff_3 = np.trapz(np.trapz(np.where(h_2 != 0, h_c*mu*w*r**2/(h_2)**2, 0), x=theta[:,0], axis=0), x=Z*l) # frictional force in cavitation zone
    
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

# Crear del archivo HDF5
file_name = 'smooth2400.h5'

with h5py.File(file_name, 'w') as filename:

    # Almacenar datos geométricos
    filename.create_dataset('radius_journal', data=r)
    filename.create_dataset('length_bearing', data=l)
    filename.create_dataset('clearance', data=c)

    # Almacenar datos de criterios de convergencia
    filename.create_dataset('error_tolerance_reynolds', data=tr)
    #filename.create_dataset('iterations_bisection', data=ib)
    #filename.create_dataset('error_tolerance_bisection', data=tb)

    # Almacenar datos de tamaño de la cuadrícula
    filename.create_dataset('divisions_theta', data=N_theta)
    filename.create_dataset('divisions_z', data=N_z)

    # Añadir atributos a los conjuntos de datos y grupos
    filename['radius_journal'].attrs['unit'] = 'm'
    filename['length_bearing'].attrs['unit'] = 'm'
    filename['clearance'].attrs['unit'] = 'm'
    filename['error_tolerance_reynolds'].attrs['unit'] = '-'
    #filename['iterations_bisection'].attrs['unit'] = '-'
    #filename['error_tolerance_bisection'].attrs['unit'] = '-'
    filename['divisions_theta'].attrs['unit'] = '-'
    filename['divisions_z'].attrs['unit'] = '-'

    filename.create_dataset('angular_velocity', data=w)
    filename['angular_velocity'].attrs['unit'] = 'radians/s'

    filename.create_dataset('coordinate_theta', data=np.rad2deg(theta).transpose())
    filename['coordinate_theta'].attrs['unit'] = '°'
    filename.create_dataset('coordinate_z', data=z.transpose())
    filename['coordinate_z'].attrs['unit'] = 'm'
    filename.create_dataset('coordinate_Z', data=Z.transpose())
    filename['coordinate_Z'].attrs['unit'] = '-'
    filename.create_dataset('dynamic_viscosity', data = mu)
    filename['dynamic_viscosity'].attrs['unit'] = 'Pa*s'

    for e in E:

        subgroup = f'Eccentricity ratio_{e}'
        subgroup = filename.create_group(subgroup)

        p, P, lcc, COF, phi_1, h, H, tau, q, Po, Po_l, fr, ft, ff_1, ff_2, ff_3, ff, er = reynolds_steady_smooth(r, c, l, w, mu, e, theta, Z, tr, 0)
        p, P, lcc, COF, phi_2, h, H, tau, q, Po, Po_l, fr, ft, ff_1, ff_2, ff_3, ff, er = reynolds_steady_smooth(r, c, l, w, mu, e, theta, Z, tr, np.deg2rad(phi_1))

        if abs(phi_2) > tol_phi:
    
            for i in range(ip):

                phi_3 = phi_1 + phi_2
                p, P, lcc, COF, phi_4, h, H, tau, q, Po, Po_l, fr, ft, ff_1, ff_2, ff_3, ff, er = reynolds_steady_smooth(r, c, l, w, mu, e, theta, Z, tr, np.deg2rad(phi_3))

                if abs(phi_4) < tol_phi:
                    phi = phi_3
                    ep = phi_4
                    break
                else:
                    phi_5 = phi_1 - phi_2
                    p, P, lcc, COF, phi_6, h, H, tau, q, Po, Po_l, fr, ft, ff_1, ff_2, ff_3, ff, er = reynolds_steady_smooth(r, c, l, w, mu, e, theta, Z, tr, np.deg2rad(phi_5))

                    if abs(phi_6) < tol_phi:
                        phi = phi_5
                        ep = phi_6
                        break
        
                    if abs(phi_4) < abs(phi_6):
                        phi_1 = phi_3
                        phi_2 = phi_4
                    else:
                        phi_1 = phi_2
                        phi_2 = phi_5
        else:
            phi = phi_1
            ep = phi_2

        if e == 0.1:
            print('{:<25} {:<25} {:<25}'.format('Eccentricity Ratio', 'Error Attitude Angle', 'Error Pressure'))
        print('{:<25} {:<25} {:<25}'.format(e, ep[0], er))
    
        subgroup.create_dataset('pressure', data=p)
        subgroup['pressure'].attrs['unit'] = 'Pa'
        subgroup.create_dataset('maximum_pressure', data=np.max(p))
        subgroup['maximum_pressure'].attrs['unit'] = 'Pa'
        subgroup.create_dataset('dimensionless_pressure', data=P)
        subgroup['dimensionless_pressure'].attrs['unit'] = '-'
        subgroup.create_dataset('load_carrying_capacity', data=np.ndarray.item(lcc))
        subgroup['load_carrying_capacity'].attrs['unit'] = 'N'
        subgroup.create_dataset('coefficient_of_friction', data=np.ndarray.item(COF))
        subgroup['coefficient_of_friction'].attrs['unit'] = '-'
        subgroup.create_dataset('attitude_angle', data=np.ndarray.item(phi))
        subgroup['attitude_angle'].attrs['unit'] = '°'
        subgroup.create_dataset('film_thickness', data=h)
        subgroup['film_thickness'].attrs['unit'] = 'm'
        subgroup.create_dataset('minimum_film_thickness', data=np.min(h))
        subgroup['minimum_film_thickness'].attrs['unit'] = 'm'
        subgroup.create_dataset('dimensionless_film_thickness', data=H)
        subgroup['dimensionless_film_thickness'].attrs['unit'] = '-'
        subgroup.create_dataset('torque', data=np.ndarray.item(tau))
        subgroup['torque'].attrs['unit'] = 'N*m'
        subgroup.create_dataset('leakage_flow_rate', data=q)
        subgroup['leakage_flow_rate'].attrs['unit'] = 'm^3/s'
        subgroup.create_dataset('eccentricity_ratio', data=e)
        subgroup['eccentricity_ratio'].attrs['unit'] = '-'
        subgroup.create_dataset('power', data=np.ndarray.item(Po))
        subgroup['power'].attrs['unit'] = 'W'
        subgroup.create_dataset('power_loss', data=Po_l)
        subgroup['power_loss'].attrs['unit'] = 'W'
        subgroup.create_dataset('radial_force', data=np.ndarray.item(fr))
        subgroup['radial_force'].attrs['unit'] = 'N'
        subgroup.create_dataset('tangential_force', data=np.ndarray.item(ft))
        subgroup['tangential_force'].attrs['unit'] = 'N'
        subgroup.create_dataset('frictional_force_due_change_pressure', data=np.ndarray.item(ff_1))
        subgroup['frictional_force_due_change_pressure'].attrs['unit'] = 'N'
        subgroup.create_dataset('frictional_force_in_fluid_zone', data=np.ndarray.item(ff_2))
        subgroup['frictional_force_in_fluid_zone'].attrs['unit'] = 'N'
        subgroup.create_dataset('frictional_force_in_cavitation_zone', data=np.ndarray.item(ff_3))
        subgroup['frictional_force_in_cavitation_zone'].attrs['unit'] = 'N'
        subgroup.create_dataset('frictional_force', data=np.ndarray.item(ff))
        subgroup['frictional_force'].attrs['unit'] = 'N'
        subgroup.create_dataset('error_reynolds', data=er)
        subgroup['error_reynolds'].attrs['unit'] = '-'

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