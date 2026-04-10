# Optimización Tribológica de Cojinetes mediante Business Intelligence 🚀

Este proyecto aplica herramientas de **Inteligencia de Negocios (BI)** para analizar y optimizar el rendimiento de cojinetes hidrodinámicos mediante la implementación de texturas superficiales láser. A través de un proceso de **ETL** y **Minería de Datos**, se transformaron simulaciones complejas en decisiones estratégicas de ingeniería.

## 📋 Resumen del Proyecto
El objetivo principal es determinar si la adición de micro-texturas en la superficie de los cojinetes mejora su eficiencia energética. Utilizando datos generados por simulaciones basadas en la literatura de **Vlădescu y Valdés**, el proyecto identifica las configuraciones geométricas que minimizan la fricción sin comprometer la capacidad de carga.

## 🛠️ Tecnologías Utilizadas
* **Análisis de Datos:** Microsoft Power BI (Power Query, DAX).
* **Simulación:** Python (Procesamiento de modelos de lubricación elastohidrodinámica).
* **Metodología:** Ciclo de vida de BI (Extracción, Transformación, Carga y Visualización).
* **Generativas:** ClaudeIA

## 🔍 Hallazgos Clave (Minería de Datos)
* **Reducción de Fricción:** Se identificó una disminución del **COF de hasta un 96.82%** en condiciones de alta carga ($E = 0.98$).
* **Escenario Dorado:** La configuración óptima identificada fue la **Familia C (Círculo)** con profundidad de $100.5 \mu m$ y densidad del $60\%$.
* **Umbral Crítico:** Se descubrió que el texturizado es ineficiente en bajas excentricidades ($E < 0.5$), donde el flujo de aceite ya es estable.

## ⚙️ Proceso de Ingeniería de Datos (ETL)
Para garantizar la integridad del análisis, se realizaron las siguientes transformaciones en Power BI:
1. **Normalización de Unidades:** Conversión de escalas científicas ($1e-06$ m) a unidades industriales ($\mu m$ y grados).
2. **Columnas Derivadas:** Creación de segmentaciones por "Familia" de textura y rangos cualitativos de capacidad de carga (LCC).
3. **Limpieza de Outliers:** Filtrado de configuraciones físicamente inviables mediante reglas de validación de negocio.
4. **Métricas DAX:** Implementación de indicadores de comparación dinámica vs. el cojinete liso convencional.

## 📊 Dashboard de Decisiones
El informe final incluye:
- **Análisis de Sensibilidad:** Comparación de formas (Círculo, Elipse, Rectángulo).
- **Mapa de Calor:** Identificación de la densidad de textura óptima.
- **KPIs de Eficiencia:** Ahorro porcentual de energía proyectado.

---
**Autor:** Luciana Velásquez Ciro
**Contexto:** Proyecto de Laboratorio de Materiales / Inteligencia de Negocios
