# Rutina_2
# 🧬 Análisis del Contenido GC y Genes Significativos en *Bartonella henselae*

Este proyecto realiza un análisis genómico de *Bartonella* utilizando Python.  
El objetivo es calcular el **contenido GC (%)** por gen, analizar la **longitud génica** y determinar qué genes son **significativamente diferentes** según su tamaño.

---

## 🔍 Descripción general

- **Lectura de archivos** `.fna` (genoma) y `.gff` (anotaciones de genes).  
- **Cálculo del contenido GC** con una función que mide el porcentaje de guanina y citosina por gen.  
- **Análisis de significancia**: se marcan como “significativos” los genes cuya longitud supera ±2 desviaciones estándar (p < 0.05).  
- **Visualizaciones**: se emplean gráficos de Seaborn y Matplotlib para representar la distribución de %GC, longitud y relaciones entre ambas variables.

---

## 📊 Gráficos principales

1. Histograma del %GC  
2. Histograma de longitudes  
3. Boxplot de longitud por rango de GC  
4. Violinplot del %GC  
5. Scatterplot longitud vs GC (con significancia)  
6. Diagrama de barras con conteo de genes por rango de GC

## ¿Que necesitas?
Para ejecutar el programa correctamente necesitas:

- **Python 3.x**
- Las siguientes librerías instaladas:
  ```bash
  pip install pandas seaborn matplotlib statistics
- Ademas los archivos
- GCF_019930925.1_ASM1993092v1_genomic.fna
- genomic.gff
