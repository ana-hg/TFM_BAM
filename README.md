# TFM_BAM
Repositorio de scripts y datos del TFM realizado por Brandon Díaz Martín, Martín Huqiao Astoreca Martín y Ana Hernández García para el máster de Bioinformática de la Universidad Internacional de la Rioja (curso 2024-2025)

Los archivos raw descargados de las bases de datos pueden descargarse utilizando el siguiente enlace: 
https://drive.google.com/drive/folders/13V5oSTvUpOY7Y3YPG0uJ6lnwNWhUL1v0?usp=sharing

## Scripts
En este directorio se encuentran los scripts empleados durante el trabajo:
- **RNA_completo.R**:  es un script en R en el que se realiza un análisis integral de datos de expresión génica de RNA para construir y evaluar modelos de aprendizaje automático. El objetivo principal es predecir la etapa FIGO (un sistema de clasificación para el cáncer) de los tumores, basándose en el perfil de expresión génica de los pacientes con cáncer de endometrio (proyecto TCGA-UCEC).
El proceso se divide en:
  1. Adquisición de datos: El script utiliza el paquete TCGAbiolinks para descargar y preparar datos brutos de expresión génica de RNA-Seq y la información clínica correspondiente del repositorio GDC.  
  2. Control de calidad y selección: Filtra los datos clínicos para seleccionar los pacientes adecuados y filtra los datos de expresión génica para eliminar genes con baja variabilidad. También se aplican métodos para normalizar los datos y corregir posibles efectos de lotes (sva).  
  3. Preprocesamiento: Reduce la dimensionalidad de los datos utilizando Análisis de Componentes Principales (PCA) y Elastic Net. Luego, divide los datos en conjuntos de entrenamiento y prueba.
  4. Entrenamiento y evaluación de modelos: Entrena y optimiza varios modelos de aprendizaje automático (SVM, XGboost, Random Forest, LDA y FDA) para clasificar las muestras de tumores en sus respectivas etapas FIGO.
  5. Análisis de resultados: Evalúa el rendimiento de los modelos, genera matrices de confusión y produce tablas con métricas de evaluación clave (exactitud, precisión, sensibilidad, etc.), que luego se guardan como imágenes para un informe.
  
- **maf_completo.R**: es un script en R un análisis completo de datos genómicos para construir y evaluar modelos de aprendizaje automático. El objetivo principal es predecir la etapa FIGO (un sistema de clasificación para el cáncer) de los tumores en pacientes con cáncer de útero, utilizando mutaciones genéticas como características predictivas.
El proceso se divide en cuatro etapas principales:
  1. Adquisición de datos: Descarga datos clínicos y de mutaciones (archivos MAF) del proyecto TCGA-UCEC de los Institutos Nacionales de Salud de EE. UU.
  2. Control de calidad y selección: Filtra y limpia los datos, eliminando pacientes y mutaciones de baja calidad. También equilibra las clases de las etapas FIGO para asegurar que los modelos no se sesguen.
  3. Preprocesamiento: Transforma los datos de mutación en formatos adecuados para el aprendizaje automático, creando cuatro conjuntos de características diferentes: conteos de mutaciones por gen, conteos de mutaciones por variante, y versiones binarias (sí/no) de ambos. El script usa Elastic Net para reducir el número de características.
  4. Entrenamiento y evaluación de modelos: Entrena varios modelos de aprendizaje automático, incluyendo SVM, XGboost, Random Forest, LDA y FDA, para clasificar las muestras. Finalmente, evalúa el rendimiento de los modelos con métricas como la exactitud, la precisión y la puntuación F1, y genera tablas y gráficos para visualizar los resultados.
 

- **batch_effect.py**: es un script de Python diseñado para detectar y visualizar efectos de lote (batch effects) en datos genéticos, específicamente datos de microarrays. Los efectos de lote son variaciones no deseadas en los datos que surgen de la forma en que las muestras fueron procesadas, en lugar de variaciones biológicas reales.

## ML_models
### DNASeq
Este directorio contiene los modelos entrenados a partir de datos de **DNA-seq** para el TFM, organizados en cuatro representaciones distintas de los datos (carpetas `1`, `2`, `3`, `4`).
- **1/** → Genes mutados (binario):  
  Indica si un paciente presenta mutación en un gen (`1`) o no (`0`).  

- **2/** → Genes mutados (conteo):  
  Indica el **número de mutaciones** que presenta cada gen en cada paciente.  

- **3/** → Tipos de mutación (binario):  
  Cada columna representa un tipo de mutación en un gen. Se indica solo **presencia/ausencia** de ese tipo de mutación en cada paciente.  

- **4/** → Tipos de mutación (conteo):  
  Similar al anterior, pero en este caso se indica el **número de veces** que aparece cada tipo de mutación en cada paciente.

### RNASeq
Este directorio contiene los modelos entrenados a partir de datos de **RNA-seq** para el TFM, organizados en 3 grupos distintos de los datos (carpetas `1`, `2`, `3`) dependiendo de cómo se hizo la correción de*Batch effect*.
- **1/** → RNASeq sin corrección de *Batch effect*  

- **2/** → RNASeq con corrección de *Batch effect* hecha en R con el paquete Combat

- **3/** → RNASeq con corrección de *Batch effect* hecha en Python con el paquete Combat y NeuroCombat

#### Notas
- Cada carpeta contiene los modelos entrenados (`.rds`) correspondientes a esa representación de los datos.  
- El sufijo del archivo (`fda`, `lda`, `rf`, `svm`, `xgb`) indica el algoritmo utilizado.  

