#***************************************************************************************************************
#----Adquisición de datos----
#***************************************************************************************************************
#install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")
#BiocManager::install("maftools")
#BiocManager::install("sva")



library(TCGAbiolinks)
#library(maftools)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(data.table)
library(tibble)
library(sva)
library(caret)

library(gt)
library(webshot2)

set.seed(123) #reproducibilidad

#****************************************************************************
## ---- Importar datos ----
#****************************************************************************
#Si no se quiere repetir el proceso completo se pueden seleccionar datos de aquí e importarlos
base_dir <- 'C:/Users/Ana/OneDrive - Universidad de Salamanca/Bioinformática/4. TFM/Datos'

#Datos RAW base de datos
maf_df <- fread(file.path(base_dir, "maf/raw/maf_data.csv"), na.strings = c("", " "))
clinical_df <- fread(file.path(base_dir, "maf/raw/clinical_data.csv"), na.strings = c("", " "))



#Datos filtrados
maf_filtrado <- fread(file.path(base_dir, "maf/filtrado/maf_filtrado.csv"), na.strings = c("", " "))
clinical_filtrado <- fread(file.path(base_dir, "maf/filtrado/clinical_filtrado.csv"), na.strings = c("", " "))
clinical_filtrado <- column_to_rownames(clinical_filtrado, var = "patient_id")

#Datos preprocesados
#Sería mejor guardar un rds
##Wider df
mc_gen_num <- fread(file.path(base_dir, "mc_gen_num.csv"), na.strings = c("", " "))
mc_gen_num <- column_to_rownames(mc_gen_num, var = "patient_id")
mc_mut_num <- fread(file.path(base_dir, "mc_mut_num.csv"), na.strings = c("", " "))
mc_mut_num <- column_to_rownames(mc_mut_num, var = "patient_id")
mc_gen_encoded <- fread(file.path(base_dir, "mc_gen_encoded.csv"), na.strings = c("", " "))
mc_gen_encoded <- column_to_rownames(mc_gen_encoded, var = "patient_id")
mc_mut_encoded <- fread(file.path(base_dir, "mc_mut_encoded.csv"), na.strings = c("", " "))
mc_mut_encoded <- column_to_rownames(mc_mut_encoded, var = "patient_id")

#Datos entrenamiento/testing
base_dir <- 'C:/Users/Ana/OneDrive - Universidad de Salamanca/Bioinformática/4. TFM/Datos/maf/elastic_net'
training_data_en <- readRDS(file.path(base_dir,"training_data_en.rds"))
training_labels <- readRDS(file.path(base_dir,"training_labels.rds"))
testing_data_en <- readRDS(file.path(base_dir,"testing_data_en.rds"))
testing_labels <- readRDS(file.path(base_dir,"testing_labels.rds"))

#Modelos ML
svm_mod_4 <- readRDS(file.path(base_dir,"svm_mod_mc_4.rds"))
xgb_mod_4 <- readRDS(file.path(base_dir,"xgb_mod_mc_4.rds"))
rf_mod_4 <- readRDS(file.path(base_dir,"rf_mod_mc_4.rds"))
lda_mod_4 <- readRDS(file.path(base_dir,"lda_mod_mc_4.rds"))
fda_mod_4 <- readRDS(file.path(base_dir,"fda_mod_mc_4.rds"))




#****************************************************************************
## ---- TCGABiolinks ----
#****************************************************************************
#setwd("C:/Users/Ana/TCGA")


### ---- 1. Descargar mutaciones (MAF) ----
query_mut <- GDCquery(
  project = "TCGA-UCEC",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open" #Muy importante = open porque no tenemos permisos para descargar los datos protegidos
)

GDCdownload(query_mut) #Descargamos los datos en la carpeta indicada
maf_df <- GDCprepare(query_mut) #Juntamos toda la información en un único dataframe



### ---- 2. Descargar datos clínicos ----
clinical_df <- GDCquery_clinic(project = "TCGA-UCEC", type = "clinical")


#***************************************************************************************************************
#----Contol de calidad y selección----
#***************************************************************************************************************

###----1. Filtrado datos clínicos----
clinical_df$patient_id <- clinical_df$submitter_id #cambiamos el nombre de la columna para poder unir los dataframes

clinical_filtrado <- clinical_df[
  !is.na(figo_stage) & 
    synchronous_malignancy == "No" & 
    prior_malignancy == "no" & 
    classification_of_tumor == "primary"
]

clinical_filtrado <- clinical_filtrado %>%
  mutate(figo_stage = str_extract(figo_stage, "Stage [IVX]+")) %>%
  filter(figo_stage!= 'Stage IV')

lista_columnas <- c('patient_id', 'figo_stage') #Seleccionamos o figo stage o tumor_grade
clinical_filtrado <- clinical_filtrado %>% select(all_of(lista_columnas))

clinical_filtrado <- drop_na(clinical_filtrado) #Eliminamos los pacientes sin figo_stage

###----2. Filtrado datos mutaciones----
maf_df$patient_id <- substr(maf_df$Tumor_Sample_Barcode, 1, 12) #Generamos una columna con los ids de los pacientes
maf_df$VAF_tumor <- maf_df$t_alt_count / maf_df$t_depth 


#Generamos histogamas de las columnas que queremos filtrar para ver la distribución
library (ggplot2)
library(scales)
library(patchwork)

maf_graficas <- maf_df[
  t_depth < 500 &
    t_alt_count < 200
]

parametros_graficar <- list('t_depth', 't_alt_count', 'VAF_tumor')

lista_graficas <- list()

for (param in parametros_graficar) {
  histograma <- ggplot(maf_graficas, aes(x = .data[[param]])) + 
    geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "cornflowerblue", color = "black") + 
    #La función after_stat(Density) es más recomendable en las verisones más modernas de R que la función ..density.., que ya quedó anticuada. 
    #Esta función permite utilizar la función geom_density en el histograma. 
    geom_density(color = "black", linewidth = 0.5) +
    labs (title = param, x = "Valor", y = "Frecuencia") + 
    theme(legend.position = "none", 
          plot.margin = margin(8,8,8,8),
          plot.title = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 7)) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma) +
    geom_density()
  lista_graficas[[param]] <- histograma
}

conjunto <- wrap_plots(lista_graficas, ncol = 3, nrow = 1, axis_titles = "collect")
conjunto


#Una vez hemos visto la distribución de los datos, realizamos un filtrado
maf_filtrado <- maf_df[
  t_depth >= 20 &
    t_alt_count >= 5 &
    VAF_tumor >= 0.05 & 
    IMPACT != 'LOW' & #Eliminamos las sinónimas y las de regions de splicing que son mutaciones de bajo impacto
    Variant_Classification != 'Silent' &
    Mutation_Status == 'Somatic' & #Nos aseguramos de utilizar sólo mutaciones somáticas
    Feature_type == 'Transcript' #Eliminamos las mutaciones que afectan a regiones reguladoras y no a transcritos
]

table(maf_filtrado$Variant_Classification)

lista_columnas <- c("Hugo_Symbol", "Variant_Classification", "patient_id")

maf_filtrado <- maf_filtrado %>% select(all_of(lista_columnas))



###----4. Selección final pacientes enfermos utilizados----
#Obtenemos la lista de pacientes seleccionados
maf_patients <- data.frame(patient_id = unique(maf_filtrado$patient_id))
clinical_patients <- data.frame(patient_id = clinical_filtrado$patient_id)

patients_id_maf_clinical <- inner_join(maf_patients, clinical_patients, by = 'patient_id')


#Filtramos los diferentes archivos
maf_filtrado <- inner_join(maf_filtrado, patients_id_maf_clinical, by = 'patient_id')
#fwrite(maf_filtrado, file.path(base_dir, "maf/filtrado/maf_filtrado.csv"))
clinical_filtrado <- inner_join(clinical_filtrado, patients_id_maf_clinical, by = 'patient_id')
#fwrite(clinical_filtrado, file.path(base_dir, "maf/filtrado/clinical_filtrado.csv"))



#***************************************************************************************************************
#----Preprocesado----
#***************************************************************************************************************

###----1. Tratamiento de NAs----
#Comprobamos si existen NAs
#any(is.na(clinical_filtrado))
#any(is.na(maf_filtrado))

#Los NAs del clinical ya se eliminaron

#Función para calcular la moda (valor más frecuente)
find_mode <- function(x) { 
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

#Funciones de imputación para datos categóricos y numéricos
impute_na_cat <- function(x) replace(x, is.na(x), find_mode(x, na.rm = TRUE))

if (any(is.na(maf_filtrado))){
  maf_filtrado <- maf_filtrado %>% 
    mutate(across(where(is.character(.)), impute_na_cat))
}


###----2. Pivotar la tabla de mutaciones----

#En vez de pivotar la tabla directamente lo hacemos en dos pasos. Primero las agrupamos y hacemos el encoding
#De esta manera queda como mutación sí/no
#Después pivotamos para que sean categóricas
#Si directamente lo pasamos a wider queda con el número y no con sí/no

####----2.1 Pivotar con número de mutaciones----
maf_wider_gen <- maf_filtrado %>%
  group_by(patient_id, Hugo_Symbol) %>%
  summarise(
    mutation_number = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = mutation_number,
    values_fill = 0
  )

#fwrite(maf_wider_gen, file.path(base_dir, "maf_wider_gen.csv"))


maf_wider_mut <- maf_filtrado %>%
  group_by(patient_id, Hugo_Symbol, Variant_Classification) %>%
  summarise(
    mutation_number = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = c(Hugo_Symbol, Variant_Classification),
    values_from = mutation_number,
    values_fill = 0,
    names_glue = "{Hugo_Symbol}_{Variant_Classification}"
  )

#fwrite(maf_wider_mut, file.path(base_dir, "maf_wider_mut.csv"))

####----2.2 Pivotar con mutación sí/no----

#Agrupamos en función de las características que nos interesan

maf_agrupado_gen <- maf_filtrado %>% 
  group_by(patient_id, Hugo_Symbol) %>%
  summarise(mutation_number = n())


maf_agrupado_mut <- maf_filtrado %>% 
  group_by(patient_id, Hugo_Symbol, Variant_Classification) %>%
  summarise(mutation_number = n())


#Hacemos encoding y pivotamos (un sólo paso)

maf_gen_encoded <- maf_agrupado_gen %>%
  ungroup() %>%
  mutate(value = 1) %>%
  select(patient_id, Hugo_Symbol, value) %>%
  pivot_wider(
    names_from = Hugo_Symbol, 
    values_from = value, 
    values_fill = 0
    )


maf_mut_encoded <- maf_agrupado_mut %>%
  ungroup() %>%
  mutate(value = 1) %>%
  select(patient_id, Hugo_Symbol, value, Variant_Classification) %>%
  pivot_wider(
    names_from = c(Hugo_Symbol, Variant_Classification), 
    values_from = value, 
    values_fill = 0,
    names_glue = "{Hugo_Symbol}_{Variant_Classification}"
    )


###----3. Unión variable predictora----
#mc es de maf_clinical pero para que no se alargue tanto el nobre

mc_gen_num <- merge(clinical_filtrado, maf_wider_gen, by="patient_id", all = FALSE)
#fwrite(mc_gen_num, file.path(base_dir, "mc_gen_num.csv"))
mc_gen_num <- column_to_rownames(mc_gen_num, var = "patient_id")

mc_mut_num <- merge(clinical_filtrado, maf_wider_mut, by="patient_id", all = FALSE)
#fwrite(mc_mut_num, file.path(base_dir, "mc_mut_num.csv"))
mc_mut_num <- column_to_rownames(mc_mut_num, var = "patient_id")



mc_gen_encoded <- merge(clinical_filtrado, maf_gen_encoded, by="patient_id", all = FALSE)
#fwrite(mc_gen_encoded, file.path(base_dir, "mc_gen_encoded.csv"))
mc_gen_encoded <- column_to_rownames(mc_gen_encoded, var = "patient_id")

mc_mut_encoded <- merge(clinical_filtrado, maf_mut_encoded, by="patient_id", all = FALSE)
#fwrite(mc_mut_encoded, file.path(base_dir, "mc_mut_encoded.csv"))
mc_mut_encoded <- column_to_rownames(mc_mut_encoded, var = "patient_id")


###----4. Filtrado de los datos----
####----4.1 Selección de genes----
#Empezamos por los genes
gen_cols <- colnames(mc_gen_encoded[,-1])

#Los niveles son los mismos para los 4 dataframes así los sacamos sólo 1 vez
niveles <- unique(mc_gen_encoded$figo_stage)
genes_validos_por_grupo <- list()

for (nivel in niveles) {
  subset_nivel <- mc_gen_encoded[mc_gen_encoded$figo_stage == nivel, gen_cols]
  p_zeros <- colSums(subset_nivel == 0) / nrow(subset_nivel)
  genes_validos <- names(p_zeros[p_zeros <= 0.95]) #Seleccionamos los que mutan en al menos un 5% de los pacientes de un grupo
  genes_validos_por_grupo[[nivel]] <- genes_validos
}

genes_finales <- Reduce(union, genes_validos_por_grupo)

mc_gen_num_fil <- mc_gen_num[, c("figo_stage",genes_finales)]
mc_gen_encoded_fil <- mc_gen_encoded[, c("figo_stage",genes_finales)]

#Continuamos con las mutaciones
mut_cols <- colnames(mc_mut_encoded[,-1])
mutaciones_validas_por_grupo <- list()

for (nivel in niveles) {
  subset_nivel <- mc_mut_encoded[mc_mut_encoded$figo_stage == nivel, mut_cols]
  p_zeros <- colSums(subset_nivel == 0) / nrow(subset_nivel)
  mutaciones_validas <- names(p_zeros[p_zeros <= 0.95]) #Seleccionamos los que mutan en al menos un 5% de los pacientes de un grupo
  mutaciones_validas_por_grupo[[nivel]] <- mutaciones_validas
}

mutaciones_finales <- Reduce(union, mutaciones_validas_por_grupo)

mc_mut_num_fil <- mc_mut_num[, c("figo_stage",mutaciones_finales)]
mc_mut_encoded_fil <- mc_mut_encoded[, c("figo_stage",mutaciones_finales)]

####----4.2 Selección de pacientes (balancear datos)----
#Como los pacientes son los mismos en todos, determinamos los pacientes con uno de ellos y después los aplicamos
#Primero comprobamos cuántos pacients hay en cada clase
table(clinical_filtrado$figo_stage)
patients_subsets <- list()
sample_size <- 50

#No necesitamos repetir esto para todos
#Lo hago con el de antes de filtrar porque ya tengo la lista de las columnas
for (nivel in niveles) {
  subset_nivel <- mc_gen_num[mc_gen_num$figo_stage == nivel, gen_cols]
  patients <- rownames(subset_nivel)
  if (nivel == 'Stage I' || nivel == 'Stage III') {
    subsample_patients <- sample(patients, size = sample_size)
    patients_subsets[[nivel]] <- subsample_patients
  } else{
    patients_subsets[[nivel]] <- patients
  }
}

patients_subset_total <- Reduce(union, patients_subsets)


mc_gen_num_bal <- mc_gen_num_fil[patients_subset_total, ]
#mc_gen_num_bal_export <- rownames_to_column(mc_gen_num_bal, var = "patient_id")
#fwrite(mc_gen_num_bal_export, file.path(base_dir, "mc_gen_num_bal.csv"))
mc_gen_encoded_bal <- mc_gen_encoded_fil[patients_subset_total, ]
#mc_gen_encoded_bal_export <- rownames_to_column(mc_gen_encoded_bal, var = "patient_id")
#fwrite(mc_gen_encoded_bal_export, file.path(base_dir, "mc_gen_encoded_bal.csv"))
mc_mut_num_bal <- mc_mut_num_fil[patients_subset_total, ]
#mc_mut_num_bal_export <- rownames_to_column(mc_mut_num_bal, var = "patient_id")
#fwrite(mc_mut_num_bal_export, file.path(base_dir, "mc_mut_num_bal.csv"))
mc_mut_encoded_bal <- mc_mut_encoded_fil[patients_subset_total, ]
#mc_mut_encoded_bal_export <- rownames_to_column(mc_mut_encoded_bal, var = "patient_id")
#fwrite(mc_mut_encoded_bal_export, file.path(base_dir, "mc_mut_encoded_bal.csv"))

### ----8. División de los datos ----
target_var <- 'figo_stage'
#Generamos la lista de pacientes seleccionados para el entrenamiento. Va a a ser para todos la misma lista así que sólo lo hacemos una vez
training_samples_mc <- mc_gen_num_bal[[target_var]] %>% createDataPartition(p = 0.8, list = FALSE) #Estratificamos por grupo. 80/20

#Empezamos con los datos de training
#Primero generamos una lista con los datos disponibles
data <- list(mc_gen_num_bal, mc_gen_encoded_bal, mc_mut_num_bal, mc_mut_encoded_bal)

#Creamos listas de dataframes vacíos
training_data_names <- c('train_mc_gen_num_bal', 'train_mc_gen_enc_bal', 'train_mc_mut_num_bal', 'train_mc_mut_enc_bal')
training_data <- setNames(lapply(training_data_names, function(x) data.frame()), training_data_names)
training_labels_names <- c('train_mc_gen_num_bal_label', 'train_mc_gen_enc_bal_label', 'train_mc_mut_mut_bal_label', 'train_mc_mut_enc_bal_label')
training_labels <- setNames(lapply(training_labels_names, function(x) data.frame()), training_labels_names)

testing_data_names <- c('test_mc_gen_num_bal', 'test_mc_gen_enc_bal', 'test_mc_mut_num_bal', 'test_mc_mut_enc_bal')
testing_data <- setNames(lapply(testing_data_names, function(x) data.frame()), testing_data_names)
testing_labels_names <- c('test_mc_gen_num_bal_label', 'test_mc_gen_enc_bal_label', 'test_mc_gen_mut_bal_label', 'test_mc_mut_enc_bal_label')
testing_labels <- setNames(lapply(testing_labels_names, function(x) data.frame()), testing_labels_names)


#Completamos la información de los dataframes
for (i in 1:4) {
  training_data[[i]] <- data[[i]][training_samples_mc,]
  training_labels[[i]] <- as.factor(training_data[[i]]$figo_stage)
  levels(training_labels[[i]]) <- c("Stage_I", "Stage_II", "Stage_III")
  training_data[[i]] <- training_data[[i]][,-1]
}

for (i in 1:4) {
  testing_data[[i]] <- data[[i]][-training_samples_mc,]
  testing_labels[[i]] <- as.factor(testing_data[[i]]$figo_stage)
  levels(testing_labels[[i]]) <- c("Stage_I", "Stage_II", "Stage_III")
  testing_data[[i]] <- testing_data[[i]][,-1]
}

#saveRDS(training_data, file = file.path(base_dir, "training_data.rds"))
#saveRDS(training_labels, file = file.path(base_dir, "training_labels.rds"))
#saveRDS(testing_data, file = file.path(base_dir, "testing_data.rds"))
#saveRDS(testing_labels, file = file.path(base_dir, "testing_labels.rds"))


#***************************************************************************************************************
#----Reducción del número de variables----
#***************************************************************************************************************

#Utilizamos ElasticNet para reducir la dimensionalidad de los datos (sólo con los de entrenamiento)
#Generamos listas para los datos reducidos
training_data_names <- c('train_mc_gen_num_bal_en', 'train_mc_gen_enc_bal_en', 'train_mc_mut_num_bal_en', 'train_mc_mut_enc_bal_en')
training_data_en <- setNames(lapply(training_data_names, function(x) data.frame()), training_data_names)

testing_data_names <- c('test_mc_gen_num_bal_en', 'test_mc_gen_enc_bal_en', 'test_mc_mut_num_bal_en', 'test_mc_mut_enc_bal_en')
testing_data_en <- setNames(lapply(testing_data_names, function(x) data.frame()), testing_data_names)

library(glmnet)

for (i in 1:4) {
  en_data <- training_data[[i]]
  en_labels <- training_labels[[i]]
  
  #Modelo de Elastic Net
  en_model <- cv.glmnet(as.matrix(en_data), en_labels, 
                        alpha = 0.5, 
                        family = "multinomial", 
                        type.measure = "class",
                        nfolds = 10)
  
  #Seleccionamos el mejor valor
  best_lambda <- en_model$lambda.min
  
  
  #Extraer coeficientes para cada clase
  coef_list <- coef(en_model, s = best_lambda)
  
  #Filtrar variables seleccionadas (no cero) en al menos una clase
  selected_genes <- unique(unlist(
    lapply(coef_list, function(coef) {
      vars <- rownames(coef)[which(coef != 0)]
      vars[vars != "(Intercept)"]
    })
  ))
  
  #Seleccionamos los datos de training y testing
  training_data_en[[i]] <- training_data[[i]][,selected_genes]
  testing_data_en[[i]] <- testing_data[[i]][,selected_genes]
}

#saveRDS(training_data_en, file = file.path(base_dir, "training_data_en.rds"))
#saveRDS(testing_data_en, file = file.path(base_dir, "testing_data_en.rds"))


#***************************************************************************************************************
#----Machine Learning----
#***************************************************************************************************************
#Lista modificable de los métodos que queremos utilizar para el ML
methods_list <- c('svmRadial', 'xgbTree', 'rf', 'lda')

methods_info <- list()

for (method in methods_list){
  info <- getModelInfo(model=method)[[method]]$parameters
  methods_info[[method]] <- info
}

#Establecemos el train control para todos los métodos
ctrl <- trainControl(
  method          = "cv",         # validación cruzada
  number          = 10,            # 10 folds
  classProbs      = TRUE,          # predecir probabilidades
  summaryFunction = multiClassSummary, # sensitivity, specificity and the area under the ROC curve
  verboseIter     = TRUE, #Mostramos información por pantalla de las iteraciones
  sampling        = 'smote'
)

#!!!!!!!!!!!!!
#Más adelante pondré un bucle para que haga todos los métodos con un bucle, pero de momento lo dejo de 1 en 1
#Es más fácil para hacer modificaciones así
#!!!!!!!!!!!!!

#===============================================================================
###---- SVM‐RBF----
#===============================================================================
#En este caso, sería necesario estandarizar los datos durante el preprocesado
tuneGrid_svm <- expand.grid(
  C = 2^(-5:7), #Penalización por error
  sigma = 2^(-15:-5) #Ancho del kernel
)

svm_mod_2 <- caret::train(
  x          = training_data_en[[2]],
  y          = training_labels[[2]],
  method     = "svmRadial",
  tuneGrid   = tuneGrid_svm,
  metric     = "Mean_F1", #funciona bien con categorías desbalanceadas
  trControl  = ctrl
)

print(svm_mod$results)
#1: Fitting sigma = 0.000122, C = 32 on full training set
#2: Fitting sigma = 0.000488, C = 8 on full training set
#3: 
#4: 

table(training_labels[[1]])

#Comprobamos si existe alguna columna consante en el modelo
constantes <- sapply(training_data_en[[1]], function(x) length(unique(x)) == 1)
names(training_data_en[[1]])[constantes]

#===============================================================================
###---- XGboost ----
#===============================================================================
#
tuneGrid_xgb <- expand.grid(
  nrounds = c(100, 200, 300), #número de boosting rounds
  max_depth = c(3, 6, 9), #máxima complejidad del árbol
  eta = c(0.01, 0.1, 0.3), #learning rate
  gamma = c(0, 1, 5), #penalización mínima para hacer split
  colsample_bytree =  c(0.5, 0.7, 1), #proporción de variables para cada árbol
  min_child_weight = c(1, 5, 10), #observaciones mínimas hijo para hacer split
  subsample = c(0.5, 0.75, 1) #Fracción de datos para cada árbol
)

xgb_mod_2 <- caret::train(
  x          = training_data_en[[2]],
  y          = training_labels[[2]],
  method     = "xgbTree",
  tuneGrid   = tuneGrid_xgb, #Modificamos estos valores
  metric     = "Mean_F1",
  trControl  = ctrl
)

#1. Fitting nrounds = 100, max_depth = 9, eta = 0.01, gamma = 0, colsample_bytree = 0.7, min_child_weight = 1, subsample = 1 on full training set
#2. Fitting nrounds = 100, max_depth = 9, eta = 0.1, gamma = 1, colsample_bytree = 1, min_child_weight = 10, subsample = 0.5 on full training set
#3. 
#4. 



#===============================================================================
###---- Random Forest----
#===============================================================================
#Este es el más recomendado para analizar datos de RNAseq
tuneGrid_rf <-  expand.grid(
  mtry = floor(seq(1, sqrt(length(training_labels[[2]])))) #La raíz cuadrada es estándar para procesos de clasificación
)
#mtry = floor(seq(1, sqrt(n_pc), length.out = 5))
#mtry = c(1, floor(sqrt(p)), floor(p/3))

rf_mod_2 <- caret::train(
  x          = training_data_en[[2]],
  y          = training_labels[[2]],
  method     = "rf",
  tuneGrid   = tuneGrid_rf,
  metric     = "Mean_F1",
  trControl  = ctrl,
  importance = TRUE
)

#1. Fitting mtry = 8 on full training set (con mensaje de warning)
#2. Fitting mtry = 4 on full training set (con mensaje de warning)
#3. 
#4.



#===============================================================================
###---- LDA/FDA ----
#===============================================================================

lda_mod_2 <- caret::train(
  x          = training_data_en[[2]],
  y          = training_labels[[2]],
  method    = "lda",
  trControl = trainControl(method = "cv")
)

library(earth)
fda_mod_2 <- caret::train(
  x          = training_data_en[[2]],
  y          = training_labels[[2]],
  method    = "fda",
  trControl = trainControl(method = "cv")
)

#Con el 2 sale un mensaje de warning

#****************************************************************************
## ---- Evaluación de métodos ----
#****************************************************************************

library(dplyr)
library(purrr)
library(tibble)
library(pROC)

# 1. Lista nombrada de los modelos entrenados
models_list <- list(
  SVM_RBF       = list(svm_mod_1, svm_mod_2, svm_mod_3, svm_mod_4),
  XGboost       = list(xgb_mod_1, xgb_mod_2, xgb_mod_3, xgb_mod_4),
  Random_Forest = list(rf_mod_1, rf_mod_2, rf_mod_3, rf_mod_4),
  LDA           = list(lda_mod_1, lda_mod_2, lda_mod_3, lda_mod_4),
  FDA           = list(fda_mod_1, fda_mod_2, fda_mod_3, fda_mod_4)
)

#Creamos una función para calcular las medias de los valores. 
calculos_metricas_binarias <- function(confusion_matrix) { 
  accuracy <- as.numeric(confusion_matrix$overall["Accuracy"])
  sensitivity <- as.numeric(mean(confusion_matrix$byClass[, "Sensitivity"]))
  specificity <- as.numeric(mean(confusion_matrix$byClass[, "Specificity"]))
  precision <- as.numeric(mean(confusion_matrix$byClass[, "Pos Pred Value"])) #La precisión es el pos pred value
  f1_score <- as.numeric(2*(precision*sensitivity)/(precision + sensitivity)) #Es necesario calcular la F1 score que
  list( #Creamos una lista con los resultados generados
    accuracy = accuracy, 
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1_score = f1_score
  )
}

#Creamos un dataframe vacío para guardar los resultados de la evaluación
metricas_evaluacion <- data.frame(
  Método = names(models_list), 
  Exactitud = numeric(length(models_list)), 
  Precisión = numeric(length(models_list)), 
  Sensibilidad = numeric(length(models_list)), 
  Especificidad = numeric(length(models_list)), 
  F1score = numeric(length(models_list)))

#Creamos dos listas vacías para guardar los resultados
conf_list <- list(
  train_gen_num_conf = list(),
  train_gen_enc_conf = list(),
  train_mut_num_conf = list(),
  train_mut_enc_conf = list()
  )

eval_list <- list(
  train_gen_num_eval = metricas_evaluacion,
  train_gen_enc_eval = metricas_evaluacion,
  train_mut_num_eval = metricas_evaluacion,
  train_mut_enc_eval = metricas_evaluacion
)


tablas_bonitas <- list()

titulos_tablas <- c(
  "**A) Métricas prueba DNAseq: número de mutaciones/gen**",
  "**B) Métricas prueba DNAseq: presencia de mutaciones/gen**",
  "**C) Métricas prueba DNAseq: número de variaciones/gen**",
  "**D) Métricas prueba DNAseq: presencia de variaciones/gen**"
)

for (j in 1:4) {
  test_data <- testing_data_en[[j]]
  test_label <- testing_labels[[j]]
  for (i in 1:length(models_list)){
    #Calculamos la preducción y la matriz de confusión
    pred <- predict(models_list[[i]][[j]], test_data)
    conf <- caret::confusionMatrix(pred, test_label)
    
    conf_list[[j]][[names(models_list)[i]]] <- conf
    
    resultado <- calculos_metricas_binarias(conf)
    
    #Añadimos los valores de las métricas al dataframe. 
    eval_list[[j]][i, "Exactitud"] <- resultado$accuracy 
    eval_list[[j]][i, "Precisión"] <- resultado$precision
    eval_list[[j]][i, "Sensibilidad"] <- resultado$sensitivity
    eval_list[[j]][i, "Especificidad"] <- resultado$specificity
    eval_list[[j]][i, "F1score"] <- resultado$f1_score
  }
  tabla_metricas <- eval_list[[j]] %>%
    gt() %>%
    cols_label(
      Método = "Método",
      Exactitud = "Exactitud", 
      Precisión = "Precisión", 
      Sensibilidad = "Sensibilidad", 
      Especificidad = "Especificidad", 
      Precisión = "Precisión")%>%
    fmt_number(
      decimals = 2,
      drop_trailing_zeros = FALSE) %>%
    tab_options( 
      table.font.size = 12,  #tamaño de fuente
      column_labels.background.color = "azure", #color de los títulos de las columnas
      column_labels.font.size = 14, #tamaño de fuente de títulos
      data_row.padding = px(10)) %>% #separación
    tab_style( #Modificamos los ajustes de las celdas del cuerpo de la tabla
      style = list(cell_text(align = "center"), #Centrar el texto
                   cell_borders(sides = "all", color = "lightblue", style = "solid", weight = px(2))), #Añadir bordes a las celdas
      locations = cells_body())%>%
    tab_style( #Modificamos los ajustes de las celdas de título de columnas
      style = list(cell_text(align = "center", weight = "bold"), 
                   cell_borders(sides = c("all"), color = "lightblue", style = "solid", weight = px(2))),
      locations = cells_column_labels()) %>%
    tab_style( #Modificamos los ajustes del título de la tabla
      style = list(cell_text(align = "left"),
                   cell_borders(sides = "bottom", color = "lightblue", style = "solid", weight = px(5))),
      locations = cells_title(groups = "title"))%>%
    tab_header(
      title = md(titulos_tablas[[j]]))
  tablas_bonitas[[j]] <- tabla_metricas
  gtsave(tabla_metricas, file.path(base_dir, paste0(names(eval_list)[j], "_tabla.png")))
}

rm(tablas_bonitas)
rm(tabla_metricas)
gc()

saveRDS(conf_list, file = file.path(base_dir, "conf_list.rds"))
saveRDS(eval_list, file = file.path(base_dir, "eval_list.rds"))

#****************************************************************************
## ---- Guardar los modelos ----
#****************************************************************************


saveRDS(svm_mod_2, file = file.path(base_dir,"svm_mod_mc_2.rds"))
saveRDS(xgb_mod_2, file = file.path(base_dir,"xgb_mod_mc_2.rds"))
saveRDS(rf_mod_2, file = file.path(base_dir,"rf_mod_mc_2.rds"))
saveRDS(lda_mod_2, file = file.path(base_dir,"lda_mod_mc_2.rds"))
saveRDS(fda_mod_2, file = file.path(base_dir,"fda_mod_mc_2.rds"))


#######################
### ***** FIN ***** ###
#######################
