#***************************************************************************************************************
#----Adquisición de datos----
#***************************************************************************************************************
#install.packages("BiocManager")
#BiocManager::install("recount3")
#BiocManager::install("maftools")
#BiocManager::install("SummarizedExperiment")
#BiocManager::install("sva")



library(recount3)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(data.table)
library(tibble)
library(sva)
library(caret)
#library(limma)

library(gt)
library(webshot2)

set.seed(123) #reproducibilidad

#****************************************************************************
## ---- Importar datos ----
#****************************************************************************
base_dir <- 'C:/Users/Ana/OneDrive - Universidad de Salamanca/Bioinformática/4. TFM/Datos'

#Datos RAW base de datos
rna_tcga <- fread(file.path(base_dir, "rna_tcga.csv"), na.strings = c("", " "))
rna_tcga <- column_to_rownames(rna_tcga, var = "Gen_Id")
rna_gtex <- fread(file.path(base_dir, "rna_gtex.csv"), na.strings = c("", " "))
rna_gtex <- column_to_rownames(rna_gtex, var = "Gen_Id")
info_df <- fread(file.path(base_dir, "info_df.csv"))




#Datos training/testing
training_data_rna <- readRDS(file.path(base_dir,"training_data_rna.rds"))
training_labels_rna <- readRDS(file.path(base_dir,"training_labels_rna.rds"))
testing_data_rna <- readRDS(file.path(base_dir,"testing_data_rna.rds"))
testing_labels_rna <- readRDS(file.path(base_dir,"testing_labels_rna.rds"))


#ML_models
svm_mod_1 <- readRDS(file.path(base_dir,"svm_mod_rna_1.rds"))
xgb_mod_2 <- readRDS(file.path(base_dir,"xgb_mod_rna_2.rds"))
rf_mod_1 <- readRDS(file.path(base_dir,"rf_mod_rna_1.rds"))
lda_mod_1 <- readRDS(file.path(base_dir,"lda_mod_rna_1.rds"))
fda_mod_1 <- readRDS(file.path(base_dir,"fda_mod_rna_1.rds"))

fda_mod_2 <- fda_mod
svm_mod_2 <- svm_mod
lda_mod_2 <- lda_mod
rf_mod_2 <- rf_mod

rm(fda_mod, lda_mod, svm_mod, rf_mod)

#****************************************************************************
## ---- Recount3 ----
#****************************************************************************
#setwd("C:/Users/Ana/Recount3")


proyectos <- available_projects()

###---- 1. RNAseq TCGA ----
rse_tcga <- recount3::create_rse_manual(
  project = "UCEC",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

data_type <- assayNames(rse_tcga) #raw counts

matrix_tcga <- recount3::transform_counts(rse_tcga, by='auc') #Normalizamos con AUC (permite comparaciones netre estudios)
#La otra opción es mapped_reads que no es tan robusto para comparar entre estudios


rna_tcga <- as.data.frame(matrix_tcga)
uuids <- c(colnames(rna_tcga))
rna_tcga <- setDT(rna_tcga, keep.rownames = 'Gen_Id')[]

#fwrite(rna_tcga, "C:\\Users\\anahe\\TCGA\\rna_tcga.csv")



#Guardamos los metadatos del rse
info <- colData(rse_tcga) #obtenemos los metadatos del archivo tipo rse

#Para filtrar por el tipo de muestra: info@listData[["tcga.cgc_sample_sample_type"]]

info_df <- data.frame( #obtenemos un dataframe de equivalencia de IDs de TCGA y de recount3.
  recount_id = info$external_id,
  tcga_id = info$tcga.tcga_barcode,
  tissue_type = info$tcga.cgc_sample_sample_type,
  stringsAsFactors = FALSE
)

#fwrite(info_df, "C:\\Users\\anahe\\TCGA\\info_df.csv")


###----2. RNAseq GTEx ----
rse_gtex <- recount3::create_rse_manual(
  project = "UTERUS",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

matrix_gtex <- recount3::transform_counts(rse_gtex)

rna_gtex <- as.data.frame(matrix_gtex)
rna_gtex <- setDT(rna_gtex, keep.rownames = 'Gen_Id')[]

#fwrite(rna_gtex, "C:\\Users\\anahe\\TCGA\\rna_gtex.csv")







#***************************************************************************************************************
#----Contol de calidad y selección----
#***************************************************************************************************************
###----1. Filtrado datos RNAseq enfermos----
#Pasamos la primera columna a nombres de fila (librería tibble). Sólo cuando importamos los datos, no cuando los descargamos
rna_tcga <- column_to_rownames(rna_tcga, var = "Gen_Id")

#Filtramos los datos de interés en el df de info_df
info_filtrado <- info_df %>%
  filter(tissue_type == "Primary Tumor") %>% #seleccionamos los tumores primarios
  mutate(tcga_id_short = substr(tcga_id, 1, 12)) %>% #acortamos el nombre de los pacientes
  distinct(recount_id, .keep_all = TRUE) #eliminamos las muestras repetidas

#Nos quedamos con los pacientes (columnas) seleccionados en el filtrado
rna_tcga_filtrado <- rna_tcga[, info_filtrado$recount_id, drop = FALSE]
colnames(rna_tcga_filtrado) <- info_filtrado$tcga_id_short

#Realizamos el filtrado
#Eliminamos las columnas (pacientes) donde la expresión es =0 en todos los genes
suma <- colSums(rna_tcga_filtrado)
cero_cols <- names(suma[suma == 0])
rna_tcga_filtrado <- rna_tcga_filtrado[ , !colnames(rna_tcga_filtrado) %in% cero_cols]

#Eliminamos las filas (genes) donde la expresión =0 en todos los pacientes
suma <- rowSums(rna_tcga)
cero_rows <- names(suma[suma == 0])
rna_tcga_filtrado <- rna_tcga_filtrado[!rownames(rna_tcga) %in% cero_rows, ]

#Eliminamos los genes que se expresan en menos de un 10% de los pacientes
p_ceros <- rowSums(rna_tcga_filtrado == 0) / ncol(rna_tcga_filtrado)
rna_tcga_filtrado <- rna_tcga_filtrado[p_ceros <= 0.9, ]

lista_genes <- rownames(rna_tcga_filtrado)


###----2. Filtrado datos RNAseq tejido sano (TCGA)----

#Seleccionamos a los pacientes sanos para utilizarlos como control para el batch effect

info_sano <- info_df %>%
  filter(tissue_type == "Solid Tissue Normal") %>% #seleccionamos el tejido sano
  mutate(tcga_id_short = substr(tcga_id, 1, 12)) %>% #acortamos el nombre de los pacientes
  distinct(recount_id, .keep_all = TRUE)

rna_tcga_sano <- rna_tcga[, info_sano$recount_id, drop = FALSE]
#Importante NO cambiarle el nombre de momento porque si no se repiten los TCGA IDs

#Realizamos el filtrado
#Eliminamos las columnas (pacientes) donde la expresión es =0 en todos los genes
suma <- colSums(rna_tcga_sano)
cero_cols <- names(suma[suma == 0])
rna_tcga_sano <- rna_tcga_sano[ , !colnames(rna_tcga_sano) %in% cero_cols]

#Selecionamos los mismos genes que en los enfermos
rna_tcga_sano <- rna_tcga_sano[rownames(rna_tcga_sano) %in% lista_genes, ]



###----3. Filtrado datos RNAseq pacientes sanos (GTEX)----

#Realizamos el mismo filtrado que para los datos de TCGA
#Eliminamos las columnas (pacientes) donde la expresión es =0 en todos los genes
suma <- colSums(rna_gtex)
cero_cols <- names(suma[suma == 0])
rna_gtex_filtrado <- rna_gtex[ , !colnames(rna_gtex) %in% cero_cols]

#Eliminamos las filas (genes) donde la expresión =0 en todos los pacientes
suma <- rowSums(rna_gtex)
cero_rows <- names(suma[suma == 0])
rna_gtex_filtrado <- rna_gtex_filtrado[!rownames(rna_gtex) %in% cero_rows, ]

#Eliminamos los genes que se expresan en menos de un 10% de los pacientes
p_ceros <- rowSums(rna_gtex_filtrado == 0) / ncol(rna_gtex_filtrado)
rna_gtex_filtrado <- rna_gtex_filtrado[p_ceros <= 0.9, ]

#Más adelante seleccionamos los mismos genes en ambos dataframes

##Exportamos los df
rna_tcga_filtrado_export <- setDT(rna_tcga_filtrado, keep.rownames = 'Gen_Id')[]
fwrite(rna_tcga_filtrado_export, file.path(base_dir, "rna_tcga_filtrado.csv"))
rna_tcga_sano_export <- setDT(rna_tcga_sano, keep.rownames = 'Gen_Id')[]
fwrite(rna_tcga_sano_export, file.path(base_dir, "rna_tcga_sano.csv"))
rna_gtex_filtrado_export <- setDT(rna_gtex_filtrado, keep.rownames = 'Gen_Id')[]
fwrite(rna_gtex_filtrado_export, file.path(base_dir, "rna_gtex_filtrado.csv"))
rm(rna_tcga_filtrado_export, rna_tcga_sano_export, rna_gtex_filtrado_export)

rna_tcga_filtrado <- column_to_rownames(rna_tcga_filtrado, var = "Gen_Id")
rna_tcga_sano <- column_to_rownames(rna_tcga_sano, var = "Gen_Id")
rna_gtex_filtrado <- column_to_rownames(rna_gtex_filtrado, var = "Gen_Id")

#***************************************************************************************************************
## ---- Preprocesado ----
#***************************************************************************************************************
###----1. Tratamiento de NAs----

#Comprobamos que no existen NAs tras el filtrado de los datos

any(is.na(rna_tcga_filtrado))
any(is.na(rna_tcga_sano))
any(is.na(rna_gtex_filtrado))


impute_na <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

if (any(is.na(rna_tcga_filtrado))){
  rna_tcga_filtrado <- rna_tcga_filtrado %>% 
    mutate(across(where(is.numeric), impute_na))
}

if (any(is.na(rna_tcga_sano))){
  rna_tcga_filtrado <- rna_tcga_sano %>% 
    mutate(across(where(is.numeric), impute_na))
}

if (any(is.na(rna_gtex_filtrado))){
  rna_gtex_filtrado <- rna_gtex_filtrado %>% 
    mutate(across(where(is.numeric), impute_na))
}



###----2. Tratamiento de datos desbalanceados----
#Seleccionamos una muestra aleatoria (mismo sample size que los sanos).

sample_size <- as.numeric(ncol(rna_gtex_filtrado))

patients_id_reduced <- data.frame(patient_id = sample(colnames(rna_tcga_filtrado), size = sample_size))

rna_tcga_reducido <- subset(rna_tcga_filtrado, select = patients_id_reduced$patient_id)



###----3. Unión de datos RNAseq----
#Filas a columnas
rna_tcga_t <- as.data.frame(t(rna_tcga_reducido))
rna_gtex_t <- as.data.frame(t(rna_gtex_filtrado))

gene_list <- intersect(colnames(rna_tcga_t), colnames(rna_gtex_t))
rna_tcga_t <- select(rna_tcga_t, all_of(gene_list))
rna_gtex_t <- select(rna_gtex_t, all_of(gene_list))

rna_seq <- bind_rows(Enfermo_Tejido_Tumoral = rna_tcga_t, Sano_Tejido_Sano = rna_gtex_t, .id = "Grupo")
#rna_seq_export <- copy(rna_seq)
#rna_seq_export <- setDT(rna_seq_export, keep.rownames = 'patient_id')
#fwrite(rna_seq_export, file.path(base_dir, "rna_seq.csv"))

#Seleccionamos esos mismos genes para los pacientes sanos de TCGA
rna_tcga_sano_t <- as.data.frame(t(rna_tcga_sano))
gene_list <- colnames(rna_seq[, -1])
rna_tcga_sano_t <- select(rna_tcga_sano_t, all_of(gene_list))
grupo_df <- data.frame(Grupo = rep("Enfermo_Tejido_Sano", times = nrow(rna_tcga_sano_t)))
rna_t_sano <- cbind(grupo_df, rna_tcga_sano_t)
#rna_t_sano_export <- copy(rna_t_sano)
#rna_t_sano_export <- setDT(rna_t_sano_export, keep.rownames = 'patient_id')[]
#fwrite(rna_t_sano_export, file.path(base_dir, "rna_t_sano.csv"))





###----4. División datos entrenamiento----
target_var <- "Grupo"
training_samples_rna <- rna_seq[[target_var]] %>% createDataPartition(p = 0.8, list = FALSE) #Estratificamos por grupo. 80/20

train_rna <- rna_seq[training_samples_rna,] #Enfermos y sanos seleccionados para el entrenamiento
train_rna_num <- train_rna[,-1]
train_rna_label <- as.factor(train_rna[[target_var]])

#Generamos otro set de entrenamiento con muestras de tejido sano de pacientes enfermos
#Las muestras de tejido sano únicamente van a encontrarse en el grupo de entrenamiento para el batch normalization
train_rna_ts <- rbind(train_rna, rna_t_sano) #ts de tejido sano
train_rna_ts_num <- train_rna_ts[,-1]
train_rna_ts_label <- as.factor(train_rna_ts[[target_var]])


test_rna <- rna_seq[-training_samples_rna,] #Datos seleccionados para la prueba
test_rna_num <- test_rna[,-1]
test_rna_label <- as.factor(test_rna[[target_var]])

###----5. Transformación y escalado----

####----5.1 Transformación logarítmica----

train_rna_num <- train_rna_num %>%
  mutate(across(everything(), ~ log(. + sqrt(.^2 + 1)))) #hyperbolic arcsine (asinh) transformation

train_rna_ts_num <- train_rna_ts_num %>%
  mutate(across(everything(), ~ log(. + sqrt(.^2 + 1)))) #hyperbolic arcsine (asinh) transformation

test_rna_num <- test_rna_num %>%
  mutate(across(everything(), ~ log(. + sqrt(.^2 + 1)))) #hyperbolic arcsine (asinh) transformation


####----5.2 Scaling ----
train_rna_num <- data.frame(scale(train_rna_num))
train_rna_scaled <- cbind(train_rna_label, train_rna_num)
colnames(train_rna_scaled)[1] <- 'Grupo'

train_rna_ts_num <- data.frame(scale(train_rna_ts_num))
train_rna_ts_scaled <- cbind(train_rna_ts_label, train_rna_ts_num)
colnames(train_rna_ts_scaled)[1] <- 'Grupo'

test_rna_num <- data.frame(scale(test_rna_num))
test_rna_scaled <- cbind(test_rna_label, test_rna_num)
colnames(test_rna_scaled)[1] <- 'Grupo'

####----5.3 Guardamos----
train_rna_scaled_export <- copy(train_rna_scaled)
train_rna_scaled_export <- setDT(train_rna_scaled_export, keep.rownames = 'patient_id')[]
fwrite(train_rna_scaled_export, file.path(base_dir, "train_rna_scaled_export.csv"))

train_rna_ts_scaled_export <- copy(train_rna_ts_scaled)
train_rna_ts_scaled_export <- setDT(train_rna_ts_scaled_export, keep.rownames = 'patient_id')[]
fwrite(train_rna_ts_scaled_export, file.path(base_dir, "train_rna_ts_scaled.csv"))

test_rna_scaled_export <- copy(test_rna_scaled)
test_rna_scaled_export <- setDT(test_rna_scaled_export, keep.rownames = 'patient_id')[]
fwrite(test_rna_scaled_export, file.path(base_dir, "test_rna_scaled_export.csv"))


#***************************************************************************************************************
##----Batch Effect----
#***************************************************************************************************************

###----1. Reducción de la dimensionalidad y representación del conjunto completo de los datos----
#Lo hacemos con los datos con el tejido sano de enfermos
pca_model_1 <- prcomp(train_rna_ts_num, center = FALSE, scale = FALSE) #Importante no usar las etiquetas
pca_df_1 <- data.frame(pca_model_1$x)

#Varianza explicada (varianza/total varianza)
varianza_explicada <- pca_model_1$sdev^2 / sum(pca_model_1$sdev^2)
varianza_acumulada <- cumsum(varianza_explicada)
n_pc_1 <- as.numeric(min(which(varianza_acumulada > 0.7)))

#Visualizamos dos dimensiones 
pca_plot_1 <- ggplot(pca_df_1, aes(x = PC1, y = PC2, color = train_rna_ts_label)) +  
  geom_point(size = 3) + 
  scale_color_manual(values = c("orange", "red", "green")) + 
  labs(title = "PCA: distribución antes de la correción del Batch Effect", 
       x = paste0(paste("PC1",round(varianza_explicada[1]*100, 2)),'%'), 
       y = paste0(paste("PC2",round(varianza_explicada[2]*100, 2)),'%'), 
       color = "Grupo") + 
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), 
        legend.position = 'bottom',
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "black"),
        #Cambiamos el tamaño de letra
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

pca_plot_1

ggsave("pca_plot_pre_BE.png", plot = pca_plot_1, width = 8, height = 6, dpi = 300)

###----2. Cuantificar----
#Calculamos las distancias de los centroides (dimensiones para conseguir al menos un 70% de la varianza):
centroides <- aggregate(. ~ train_rna_ts_label, data = pca_df_1[,1:n_pc_1], FUN = mean)
coords <- centroides[,-1]
rownames(coords) <- centroides$train_rna_ts_label
distancias <- as.matrix(dist(coords))
print(distancias)
#write.csv(distancias, file = file.path(base_dir, "distancias_centroides_preBE.csv"), row.names = TRUE)



###----3.Solucionar Batch Effect (sva::ComBat)----
#Necesitamos el conjunto completo de los datos con las correspondientes etiquetas.
#Transponemos porque necesitamos que cada paciente sea una columna
#En la matriz de expresión NO pueden aparecer las etiquetas
expr_matrix <- as.matrix(t(rbind(train_rna_ts_num, test_rna_num)))
labels <- c(train_rna_ts_label, test_rna_label)

#Lo que vamos a consguir en este procedimiento es separar entre condición y batch.
#De esta manera, como disponemos de datos sanos en la base de datos de enfermos, podemos equiparar los datos de sanos de la otra base de datos a estos.
#Si no tuviéramos muestras de tejido sano en pacientes enfermos no podríamos eliminar el batch effect. 

#Separamos las etiquetas en dos condiciones: enfermos y sanos. Si son Enfermo_Tejido_Tumoral son enfermos y si no, sanos
condition <- factor(ifelse(labels == "Enfermo_Tejido_Tumoral", "enfermo", "sano"))
#Separamos las etiquetas en dos batchs: TCGA y GTEX. Si son Sano_Tejido_Sano son de GTEX y si no, TCGA.
batch <- factor(ifelse(labels == "Sano_Tejido_Sano", "GTEX", "TCGA"))

#Creamos la matriz modelo para utilizarla en ComBat
mod <- model.matrix(~ condition) 

#Aplicamos el modelo de batch al conjunto de los datos
combat_expr <- ComBat(dat = expr_matrix, batch = batch, mod = mod)
rna_seq_nb_R <- data.frame(t(combat_expr))

#Dividimos de nuevo en datos de training y testing. 
#Como se mantinen en el mismo orden que al juntoarlos, con mantener el número de filas es suficiente para que sean iguales.
train_rna_nbe_R <- rna_seq_nb_R[1:nrow(train_rna_ts_num),]
test_rna_nbe_R <- rna_seq_nb_R[(nrow(train_rna_ts_num) + 1):nrow(rna_seq_nb_R),]




####----4. Representación distribución corregida----
pca_model_2 <- prcomp(train_rna_nbe_R, center = FALSE, scale = FALSE)
pca_df_2 <- data.frame(pca_model_2$x)

#Varianza explicada (varianza/total varianza)
varianza_explicada <- pca_model_2$sdev^2 / sum(pca_model_2$sdev^2)
varianza_acumulada <- cumsum(varianza_explicada)
n_pc_2 <- as.numeric(min(which(varianza_acumulada > 0.7)))

#Visualizamos dos dimensiones 
pca_plot_2 <- ggplot(pca_df_2, aes(x = PC1, y = PC2, color = train_rna_ts_label)) +  
  geom_point(size = 3) + 
  scale_color_manual(values = c("orange", "red", "green")) + 
  labs(title = "PCA: distribución después de la correción del Batch Effect (R)", 
       x = paste0(paste("PC1",round(varianza_explicada[1]*100, 2)),'%'), 
       y = paste0(paste("PC2",round(varianza_explicada[2]*100, 2)),'%'), 
       color = "Grupo") + 
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), 
        legend.position = 'bottom',
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "black"),
        #Cambiamos el tamaño de letra
        plot.title = element_text(size = 18, face = "bold", hjust = 1),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

pca_plot_2

ggsave("pca_plot_post_BE_R.png", plot = pca_plot_2, width = 8, height = 6, dpi = 300)

####----5. Cuantificar distancias de nuevo----
centroides_r <- aggregate(. ~ train_rna_ts_label, data = pca_df_2[,1:n_pc_2], FUN = mean)
coords_r <- centroides_r[,-1]
rownames(coords_r) <- centroides_r$train_rna_ts_label
distancias_r <- as.matrix(dist(coords_r))
print(as.matrix(distancias_r))

#write.csv(distancias_r, file = file.path(base_dir, "distancias_centroides_postBE_R.csv"), row.names = TRUE)



###----6.Importar datos Batch effect corregido (python)----
train_rna_nbe_python <- fread(file.path(base_dir, "RNA_train_batchCombat_corrected.csv"), na.strings = c("", " "))
train_rna_nbe_python <- column_to_rownames(train_rna_nbe_python, var = "patient_id")
train_rna_num_nbe_python <- train_rna_nbe_python[,-1]
train_rna_label_nbe_python <- as.factor(train_rna_nbe_python[[1]])


test_rna_nbe_python <- fread(file.path(base_dir, "RNA_test_batchCombat_corrected.csv"), na.strings = c("", " "))
test_rna_nbe_python <- column_to_rownames(test_rna_nbe_python, var = "patient_id")
test_rna_num_nbe_python <- test_rna_nbe_python[,-1]
test_rna_label_nbe_python <- as.factor(test_rna_nbe_python[[1]])


####----7. Representación distribución corregida----
pca_model_3 <- prcomp(train_rna_num_nbe_python, center = FALSE, scale = FALSE)
pca_df_3 <- data.frame(pca_model_3$x)

#Varianza explicada (varianza/total varianza)
varianza_explicada <- pca_model_3$sdev^2 / sum(pca_model_3$sdev^2)
varianza_acumulada <- cumsum(varianza_explicada)
n_pc_3 <- as.numeric(min(which(varianza_acumulada > 0.7)))

#Visualizamos dos dimensiones 
pca_plot_3 <- ggplot(pca_df_3, aes(x = PC1, y = PC2, color = train_rna_label_nbe_python)) +  #Las etiquetas que importo por si el orden de los pacientes es diferente
  geom_point(size = 3) + 
  scale_color_manual(values = c("orange", "red", "green")) + 
  labs(title = "PCA: distribución después de la correción del Batch Effect (Python)", 
       x = paste0(paste("PC1",round(varianza_explicada[1]*100, 2)),'%'), 
       y = paste0(paste("PC2",round(varianza_explicada[2]*100, 2)),'%'), 
       color = "Grupo") + 
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), 
        legend.position = 'bottom',
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "black"),
        #Cambiamos el tamaño de letra
        plot.title = element_text(size = 17, face = "bold", hjust = 1),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

pca_plot_3

ggsave("pca_plot_post_BE_Python.png", plot = pca_plot_3, width = 8, height = 6, dpi = 300)

####----8. Cuantificar distancias de nuevo----
centroides_p <- aggregate(. ~ train_rna_ts_label, data = pca_df_3[,1:n_pc_3], FUN = mean)
coords_p <- centroides_p[,-1]
rownames(coords_p) <- centroides_p$train_rna_ts_label
distancias_p <- as.matrix(dist(coords_p))
print(as.matrix(distancias_p))

#write.csv(distancias_p, file = file.path(base_dir, "distancias_centroides_postBE_python.csv"), row.names = TRUE)


####----9. Eliminar datos tejido sano pacientes enfermos----

#Datos obtenidos en R:
train_rna_nbe_R <- cbind(data.frame("Grupo" = train_rna_ts_label), train_rna_nbe_R)
#test_rna_nbe_R <- cbind(data.frame("Grupo" = test_rna_label), test_rna_nbe_R)
#Esto último es sólo para tenerlo unido a las etiquetas, que realmente no hace falta

train_rna_nbe_R <- train_rna_nbe_R %>%
  dplyr::filter(Grupo != "Enfermo_Tejido_Sano")

target_var = 'Grupo'
train_rna_nbe_num_R <- train_rna_nbe_R[,-1]
train_rna_nbe_label_R <- as.factor(train_rna_nbe_R[[target_var]])

#train_rna_nbe_R_export <- copy(train_rna_nbe_R)
#train_rna_nbe_R_export <- setDT(train_rna_nbe_R_export, keep.rownames = 'patient_id')[]
#fwrite(train_rna_nbe_R_export, file.path(base_dir, "train_rna_nbe_R.csv"))




#Datos obtenidos en Python:

train_rna_nbe_python <- train_rna_nbe_python %>%
  dplyr::filter(Grupo != "Enfermo_Tejido_Sano")

target_var = 'Grupo'
train_rna_nbe_num_python <- train_rna_nbe_python[,-1]
train_rna_nbe_label_python <- as.factor(train_rna_nbe_python[[target_var]])

#train_rna_nbe_python_export <- copy(train_rna_nbe_python)
#train_rna_nbe_python_export <- setDT(train_rna_nbe_python_export, keep.rownames = 'patient_id')[]
#fwrite(train_rna_nbe_python_export, file.path(base_dir, "train_rna_nbe_python.csv"))

#***************************************************************************************************************
##----Reducción de la dimensionalidad - PCA ----
#***************************************************************************************************************

###----1. Datos a los que no se les ha corregido el batch effect----
pca_model <- prcomp(train_rna_num, center = FALSE, scale = FALSE) #Construimos el modelo de PCA con los datos de entrenamiento
train_rna_pca <- as.data.frame(predict(pca_model, train_rna_num))
test_rna_pca <- as.data.frame(predict(pca_model, test_rna_num))

pca_genes <- pca_model$rotation #Esto nos vendrá bien para la interpretación posterior
genes <- data.frame(genes = rownames(pca_genes))


#Varianza explicada (varianza/total varianza)
varianza_explicada <- pca_model$sdev^2 / sum(pca_model$sdev^2)
varianza_acumulada_df <- data.frame( #Creamos una tabla con la varianza acumulada
  Component = 1:length(varianza_explicada),
  Variance = cumsum(varianza_explicada)
)
varianza_acumulada <- cumsum(varianza_explicada)
# Tomamos el numero de componentes principales que explican el 70% de la varianza
n_pc <- as.numeric(min(which(varianza_acumulada > 0.7)))

#Nos quedamos únicamente con el número de componentes principales establecido
train_rna_pca <- train_rna_pca[,1:n_pc]
test_rna_pca <- test_rna_pca[,1:n_pc]


#saveRDS(pca_model, file = file.path(base_dir,"pca_model.rds"))
#train_rna_pca_export <- copy(train_rna_pca)
#train_rna_pca_export <- setDT(train_rna_pca_export, keep.rownames = 'patient_id')[]
#fwrite(train_rna_pca_export, file.path(base_dir, "train_rna_pca.csv"))
#test_rna_pca_export <- copy(test_rna_pca)
#test_rna_pca_export <- setDT(test_rna_pca_export, keep.rownames = 'patient_id')[]
#fwrite(test_rna_pca_export, file.path(base_dir, "test_rna_pca.csv"))

#Obtenemos lo que influye cada gen a los componentes
genes_mean <- data.frame(
  Gene = rownames(pca_genes),
  MeanAbs = apply(pca_genes[,1:n_pc], 1, function(x) mean(abs(x)))
)

genes_mean <- genes_mean[order(-genes_mean$MeanAbs), ]
rownames(genes_mean) <- NULL
genes_relevantes <- genes_mean[1:5,]

###----2. Datos con Batch effect corregido (R)----
pca_model_be_R <- prcomp(train_rna_nbe_num_R, center = FALSE, scale = FALSE) #Construimos el modelo de PCA con los datos de entrenamiento

train_rna_nbe_pca_R <- as.data.frame(predict(pca_model_be_R, train_rna_nbe_num_R))
test_rna_nbe_pca_R <- as.data.frame(predict(pca_model_be_R, test_rna_nbe_R))

pca_genes_R <- pca_model_be_R$rotation #Esto nos vendrá bien para la interpretación posterior
genes_R <- data.frame(genes = rownames(pca_genes_R))

#Varianza explicada (varianza/total varianza)
varianza_explicada_R <- pca_model_be_R$sdev^2 / sum(pca_model_be_R$sdev^2)
varianza_acumulada_df_R <- data.frame( #Creamos una tabla con la varianza acumulada
  Component = 1:length(varianza_explicada_R),
  Variance = cumsum(varianza_explicada_R)
)
varianza_acumulada_R <- cumsum(varianza_explicada_R)
# Tomamos el numero de componentes principales que explican el 70% de la varianza
n_pc_nbe_R <- as.numeric(min(which(varianza_acumulada > 0.7))) #73

#Nos quedamos únicamente con el número de componentes principales establecido
train_rna_nbe_pca_R <- train_rna_nbe_pca_R[,1:n_pc_nbe_R]
test_rna_nbe_pca_R <- test_rna_nbe_pca_R[,1:n_pc_nbe_R]


saveRDS(pca_model_be_R, file = file.path(base_dir,"pca_model_be_R.rds"))

train_rna_nbe_pca_R_export <- copy(train_rna_nbe_pca_R)
train_rna_nbe_pca_R_export <- setDT(train_rna_nbe_pca_R_export, keep.rownames = 'patient_id')[]
fwrite(train_rna_nbe_pca_R_export, file.path(base_dir, "train_rna_nbe_pca_R.csv"))

test_rna_nbe_pca_R_export <- copy(test_rna_nbe_pca_R)
test_rna_nbe_pca_R_export <- setDT(test_rna_nbe_pca_R_export, keep.rownames = 'patient_id')[]
fwrite(test_rna_nbe_pca_R_export, file.path(base_dir, "test_rna_nbe_pca_R.csv"))


###----3. Datos con Batch effect corregido (Python)----
pca_model_be_python <- prcomp(train_rna_nbe_num_python, center = FALSE, scale = FALSE) #Construimos el modelo de PCA con los datos de entrenamiento

train_rna_nbe_pca_python <- as.data.frame(predict(pca_model_be_python, train_rna_nbe_num_python))
test_rna_nbe_pca_python <- as.data.frame(predict(pca_model_be_python, test_rna_nbe_python))

pca_genes_python <- pca_model_be_python$rotation #Esto nos vendrá bien para la interpretación posterior
genes_python <- data.frame(genes = rownames(pca_genes_python))

#Varianza explicada (varianza/total varianza)
varianza_explicada_python <- pca_model_be_python$sdev^2 / sum(pca_model_be_python$sdev^2)
varianza_acumulada_df_python <- data.frame( #Creamos una tabla con la varianza acumulada
  Component = 1:length(varianza_explicada_python),
  Variance = cumsum(varianza_explicada_python)
)
varianza_acumulada_python <- cumsum(varianza_explicada_python)
# Tomamos el numero de componentes principales que explican el 70% de la varianza
n_pc_nbe_python <- as.numeric(min(which(varianza_acumulada_python > 0.7))) 

#Nos quedamos únicamente con el número de componentes principales establecido
train_rna_nbe_pca_python <- train_rna_nbe_pca_python[,1:n_pc_nbe_python]
test_rna_nbe_pca_python <- test_rna_nbe_pca_python[,1:n_pc_nbe_python]


saveRDS(pca_model_be_python, file = file.path(base_dir,"pca_model_be_python.rds"))

train_rna_nbe_pca_python_export <- copy(train_rna_nbe_pca_python)
train_rna_nbe_pca_python_export <- setDT(train_rna_nbe_pca_python_export, keep.rownames = 'patient_id')[]
fwrite(train_rna_nbe_pca_python_export, file.path(base_dir, "train_rna_nbe_pca_python.csv"))

test_rna_nbe_pca_python_export <- copy(test_rna_nbe_pca_python)
test_rna_nbe_pca_python_export <- setDT(test_rna_nbe_pca_python_export, keep.rownames = 'patient_id')[]
fwrite(test_rna_nbe_pca_python_export, file.path(base_dir, "test_rna_nbe_pca_python.csv"))


###----4. Guardamos en dos listas ----
training_data_rna <- list(train_rna_pca, train_rna_nbe_pca_R, train_rna_nbe_pca_python)
training_labels_rna <- list(train_rna_label, train_rna_label, train_rna_label)

testing_data_rna <- list(test_rna_pca, test_rna_nbe_pca_R, test_rna_nbe_pca_python)
testing_labels_rna <- list(test_rna_label, test_rna_label, test_rna_label)

saveRDS(training_data_rna, file = file.path(base_dir, "training_data_rna.rds"))
saveRDS(training_labels_rna, file = file.path(base_dir, "training_labels_rna.rds"))
saveRDS(testing_data_rna, file = file.path(base_dir, "testing_data_rna.rds"))
saveRDS(testing_labels_rna, file = file.path(base_dir, "testing_labels_rna.rds"))


#***************************************************************************************************************
## ---- Machine Learning ----
#***************************************************************************************************************

#Utilizamos la libreía caret para todos los métodos para "estandarizar el proceso"

#Vemos los parámetros que se pueden personalizar

methods_list <- c('svmRadial', 'glmnet', 'xgbTree', 'rf', 'lda')

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
  summaryFunction = twoClassSummary, # sensitivity, specificity and the area under the ROC curve
  verboseIter     = TRUE #Mostramos información por pantalla de las iteraciones
)


#===============================================================================
###---- SVM‐RBF----
#===============================================================================
#En este caso, sería necesario estandarizar los datos durante el preprocesado
tuneGrid_svm <- expand.grid(
  C = 2^(-5:7), #Penalización por error
  sigma = 2^(-15:-5) #Ancho del kernel
)

#En verdad se quedan en el mínimo así que debería ampliar el rango (puedo hacerlo porque no tarda mucho)

svm_mod <- caret::train(
  x          = training_data_rna[[3]],
  y          = training_labels_rna[[3]],
  method     = "svmRadial",
  tuneGrid   = tuneGrid_svm,
  metric     = "ROC",
  trControl  = ctrl
)

#1. (Lo hizo Martín)
#2. Fitting sigma = 0.0312, C = 0.0312 on full training set
#3. Fitting sigma = 0.0312, C = 2 on full training set

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

xgb_mod <- caret::train(
  x          = training_data_rna[[3]],
  y          = training_labels_rna[[3]],
  method     = "xgbTree",
  tuneGrid   = tuneGrid_xgb, #Modificamos estos valores
  metric     = "ROC",
  trControl  = ctrl
)


#1. (Martín)
#2. (Martín)
#3. Fitting nrounds = 100, max_depth = 6, eta = 0.1, gamma = 0, colsample_bytree = 0.7, min_child_weight = 1, subsample = 0.75 on full training set


#===============================================================================
###---- Random Forest----
#===============================================================================
#Este es el más recomendado para analizar datos de RNAseq
tuneGrid_rf <-  expand.grid(
  mtry = floor(seq(1, sqrt(length(training_labels_rna[[2]])))) #La raíz cuadrada es estándar para procesos de clasificación
)
#mtry = floor(seq(1, sqrt(n_pc), length.out = 5))
#mtry = c(1, floor(sqrt(p)), floor(p/3))

rf_mod <- caret::train(
  x          = training_data_rna[[3]],
  y          = training_labels_rna[[3]],
  method     = "rf",
  tuneGrid   = tuneGrid_rf,
  metric     = "ROC",
  trControl  = ctrl,
  importance = TRUE
)

#1.
#2.Fitting mtry = 1 on full training set
#3.Fitting mtry = 15 on full training set



#===============================================================================
###---- LDA/FDA ----
#===============================================================================

lda_mod <- caret::train(
  x          = training_data_rna[[3]],
  y          = training_labels_rna[[3]],
  method    = "lda",
  trControl = trainControl(method = "cv")
)
library('earth')

fda_mod <- caret::train(
  x          = training_data_rna[[3]],
  y          = training_labels_rna[[3]],
  method    = "fda",
  trControl = trainControl(method = "cv")
)



#****************************************************************************
## ---- Guardar los modelos ----
#****************************************************************************

saveRDS(svm_mod, file = file.path(base_dir,"svm_mod_rna_3.rds"))
saveRDS(xgb_mod, file = file.path(base_dir,"xgb_mod_rna_3.rds"))
saveRDS(rf_mod, file = file.path(base_dir,"rf_mod_rna_3.rds"))
saveRDS(lda_mod, file = file.path(base_dir,"lda_mod_rna_3.rds"))
saveRDS(fda_mod, file = file.path(base_dir,"fda_mod_rna_3.rds"))

#Para importar los modelos mirar el código arriba del todo del script


#***************************************************************************************************************
## ---- Validación con nuevos datos (base de datos diferente) ----
#***************************************************************************************************************


#****************************************************************************
#----Adquisición de datos----
#****************************************************************************

#BiocManager::install("GEOquery")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("biomaRt")


library(GEOquery)
library(biomaRt)
library(dplyr)
library(tibble)
library(data.table)

set.seed(123) #reproducibilidad

#****************************************************************************
## ---- Importar datos ----
#****************************************************************************
#Datos que necesito:
rna_seq <- fread(file.path(base_dir, "rna_seq.csv"), na.strings = c("", " "))
rna_seq <- column_to_rownames(rna_seq, var = "patient_id")

#Datos que genero:
rna_val <- fread(file.path(base_dir, "rna_val.csv"), na.strings = c("", " "))
rna_val <- column_to_rownames(rna_val, var = "patient_id")

rna_val_filtered <- fread(file.path(base_dir, "rna_val_filtered.csv"), na.strings = c("", " "))
rna_val_filtered <- column_to_rownames(rna_val_filtered, var = "patient_id")

#****************************************************************************
## ---- Descargar datos ----
#****************************************************************************

gds858 <- getGEO('GDS4589', destdir=base_dir) #Con este no nos descargamos el full

#Me lo descargué de la web y lo importé con la API
gds4589 <- getGEO(filename='C:/Users/Ana/OneDrive - Universidad de Salamanca/Bioinformática/4. TFM/Datos/rna_seq/raw/GDS4589_full.soft')

#Visualizamos las diferentes partes del archivo
str(gds4589)

#Extraemos la matriz de expresión y otros metadatos
exprs_matrix <- gds4589@dataTable@table
samples <- gds4589@dataTable@columns
metadata <- gds4589@header

gene_names <- exprs_matrix$IDENTIFIER

#****************************************************************************
#----Selección de datos Validación----
#****************************************************************************
###----1. Seleccionar genes únicos----
#Para los genes que están repetidos, seleccionamos los que presentan expresión máxima

unique_exprs_matrix <- exprs_matrix %>%
  group_by(IDENTIFIER) %>%
  summarise(across(everything(), max)) %>%
  filter(IDENTIFIER != ID_REF) %>%
  filter(IDENTIFIER != "--Control")

unique_gene_names <- unique_exprs_matrix$IDENTIFIER

###----2. Cambiar el nombre de los genes----

#Tengo que buscar y entender la documentación de la API. 
# Conectar al biomart de Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Consultar la tabla
mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id_version"),
  filters = "external_gene_name",
  values = unique_gene_names,
  mart = mart
)


colnames(mapping) <- c('IDENTIFIER', 'ENSEMBL_ID')
#Aquí vuelven a aparecer genes repetidos (para los que encuentra más de 1 ID) y hay algunos que desaparecen.


#Continuamos con la selección de pacientes (más tarde continúo con la selección de genes)

###----2. Selección de pacientes----
#Hago un inner join de las tablas y elimino las columnas que no me interesan. 
full_df <- merge(mapping, unique_exprs_matrix, by ="IDENTIFIER", all = FALSE)
full_df <- subset(full_df, select = -c(IDENTIFIER, ID_REF))
full_df <- full_df[1:104]


#Cambio las filas por columnas
full_t_df <- as.data.frame(t(full_df))
colnames(full_t_df) <- full_t_df[1, ]
full_t_df <- full_t_df[-1, ]
full_t_df <- full_t_df %>% rownames_to_column(var = "sample")




#Seleccionamos los pacientes que nos interesan
classification <- samples[c('sample', 'tissue')]

#Juntamos la clasificación con los datos de expresión
exp_type <- merge(classification, full_t_df, by = 'sample')

rna_val <- exp_type %>%
  filter(tissue != "papillary serous tumor") %>%
  mutate(tissue = case_when(
    tissue == "endometrioid tumor" ~ "Enfermo_Tejido_Tumoral",
    tissue == "inactive endometrium" ~ "Sano_Tejido_Sano",
  ))

rna_val <- column_to_rownames(rna_val, var = "sample")
colnames(rna_val)[1] <- "Grupo"

#Guardamos el dataframe en un csv
rna_val_export <- copy(rna_val)
rna_val_export <- setDT(rna_val_export, keep.rownames = 'patient_id')[]
fwrite(rna_val_export, file.path(base_dir, "rna_val.csv"))







###----3. Selección de los mismos genes usados en el entrenamiento----

train_genes <- colnames(rna_seq[2:ncol(rna_seq)])
train_genes_clean <- sub("\\..*", "", train_genes) #Eliminamos la versión después de los puntos para que coincidan más

train_df <- data.frame(
  train_genes = train_genes,
  genes_clean = train_genes_clean
)


val_df <- mapping %>%
  mutate(genes_clean = sub("\\..*", "", ENSEMBL_ID))
colnames(val_df)[2] <- "val_genes"

genes_df <- merge(train_df, val_df, by = "genes_clean")

names(table(genes_df$IDENTIFIER)[table(genes_df$IDENTIFIER) > 1])
#Vemos que algunos están repetidos pero no pasa nada porque vamos a seleccionarlos sólo una vez

rna_val_filtered <- dplyr::select(rna_val, all_of(genes_df$val_genes))
#Una vez seleccionados, cambiamos los nombres por los mismos del training. 
#Es un poco una guarrada pero es lo único que se me ocurre para no tener que repetir el entrenamiento con los nombres sin el .XX
colnames(rna_val_filtered) <- genes_df$train_genes
#Como están en el mismo orden, lo cambiamos por la notación de los de training para que coincida

rna_val_filtered_export <- copy(rna_val_filtered)
rna_val_filtered_export <- setDT(rna_val_filtered_export, keep.rownames = 'patient_id')[]
fwrite(rna_val_filtered_export, file.path(base_dir, "rna_val_filtered.csv"))

#****************************************************************************
## ---- Preprocesado ----
#****************************************************************************
#Después de seleccionar los genes hacemos las transformaciones (escalar). 
#Las transformaciones tienen que ser iguales que en los datos de entrenamiento
#Estos datos están normalizados en escala lineal
#https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2013.00139/full
#Como en el archivo de metadatos pone que vlaue type es count no tiene ningún tipo de normalización logarítmica.
#Lo lineal concuerda con los valores. No son raw counts porque superarían los miles de lecturas y aquí superan los 100 pero no llegan a tanto

###----1. Transformación y escalado----

####----1.1 Transformación logarítmica----
#Ya nos hemos quedado antes únicamente con las columnas numéricas (son las de los genes seleccionados)

rna_val_filtered[] <- lapply(rna_val_filtered, as.numeric)

rna_val_transformed <- rna_val_filtered %>%
  mutate(across(everything(), ~ log(. + sqrt(.^2 + 1)))) #hyperbolic arcsine (asinh) transformation

####----1.2 Scaling ----
rna_val_num <- data.frame(scale(rna_val_transformed))
rna_val_label <- rna_val[1]
rna_val_scaled <- cbind(rna_val_label, rna_val_num)

###----2. Añadimos columnas faltantes----
#Para el PCA hacen falta exactamente las mismas columnas que en los datos usados para el entrenamiento
#En nuestro caso tenemos menos así que las rellenamos con 0
#Si eliminaros genes del entrenamiento sería contraproducente

missing_cols <- setdiff(colnames(rna_seq[2:ncol(rna_seq)]), colnames(rna_val_scaled))
missing_cols_2 <- setdiff(mapping$ENSEMBL_ID, genes_df$val_genes)

rna_val_scaled[missing_cols] <- 0

# reordenar para que el orden sea idéntico al de los datos de entrenamiento
rna_val_scaled <- rna_val_scaled[, train_genes]


rna_val_scaled_export <- copy(rna_val_scaled)
rna_val_scaled_export <- setDT(rna_val_scaled_export, keep.rownames = 'patient_id')[]
fwrite(rna_val_scaled_export, file.path(base_dir, "rna_val_scaled.csv"))

#****************************************************************************
##----Reducción de la dimensionalidad----
#****************************************************************************

###----1. Datos sin eliminar Batch effect----
rna_val_pca <- as.data.frame(predict(pca_model, rna_val_scaled))
#Nos quedamos con el número de componentes establecido en la fase anterior
rna_val_pca <- rna_val_pca[,1:n_pc]

rna_val_pca_export <- copy(rna_val_pca)
rna_val_pca_export <- setDT(rna_val_pca_export, keep.rownames = 'patient_id')[]
fwrite(rna_val_pca_export, file.path(base_dir, "rna_val_pca.csv"))

###----2. Datos eliminando Batch effect (R)----
rna_val_pca_nbe_R <- as.data.frame(predict(pca_model_be_R, rna_val_scaled))
#Nos quedamos con el número de componentes establecido en la fase anterior
rna_val_pca_nbe_R <- rna_val_pca_nbe_R[,1:n_pc_nbe_R]

#rna_val_pca_nbe_R_export <- copy(rna_val_pca_nbe_R)
#rna_val_pca_nbe_R_export <- setDT(rna_val_pca_nbe_R_export, keep.rownames = 'patient_id')[]
#fwrite(rna_val_pca_nbe_R_export, file.path(base_dir, "rna_val_pca_nbe_R.csv"))

###----3. Datos eliminando Batch effect (Python)----
rna_val_pca_nbe_python <- as.data.frame(predict(pca_model_be_python, rna_val_scaled))
#Nos quedamos con el número de componentes establecido en la fase anterior
rna_val_pca_nbe_python <- rna_val_pca_nbe_python[,1:n_pc_nbe_python]

#rna_val_pca_nbe_python_export <- copy(rna_val_pca_nbe_python)
#rna_val_pca_nbe_python_export <- setDT(rna_val_pca_nbe_python_export, keep.rownames = 'patient_id')[]
#fwrite(rna_val_pca_nbe_python_export, file.path(base_dir, "rna_val_pca_nbe_python.csv"))

###----4. Guardar los datos en una lista----
target_var <- "Grupo"
val_labels <- as.factor(rna_val[[target_var]])
validating_data_rna <- list(rna_val_pca, rna_val_pca_nbe_R, rna_val_pca_nbe_python)
validating_labels_rna <- list(val_labels, val_labels, val_labels)











#***************************************************************************************************************
## ---- Machine Learning - Testing y Validación----
#***************************************************************************************************************

new_data_rna <- list(testing_data_rna, validating_data_rna)
new_labels_rna <- list(testing_labels_rna, validating_labels_rna)

###----Matrices de confusión y métricas----

library(dplyr)
library(purrr)
library(tibble)
library(pROC)

# 1. Lista nombrada de tus modelos entrenados
models_list <- list(
  SVM_RBF       = list(svm_mod_1, svm_mod_2, svm_mod_3),
  XGboost       = list(xgb_mod_1, xgb_mod_2, xgb_mod_3),
  Random_Forest = list(rf_mod_1, rf_mod_2, rf_mod_3),
  LDA           = list(lda_mod_1, lda_mod_2, lda_mod_3),
  FDA           = list(fda_mod_1, fda_mod_2, fda_mod_3)
)


#Creamos una función para calcular las medias de los valores. 
calculos_metricas_binarias_rna <- function(confusion_matrix) { #Creamos una función para calcular las medias de los valores. 
  accuracy <- as.numeric(confusion_matrix$overall[["Accuracy"]])
  sensitivity <- as.numeric(confusion_matrix$byClass[["Sensitivity"]])
  specificity <- as.numeric(confusion_matrix$byClass[["Specificity"]])
  precision <- as.numeric(confusion_matrix$byClass[["Pos Pred Value"]]) #La precisión es el pos pred value
  f1_score <- as.numeric(2*(precision*sensitivity)/(precision + sensitivity)) #Es necesario calcular la F1 score que depende de la prec y sensi
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






#Creamos dos listas vacías para guardar los resultados (de testing y de evaluación)
conf_list <- list(
  confusion <- list(
    train_rna_conf = list(),
    train_rna_nbe_R_conf = list(),
    train_rna_nbe_python_conf = list()),
  confusion_eval <- list(
    rna_eval_conf = list(),
    rna_eval_nbe_R_conf = list(),
    rna_eval_nbe_python_conf = list())
)

metrics_list <- list(
  metrics <- list(
    train_rna_met = metricas_evaluacion,
    train_rna_nbe_R_met = metricas_evaluacion,
    train_rna_nbe_python_met = metricas_evaluacion),
  metrics_eval <- list(
    rna_eval_met = metricas_evaluacion,
    rna_eval_nbe_R_met = metricas_evaluacion,
    rna_eval_nbe_python_met = metricas_evaluacion)
)


tablas_bonitas <- list(
  tablas_bonitas_testing <- list(),
  tablas_bonitas_validating <- list()
)


titulos_tablas <- list(
  titulos_tablas_testing <- list(
    "**A) Métricas prueba: RNAseq**",
    "**B) Métricas prueba: RNAseq sin batch effect (R)**",
    "**C) Métricas prueba: RNAseq sin batch effect (Python)**"),
  titulos_tablas_validating <- list(
    "**A) Métricas evaluación: RNAseq**",
    "**B) Métricas evaluación: RNAseq sin batch effect (R)**",
    "**C) Métricas evaluación: RNAseq sin batch effect (Python)**")
)


nombres_tablas <- list(
  nombres_tablas_testing <- list("testing","testing_nbe_R", "testing_nbe_python"),
  nombres_tablas_validating <- list("val", "val_nbe_R", "val_nbe_python")
)


for (k in 1:2){ #Testing y validating
  for (j in 1:3) { #normal/batch con R/batch con python
    new_data <- new_data_rna[[k]][[j]]
    new_label <- new_labels_rna[[k]][[j]]
    for (i in 1:length(models_list)){#Los diferentes métodos
      #Calculamos la preducción y la matriz de confusión
      pred <- predict(models_list[[i]][[j]], new_data)
      
      
      #print(k)
      #print(j)
      #print(i)
      #print(length(pred))
      #print(length(new_label))
      
      conf <- caret::confusionMatrix(pred, new_label)
      
      conf_list[[k]][[j]][[names(models_list)[i]]] <- conf
      
      resultado <- calculos_metricas_binarias_rna(conf)
      
      #Añadimos los valores de las métricas al dataframe. 
      metrics_list[[k]][[j]][i, "Exactitud"] <- resultado$accuracy 
      metrics_list[[k]][[j]][i, "Precisión"] <- resultado$precision
      metrics_list[[k]][[j]][i, "Sensibilidad"] <- resultado$sensitivity
      metrics_list[[k]][[j]][i, "Especificidad"] <- resultado$specificity
      metrics_list[[k]][[j]][i, "F1score"] <- resultado$f1_score
    }
    tabla_metricas <- metrics_list[[k]][[j]] %>%
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
        table.font.size = 13,  #tamaño de fuente
        column_labels.background.color = "azure", #color de los títulos de las columnas
        column_labels.font.size = 14, #tamaño de fuente de títulos
        data_row.padding = px(5)) %>% #separación
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
        title = md(titulos_tablas[[k]][[j]]))
    tablas_bonitas[[j]] <- tabla_metricas
    gtsave(tabla_metricas, file.path(base_dir, paste0(nombres_tablas[[k]][[j]], "_tabla.png")))
  }
}


saveRDS(conf_list, file = file.path(base_dir, "conf_list_rna.rds"))
saveRDS(metrics_list, file = file.path(base_dir, "metrics_list_rna.rds"))




###----Curvas ROC----

roc_list <- list(
  roc_test <- list(
    roc_rna <- list(),
    roc_rna_nbe_R <- list(),
    roc_rna_nbe_python <- list()),
  roc_eval <- list(
    roc_rna_val <- list(),
    roc_rna_nbe_R_val <- list(),
    roc_rna_nbe_python_val <- list())
)



auc_list <- list(
  auc_test <- list(
    auc_rna <- list(),
    auc_rna_nbe_R <- list(),
    auc_rna_nbe_python <- list()),
  aunc_eval <- list(
    auc_rna_eval <- list(),
    auc_rna_nbe_R_eval <- list(),
    auc_rna_nbe_python_eval <- list())
)


ROC_names <- list(
  testing <- list(
    'curvas_ROC_rna_testing.png', 
    'curvas_ROC_rna_nbe_R_testing.png',
    'curvas_ROC_rna_nbe_python_testing.png'),
  eval <- list(
    'curvas_ROC_rna_validating.png', 
    'curvas_ROC_rna_nbe_R_validating.png',
    'curvas_ROC_rna_nbe_python_validating.png')
)



ROC_titles <- list(
  testing <- list(
    'Curvas ROC prueba: RNAseq', 
    'Curvas ROC prueba: RNAseq sin Batch Effect (R)',
    'Curvas ROC prueba: RNAseq sin Batch Effect (Python)'),
  eval <- list(
    'Curvas ROC evaluación: RNAseq', 
    'Curvas ROC evaluación: RNAseq sin Batch Effect (R)',
    'Curvas ROC evaluación: RNAseq sin Batch Effect (Python)')
)



# Obtener etiquetas de clase
#Como son las mismas clases en testing y evaluación, conque lo hagamos una vez es suficiente
clases <- levels(test_rna_label)


for (k in 1:2) { #k es si es testing o validación
  for (j in 1:3) { #j es si es con batch o sin él
    new_data <- new_data_rna[[k]][[j]]
    new_label <- new_labels_rna[[k]][[j]]
    for (i in 1:length(models_list)) { #i es el modelo exacto de machine learning
      
      modelo <- models_list[[i]][[j]]
      
      probs <- predict(modelo, new_data, type = "prob")[, 1] #La clase positiva es la primera
      
      roc_obj <- roc(new_label, probs, quiet = TRUE)
      auc_obj <- auc(roc_obj)
      
      roc_list[[k]][[j]][[names(models_list)[i]]] <- roc_obj
      auc_list[[k]][[j]][[names(models_list)[i]]] <- auc_obj
    }
    
    #Hacemos la gráfica
    
    png(ROC_names[[k]][[j]], width = 800, height = 600) #Para guardar la gráfica (carpeta Documentos, no en basedir)
    
    colores <- c('purple', 'lightblue', 'pink', 'green', 'orange')[1:length(roc_list[[k]][[j]])]
    
    for (i in 1:length(roc_list[[k]][[j]])){
      if (i == 1){
        plot(roc_list[[k]][[j]][[i]], 
             col = colores[i], 
             main = ROC_titles[[k]][[j]], 
             lwd = 4,
             cex.main = 1.8,
             cex.axis = 1.4,
             cex.lab  = 1.6)
      } else{
        plot(roc_list[[k]][[j]][[i]], 
             col = colores[i], 
             add = TRUE, 
             lwd = 4)
      }
    }
    
    
    abline(a = 0, b = 1, lty = 2, col = "gray")
    
    legend("bottomright",
           legend = paste0(names(roc_list[[k]][[j]]), " (AUC = ", round(unlist(auc_list[[k]][[j]]), 3), ")"),
           col = colores, 
           lwd = 3,
           cex = 1.5)
    
    
    dev.off()
  }
}


#######################
### ***** FIN ***** ###
#######################
