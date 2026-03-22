# Directorio de trabajo Daniel
getwd()
setwd("C:/Users/driao/Desktop/Proyecto Vanne")

# Directorio de trabajo Vannesa
setwd('C:/Users/lvmgu/OneDrive/Documentos/UNIVERSIDAD MAESTRIA/TESIS/METODO MIRDEEP2')

# Librerias
BiocManager::install("apeglm")
BiocManager::install("multiMiR")

library(tidyverse)
library(reshape2)
library(vsn)
library(DESeq2)
library(corrplot)
library(pheatmap)
library(apeglm)
library(RColorBrewer)
library(multiMiR)


# Objetos -----------------------------------------------------------------

# data_filt: conteos de miRNA filtrados filas y columnas
# DEG_c1: tabla de los degs de NIC2/NIC3 vs SL/NIC1, tiene el padj, logFoldChance,ect
# DEG_c2_nic2: tabla de los degs de NIC2 vs SL/NIC1
# DEG_c2_nic3: tabla de los degs de NIC3 vs SL/NIC1
# tab1: Tabla 1 analoga al articulo. degs en comun entre que NIC2/NIC3 y NIC3 contra SL/NIC1
# vph: Variables clinicas/fenotipicas de las muestras
# list_miRNA_targets: lista de genes para cada miRNA de la comp1 : NIC2/3 vs SL/NIC1
# mat: conteos transformados mediante rlog() depues de la normalizacion de DESeq
# tab2: tabla de contigencias del pca de la 1° comparacion
# tab2: tabla de contigencias del pca de la 2° comparacion

# 1. Cargar Datos ---------------------------------------------------------
library(tidyverse)
library(readxl) # install.packages("readxl")

# Datos de miRNA
data <- read.table('tabla_conteos_mirdeep2_tesis.tsv', sep='\t', header=TRUE) # leer tabla
dim(data) # 2888 x 74

# Sumar los miRNAs que pertenecen a la misma muestra
data <- data %>% group_by(miRNA) %>% summarize_if(is.numeric, sum) %>% 
  column_to_rownames(var = "miRNA")  # tomar nombre de filas de la 1° columna

dim(data) # 2656 x 73
colnames(data) # nombre de las muestras
rownames(data) # nombre de los miRNA

# Datos de VPH
vph <- read_excel("BASE DE DATOS RESULTADOS MICRORNAs.xlsx", range = "A2:J47",
                  col_types = c(rep("guess",3),"numeric",rep("guess",3),"text",rep("guess",2)))
vph2 <- read_excel("BASE DE DATOS RESULTADOS MICRORNAs.xlsx", range = "A51:J84",
                   col_types = c(rep("guess",3),"numeric",rep("guess",3),"text",rep("guess",2)))
vph <- rbind(vph,vph2)
rm(vph2)

head(vph)


# 2. Modificar nombres de columnas ----------------------------------------

# Crear una función para extraer el tipo de condición (NIC o SL) 
# y el número de réplica de los nombres de las columnas
extract_info <- function(colname) {
  if(grepl("NIC.3", colname)) {
    return(c("NIC3", gsub("X(\\d+)_NIC.3.*", "\\1", colname)))
  } else if (grepl("NIC.2", colname)) {
    return(c("NIC2", gsub("X(\\d+)_NIC.2.*", "\\1", colname)))
  } else if (grepl("NIC.1", colname)) {
    return(c("NIC1", gsub("X(\\d+)_NIC.1.*", "\\1", colname)))
  } else if (grepl("SL", colname)) {
    return(c("SL", gsub("X(\\d+)_SL.*", "\\1", colname)))
  } else {
    return(c("Otro", "NA"))
  }
}

# Obtener los nombres de las columnas existentes
nombres_actuales <- colnames(data)

# Crear un nuevo vector de nombres usando la función extract_info
nuevos_nombres <- sapply(nombres_actuales, function(colname) {
  info <- extract_info(colname)
  paste(info[1], info[2], sep = "_")
})

# Asignar los nuevos nombres a las columnas
colnames(data) <- nuevos_nombres


# Tarea: crear un codigo mas corto y mas entendible para que se pueda replicar
# para otros analisis o proyectos en caso que el formato de los primeros nombres cambie

# Obtener los nombres de las columnas existentes
nombres_actuales <- colnames(data)

new_names <- sapply(nombres_actuales, function(nombre) {
  partes = unlist(strsplit(nombre, "_"))       # Dividir el nombre por el guion "_"
  lesion = gsub("\\.", "", partes[2], fixed=F) # Eliminar el punto de la segunda parte
  num_muestra = gsub("X", "", partes[1])       # Extraer el número después de la 'X'
  nuevo_nombre = paste(lesion, num_muestra, sep = "_")  # Crear el nuevo nombre
  return(nuevo_nombre)})

colnames(data) <- new_names
identical(new_names, nuevos_nombres) # comparar con la fucion de Vanne, da lo mismo

# Ordenar las columnas: NIC3, NIC2, NIC1, SL
columnas_ordenadas <- c(
  grep("^NIC3_", colnames(data), value = TRUE),
  grep("^NIC2_", colnames(data), value = TRUE),
  grep("^NIC1_", colnames(data), value = TRUE),
  grep("^SL_", colnames(data), value = TRUE)
)

data <- data %>%
  dplyr::select(all_of(columnas_ordenadas))

colnames(data)


# 3. Cargar datos de Lesion y VPH ---------------------------------------------------------------

lesion <- sub("_.*", "", colnames(data))
table(lesion)
sum(table(lesion))
# NIC1 NIC2 NIC3   SL 
# 29   10   27    7 

num_muestra <- sub(".*_", "", colnames(data))
num_muestra <- as.numeric(num_muestra)
length(num_muestra) # 73 muestras
dim(vph) # 78 x 10

# Ordenar tabla vph segun el orden de las columnas de data
vph <- vph[match(num_muestra, vph$`# Muestra`),]
vph <- vph %>% dplyr::select(-c("# Secuencia"))

# Varibales Fenotipicas
coldata <- data.frame(Muestra = colnames(data), Lesion=lesion,
                      VPH=vph$VPH, Riesgo_lesion=vph$`RIESGO POR LESION`,
                      Riesgo_VPH=vph$`RIESGO POR VPH`)
coldata$VPH_16_18 <- factor(coldata$VPH_16_18, levels = c("16/18","Otros"))

# VPH_16_18: si tiene VPH tipo 16 o 18
coldata$VPH_16_18 <- ifelse(grepl("16|18", coldata$VPH), "16/18", "Otros")

# VPH_16_18: si tiene VPH tipo 16 o 18, tomando por aparte el tipo 16 y 18
coldata$VPH_16y18 <- ifelse(grepl("16", coldata$VPH), "16",
                            ifelse(grepl("18", coldata$VPH), "18", "Otros"))

coldata[coldata$Muestra == "SL_37", ]$VPH_16y18 <- "16/18"
coldata$VPH_16y18 <- factor(coldata$VPH_16y18, levels=c("16","18","16/18","Otros"))

# NOTA: cambia el nombre de estas variables si quieres, para ello
#       pones el nombre que quiera despues del signo $
                        
# ¿Cuantas personas tienen VPH tipo 16 o 18?
table(coldata$VPH_16_18) # 23 muestras

# ¿Cuantas personas tienen VPH tipo 16 o 18, por separado?
table(coldata$VPH_16y18) 

# Tabla de contigencias
table(coldata$Lesion, coldata$VPH_16_18)
table(coldata$Lesion, coldata$VPH_16y18)

# NOTA: Consideradno el tipo 16 y 18 por aparte hay pocas replicas por
#       lesion: SL, NIC1, NIC2, NIC3
# NOTA: Es mejor considerar tipo 16 y 18 unidos


# 4. Paleta de colores ----------------------------------------------------
library(RColorBrewer)
nb.cols <- dim(as.matrix(data))[2]
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

# Colores: SL, NIC1, NIC2, NIC3
colores <- c("lightblue","violet","violetred","red")

#######################################################################
###################### 5. Contro de Calidad #########################

# 5.1 Corrplot ----------------------------------------------------------
library(corrplot)

# NOTA: puedes cambiar la escala de los colores en colorRampPalette()

# Todas las muestras
corrplot(cor(data_filt), cex.main = 1, tl.cex = 0.5, tl.pos = "l",
         title="Corrplot todas las muestras",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,1,0))

# solo NIC3
corrplot(cor(data[1:27]), cex.main = 1, tl.cex = 0.65,
         title="Corrplot muestras NIC 3",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,1.5,0))

# solo NIC2
corrplot(cor(data[28:37]),cex.main = 1, tl.cex = 0.85,
         title="Corrplot muestras NIC 2",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,2,0))

# solo NIC1
corrplot(cor(data[38:66]),cex.main = 1, tl.cex = 0.65,
         title="Corrplot muestras NIC 1",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,1.5,0))

# solo SL
corrplot(cor(data[67:73]),cex.main = 1.2, tl.cex = 0.85,
         title="Corrplot muestras Sin Lesion",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,2,0))


# 5.1 Grafico de barras - promedio de correlaciones -----------------------
cor_matrix <- cor(data)
promedios_muestras <- rowMeans(cor_matrix, na.rm = TRUE)
df_promedios <- data.frame(Muestra= rownames(cor_matrix),
                           Promedio= promedios_muestras,
                           Lesion = lesion)

ggplot(df_promedios, aes(x= Muestra, y= Promedio, fill= Promedio <= 0.5)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values= c("TRUE" = "red", "FALSE" = "blue")) + # Asignar colores
  facet_wrap(vars(Lesion), scales = "free_x") +
  labs(title = "Promedios de Correlación de Muestras",
       x = "Muestras", y = "Promedio de Correlación", fill="corr promedio < 0.5") +
  theme(
    # opciones Ejes
    axis.text.x = element_text(angle = 45, hjust = 1, size=8), # Rotar etiquetas en el eje X para mayor legibilidad
    # opciones Titulo
    plot.title = element_text(vjust=0.5, hjust = 0.5, size=15, face="bold",
                              margin = margin(t=0,r=0,b=15,l=0)),
    # opciones Leyenda
    legend.position = "bottom"
    )

# NOTA: los promedios de correlacion estan calculados contra todas las muestras

# 5.2 Boxplots ------------------------------------------------------------
library(tidyverse)

df <- data %>% pivot_longer(everything(), names_to="Muestra" , values_to="value")
df <- df %>% mutate(Lesion = sub("_.*", "", Muestra))
df$Lesion <- factor(df$Lesion, levels = c("SL","NIC1","NIC2","NIC3"))
df <- df %>% arrange(Lesion)

ggplot(df, aes(x = Muestra, y = value, fill=Lesion)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = colores,
                    labels = c("Sin lesion", "NIC1", "NIC2", "NIC3"),
                    limits = c("SL","NIC1","NIC2","NIC3")) +
  geom_hline(yintercept = 1e5, col="red") + # recta en 100_000
  #scale_x_discrete(paste0("SL_",)) +
  labs(title = "Boxplots todas las muestras",
       x = "Muestra", y="conteos de miRNA", fill="Lesion") +
  theme_light() +
  theme(plot.title = element_text(family="sans",face="bold",size=15,vjust=0.5,hjust=0.5,color="black",angle=0),
        axis.text.x = element_text(angle = 30, hjust = 1, size=9),
        axis.text.y = element_text(size=12),
        legend.position = "right",
        legend.title = element_blank())

boxplot(data, pch=20)
abline(h=1e5, col="red")


###MEANSPLOT###
data_matrix<- as.matrix(data)
meanSdPlot(data_matrix)

#####densidad crudo
data_matrix_melt = melt(data_matrix)[,-1]
colnames(data_matrix_melt) = c("Sample", "Value")
##densidad
ggplot(data_matrix_melt, aes(x=Value)) +
  theme_bw() +
  theme(legend.position="bottom") +
  geom_density(aes(group=Sample, colour=Sample)) +
  labs(x="Intensity", y="Density") +
  scale_color_manual(values = mycolors) +labs(title = "DATOS CRUDOS") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))



# 5.2 Filtros de Muestras -------------------------------------------------

data_filt <- subset(data, select = -c(NIC3_277, NIC3_411))
dim(data_filt) # 2656 x 71

data_filt <- subset(data, select = -c(NIC3_277, NIC3_411))


# 5.3 Filtros de Filas ----------------------------------------------------

# Filtro de filas:
# 1) eliminar filas con suma de conteos < 50 
# 2) eliminar filas con varianza = 0
# 3) eliminar filas con conteos superiores > 100_000

library(matrixStats)
dim(data_filt)

suma_fila <- rowSums(data_filt)
boxplot.stats(suma_fila)$stats # bigote inf, Q1, Q2, Q3, bigote sup
# Q3 = 56
sum(suma_fila == 0)

var_fila <- rowVars(as.matrix(data_filt))
boxplot.stats(var_fila)$stats # Q1=0, Q3 =2.09698189
sum(var_fila == 0)

cv_fila <- sqrt(var_fila) / rowMeans(data_filt)
boxplot.stats(cv_fila)$stats # Q1=1.5752346, Q3= 4.7948385

conteo_max <- apply(data_filt, 1, max)
boxplot.stats(conteo_max)$stats
sum(conteo_max >= 1e5) # 18 genes con valor max superior a 100_000
sum(conteo_max >= 1e6) # 2 genes con valor max superior a 1'000_000
sum(conteo_max == 1e5)

# Filtro
library(matrixStats)
suma_fila <- rowSums(data_filt)
var_fila <- rowVars(as.matrix(data_filt))
conteo_max <- apply(data_filt, 1, max)

data_filt <- data_filt[suma_fila >= 50 &
             var_fila > 0 &
             conteo_max < 1e5,]

dim(data)      # 2656 x 73
dim(data_filt) # 669 x 71, nos quedamos con 669 miRNA y 71 muestras

# boxplot filtrado
boxplot(data_filt, pch=20)

# Hay unos boxplot sin datos atipicos como las otras muestras

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. CONTROL DE CALIDAD despues de FILTROS ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(corrplot)

range(cor(data_filt)) # 0 a 1
dim(data_filt)

lesion <- sub("_.*", "", colnames(data_filt))
table(lesion)
sum(table(lesion))


# Todas las muestras
corrplot(cor(data_filt), cex.main = 1, tl.cex = 0.5, tl.pos = "l",
         title="Corrplot todas las muestras",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,1,0))

# solo NIC3
corrplot(cor(data_filt[1:25]), cex.main = 1, tl.cex = 0.65,
         title="Corrplot muestras NIC 3",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,1.5,0))

# solo NIC2
corrplot(cor(data_filt[26:35]),cex.main = 1, tl.cex = 0.85,
         title="Corrplot muestras NIC 2",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,2,0))

# solo NIC1
corrplot(cor(data_filt[36:64]),cex.main = 1, tl.cex = 0.65,
         title="Corrplot muestras NIC 1",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,1.5,0))

# solo SL
corrplot(cor(data_filt[65:71]),cex.main = 1.2, tl.cex = 0.85,
         title="Corrplot muestras Sin Lesion",
         is.corr = F, tl.col = "black",
         col = colorRampPalette(c("red","white","blue"))(200),
         mar=c(0,0,2,0))


# Boxplot crudo
df2 <- data_filt %>% pivot_longer(everything(), names_to="Muestra" , values_to="value")
df2 <- df2 %>% mutate(Lesion = sub("_.*", "", Muestra))
df2$Lesion <- factor(df2$Lesion, levels = c("SL","NIC1","NIC2","NIC3"))
df2 <- df %>% arrange(Lesion)

ggplot(df2, aes(x = Muestra, y = value, fill=Lesion)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = colores,
                    labels = c("Sin lesion", "NIC1", "NIC2", "NIC3"),
                    limits = c("NIC3","NIC2","NIC1","SL")) +
  geom_hline(yintercept = 1e5, colour="red", linetype="dashed") +
  labs(title = "Boxplots datos filtrados todas las muestras",
       x = "Muestra", y="conteos de miRNA", fill="Lesion") +
  theme_light() +
  theme(plot.title = element_text(family="sans",face="bold",size=15,vjust=0.5,hjust=0.5,color="black",angle=0),
        axis.text.x = element_text(angle = 30, hjust = 1, size=7),
        axis.text.y = element_text(size=12),
        legend.position = "right",
        legend.title = element_blank())

# Columnas con conteos bajos
conteo_max2 <- apply(data_filt, 2, max)
which(conteo_max2 < 10000)


# Boxplot en escala logaritmica: FALTA MIRAR CUAL ES LA TRANSFORMACION
library(DESeq2)
ggplot(df2, aes(x = Muestra, y = , fill=Lesion)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = colores,
                    labels = c("Sin lesion", "NIC1", "NIC2", "NIC3"),
                    limits = c("NIC3","NIC2","NIC1","SL")) +
  geom_hline(yintercept = 1e5, colour="red", linetype="dashed") +
  labs(title = "Boxplots todas las muestras",
       x = "Muestra", y="conteos de miRNA", fill="Lesion") +
  theme_light() +
  theme(plot.title = element_text(family="sans",face="bold",size=15,vjust=0.5,hjust=0.5,color="black",angle=0),
        axis.text.x = element_text(angle = 30, hjust = 1, size=8),
        axis.text.y = element_text(size=12),
        legend.position = "right",
        legend.title = element_blank())


# Grafico de densidad
ggplot(df2, aes(x = value, group = Muestra, color=Lesion)) +
  geom_density(position = "identity", linewidth=0.8) +
  scale_color_manual(values = colores,
                     labels = c("SL", "NIC1","NIC2","NIC3")) +
  labs(title= "Grafico de densidad de las muestras",
       x="Conteos de miRNA", y="Densidad") +
  theme_minimal()

# mean vs SD plot
library(vsn)
meanSdPlot(as.matrix(data_filt))

# Heat map
library(pheatmap)
pheatmap(as.matrix(data_filt[1:50,]))


# 7. Comparar NIC< contra NIC2/NIC3 --------------------------------------------------
library(DESeq2)

dim(data_filt)   # 669 x 71

lesion <- sub("_.*", "", colnames(data_filt))
table(lesion)
sum(table(lesion))

# Creamos la tabla agrupando las condiciones
# condition 1: 'SL/NIC1' y 'NIC2/3'
# condition 2: 'SL/NIC1', 'NIC2', 'NIC3'

coldata2 <- coldata %>% filter(Muestra %in% colnames(data_filt))
coldata2$condition1 <- factor( c(rep('NIC2/3', 35), rep('SL/NIC1', 36)))
coldata2$condition2 <- factor( c(rep('NIC3',25), rep('NIC2',10), rep('SL/NIC1', 36)))

coldata2$condition1 <- relevel(coldata2$condition1, ref = 'SL/NIC1') # nivel de ref
coldata2$condition2 <- relevel(coldata2$condition2, ref = 'SL/NIC1') # nivel de ref
rownames(coldata2) <- colnames(data_filt)

# Tabla de contigenicas
table(coldata2$condition1, coldata2$VPH_16_18)
#            16/18  Otros
# SL/NIC1    12     24
# NIC2/3     11     24

table(coldata2$condition2, coldata2$VPH_16_18)
#            16/18  Otros
# SL/NIC1    12     24
# NIC2        4      6
# NIC3        7     18


## 7.1 Primera Comparacion -------------------------------------------------

# Crear el objeto DESeq
dds_c1 <- DESeqDataSetFromMatrix(
  countData = data_filt,
  colData = coldata2,
  design= ~ 1 + condition1)

table(coldata2$condition1)

# Estimar los sizefactors
dds_c1 <- estimateSizeFactors(dds_c1)
sizeFactors(dds_c1)

sum(sizeFactors(dds_c1) >= 5) # posiblemente eliminar 7 muestras
sum(sizeFactors(dds_c1) <= 0.1) # posiblemente eliminar 7 muestras

# Estimacion de las dispersiones
dds_c1 <- estimateDispersions(dds_c1)
plotDispEsts(dds_c1)

# ¿Cuales son esas muestras son sizefactors raros?
which(sizeFactors(dds_c1) >= 5)   # esto te bota el nombre y la posicion de la muestra
which(sizeFactors(dds_c1) <= 0.1)

conteos_columna <- apply(data_filt, 2, max)
which(conteos_columna < 10000)


head(counts(dds0))
head(counts(dds0, normalized=TRUE))


## 7.2  DEGs y Volcano comparacion 1  -------------------------------------------------
dds_c1 <- DESeq(dds_c1, test="Wald")
resultsNames(dds_c1) # "condition1_NIC2.3_vs_SL.NIC1"
colData(dds_c1)

res_c1 <- results(dds_c1, contrast=list("condition1_NIC2.3_vs_SL.NIC1"),
                  pAdjustMethod = "BH", alpha=0.05, tidy=F) # modificar pvalor

# MAplot
plotMA(res_c1, alpha=0.05,  colNonSig ="gray60", colSig = "red", cex=0.7)
summary(res_c1)

# Volcano Plot: Creo que este es mejor
c1 <- as.data.frame(res_c1) %>% arrange(padj)

c1$expressed <- "NO"
c1$expressed[c1$log2FoldChange > 1   & c1$padj < 0.05] <- "UP"
c1$expressed[c1$log2FoldChange < -1   & c1$padj < 0.05] <- "DOWN"

c1 %>% filter(log2FoldChange < 1 & padj < 0.05)

# linea horizontal
limite_h <- -log10(max(c1[c1$log2FoldChange > 1   & c1$padj < 0.05,]$pvalue, na.rm=T))

# Nombres de miRNA expresados
c1$label <- NA
c1$label[c1$expressed != "NO"] <- rownames(c1)[c1$expressed != "NO"]

# Grafico Volcano Plot
library(ggrepel)
ggplot(c1, aes(x=log2FoldChange, y= -log10(pvalue), col = expressed, label=label)) +
  geom_point(size=2) +
  scale_color_manual(values = c("blue","black","red")) +
  geom_vline(xintercept = c(-1,1), col="gray40", linetype="dashed", linewidth=1) +
  geom_hline(yintercept = limite_h, col="gray40", linetype="dashed", linewidth=1) +
  geom_text_repel() +
  labs(title="Volcano plot - NIC2/3 vs SL/NIC1", color="Diferentially") +
  theme_light() +
  theme(
    # Titulo
    plot.title = element_text(hjust=0.5, vjust=0.5, size=18),
    axis.title = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Observaciones
# cambiar tamaño de puntos: geom_point(size= 3)
# cambiar el color de los: scale_color_manual(values = c("DOWN","NO","UP"))
# cambiar color de las lineas punteadas: geom_vline( col=), geom_hline( col=)
# cambiar tamaño de la leyenda: gides( guide_legend( ...., size=))

## 7.4 Tabla de DEGs comparacion 1 -----------------------------------------

DEG_c1 <- c1 %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
DEG_c1$FC <- 2^(DEG_c1$log2FoldChange) # agregar FC
dim(DEG_c1) # 20
rownames(DEG_c1)

DEG_c1 <- DEG_c1 %>% dplyr::select("baseMean","FC","log2FoldChange","lfcSE",
                            "stat","pvalue","padj","expressed","label")

CV <- function(x){100*sd(x)/mean(x)}
DEG_c1$CV <- apply(data_filt[rownames(DEG_c1),], 1, CV)

#install.packages("openxlsx")
library(openxlsx)

# Guardar DEGs comparacion 1
write.xlsx(DEG_c1, "DEG_comp1.xlsx", rowNames=T)

## 7.5 Con Shrinkage - NO hay diferencia -----------------------------------

# Con shrinkage
library("apeglm")
resultsNames(dds_c1)
res_c1_sh <- lfcShrink(dds_c1, coef="condition1_NIC2.3_vs_SL.NIC1", type="apeglm")
res_c1_sh
class(res_c1_sh)

# valores P sin y con shrinkage
Padj = data.frame(padj = res_c1$padj, sh = "NO")
Padj2 = data.frame(padj = res_c1_sh$padj, sh = "YES")
Padj = rbind(Padj, Padj2)

ggplot(Padj, aes(x=padj)) + 
  geom_histogram() +
  facet_wrap(~sh)  +
  labs(y='frecuencia', x="p valor ajustado") 

plotCounts(dds_c1, gene=which.min(res_c1$padj), intgroup="condition1")
plotCounts(dds_c1, gene=which.min(res_c1_sh$padj), intgroup="condition1")


## 7.6 Heatmap -------------------------------------------------------------

# Transformacion rlog
rld <-  rlog(dds_c1, blind = TRUE)

mat <- assay(rld) # obtener matriz de datos
class(mat)

up_c1 <- DEG_c1 %>% filter(expressed == "UP") %>% rownames()
down_c1 <- DEG_c1 %>% filter(expressed == "DOWN") %>% rownames()

anot <- colData(dds_c1)[, c("condition2","VPH_16_18")]
anot <- as.data.frame(anot)
rownames(anot) <- colnames(mat)
colnames(anot) <- c("Condicion", "Tipo de VPH")

anot_colors <- list(
  "Tipo de VPH" = c("16/18"="limegreen",
                "Otros"="navy"),
  Condicion = c("SL/NIC1"="lightblue",
                 "NIC2"="violet",
                 "NIC3"="red")
)

library(pheatmap)
pheatmap(mat[up_c1,] , annotation_col = anot, show_colnames = T,
         fontsize_col=7, annotation_colors= anot_colors,
         main="Heatmap miRNA UP from NIC2/3 vs SL/NIC1")

pheatmap(mat[down_c1,] , annotation_col = anot, show_colnames = T,
         fontsize_col=7, annotation_colors= anot_colors,
         main="Heatmap miRNA DOWN from NIC2/3 vs SL/NIC1")

# Tipo VPH
# Condicion


# 8. Comparar NIC1< contra NIC2 y NIC3 por aparte ----------------------------------------------------------
dds_c2 <- DESeqDataSetFromMatrix(
  countData = data_filt,
  colData = coldata2,
  design= ~ condition2)

table(coldata2$condition2)

dds_c2 <- estimateSizeFactors(dds_c2)
sizeFactors(dds_c2)

# son los mismo sizefactors de la primera comparacion
identical(sizeFactors(dds_c1), sizeFactors(dds_c2))

# plot mean vs Sd
dds_c2 <- estimateDispersions(dds_c2)
plotDispEsts(dds_c2) # estimaciones dispersiones


## 8.1 DEGs y Volcano ----------------------------------------------------------

dds_c2 <- DESeq(dds_c2, test="LRT", reduced= ~ 1)
resultsNames(dds_c2) # coeficientes

# Resultados
res_c1_nic3 <- results(dds_c2, contrast=list("condition2_NIC3_vs_SL.NIC1"),
                       pAdjustMethod = "BH", alpha=0.05, tidy=F)
res_c1_nic2 <- results(dds_c2, contrast=list("condition2_NIC2_vs_SL.NIC1"),
                       pAdjustMethod = "BH", alpha=0.05, tidy=F)

# MAplot
plotMA(res_c1_nic3, alpha=0.05,  colNonSig ="gray60", colSig = "red", cex=0.7)
summary(res_c1_nic3)

plotMA(res_c1_nic2, alpha=0.05,  colNonSig ="gray60", colSig = "red", cex=0.7)
summary(res_c1_nic2)


# Volcano
c2_nic3 <- as.data.frame(res_c1_nic3) %>% arrange(padj)
c2_nic2 <- as.data.frame(res_c1_nic2) %>% arrange(padj)

# contra NIC3
c2_nic3$expressed <- "NO"
c2_nic3$expressed[c2_nic3$log2FoldChange > 1    & c2_nic3$padj < 0.05] <- "UP"
c2_nic3$expressed[c2_nic3$log2FoldChange < -1   & c2_nic3$padj < 0.05] <- "DOWN"

# contra NIC2
c2_nic2$expressed <- "NO"
c2_nic2$expressed[c2_nic2$log2FoldChange > 1    & c2_nic2$padj < 0.05] <- "UP"
c2_nic2$expressed[c2_nic2$log2FoldChange < -1   & c2_nic2$padj < 0.05] <- "DOWN"
table(c2_nic2$expressed)


# Grafico Volcano: NIC3 vs SL/NIC1
limite_h <- -log10(max(c2_nic3[c2_nic3$log2FoldChange > 1   & c2_nic3$padj < 0.05,]$pvalue, na.rm=T))

# Nombres de miRNA expresados
c2_nic3$label <- NA
c2_nic3$label[c2_nic3$expressed != "NO"] <- rownames(c2_nic3)[c2_nic3$expressed != "NO"]
table(c2_nic3$expressed)

ggplot(c2_nic3, aes(x=log2FoldChange, y= -log10(pvalue), col = expressed, label=label)) +
  geom_point(size=2) +
  scale_color_manual(values = c("blue","black","red")) +
  geom_vline(xintercept = c(-1,1), col="darkblue", linetype="dashed", linewidth=1) +
  geom_hline(yintercept = limite_h, col="darkblue", linetype="dashed", linewidth=1) +
  geom_text_repel() +
  labs(title="Volcano plot - NIC3 vs SL/NIC1", color="Diferentially") +
  theme_light() +
  theme(
    #Titulo
    plot.title = element_text(hjust=0.5, vjust=0.5, size=18),
    axis.title = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Grafico Volcano: NIC2 vs SL/NIC1
limite_h <- -log10(max(c2_nic2[c2_nic2$log2FoldChange < - 1   & c2_nic2$padj < 0.05,]$pvalue, na.rm=T))

# Nombres de miRNA expresados
c2_nic2$label <- NA
c2_nic2$label[c2_nic2$expressed != "NO"] <- rownames(c2_nic2)[c2_nic2$expressed != "NO"]
table(c2_nic2$expressed)


ggplot(c2_nic2, aes(x=log2FoldChange, y= -log10(pvalue), col = expressed, label=label)) +
  geom_point(size=2) +
  scale_color_manual(values = c("blue","black","red")) +
  geom_vline(xintercept = c(-1,1), col="gray40", linetype="dashed", linewidth=1) +
  geom_hline(yintercept = limite_h, col="gray40", linetype="dashed", linewidth=1) +
  geom_text_repel() +
  labs(title="Volcano plot - NIC2 vs SL/NIC1", color="Diferentially") +
  theme_light() +
  theme(
    #Titulo
    plot.title = element_text(hjust=0.5, vjust=0.5, size=18),
    axis.title = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 5)))


## 8.2 Tabla de DEGs comparacion 2 -----------------------------------------

# Crear tabla de resultados y agregar Fold Change
DEG_c2_nic3 <- c2_nic3 %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
dim(DEG_c2_nic3) # 23
DEG_c2_nic3$FC <- 2^(DEG_c2_nic3$log2FoldChange) # agregar FC

DEG_c2_nic2 <- c2_nic2 %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
DEG_c2_nic2$FC <- 2^(DEG_c2_nic2$log2FoldChange) # agregar FC
dim(DEG_c2_nic2) # 15

# Ordenar columnas
DEG_c2_nic2 <- DEG_c2_nic2 %>% dplyr::select("baseMean","FC","log2FoldChange","lfcSE",
                                   "stat","pvalue","padj","expressed","label")

DEG_c2_nic3 <- DEG_c2_nic3 %>% dplyr::select("baseMean","FC","log2FoldChange","lfcSE",
                                             "stat","pvalue","padj","expressed","label")

# Guardar Tablas
write.xlsx(DEG_c2_nic2, "DEG_comp2_NIC2.xlsx", rowNames=T)
write.xlsx(DEG_c2_nic3, "DEG_comp2_NIC3.xlsx", rowNames=T)


# Intersecciones
intersect(rownames(DEG_c1), rownames(DEG_c2_nic3)) # 20 en comun
intersect(rownames(DEG_c1), rownames(DEG_c2_nic2)) # 6 en comun

intersect(rownames(DEG_c2_nic2), rownames(DEG_c2_nic3)) # 6 en comun



## 8.3 Heatmap  ----------------------------------------------------------------

rld2 <- rlog(dds_c2, blind = T)
mat2 <- assay(rld2)

identical(mat, mat2) # es lo mismo

# up y down genes
up_c2_nic2 <- DEG_c2_nic2 %>% filter(expressed == "UP") %>% rownames() # no hay genes a la alta
down_c2_nic2 <- DEG_c2_nic2 %>% filter(expressed == "DOWN") %>% rownames()

anot <- colData(dds_c1)[, c("condition2","VPH_16_18")]
anot <- as.data.frame(anot)
rownames(anot) <- colnames(mat)
colnames(anot) <- c("Condicion", "Tipo de VPH")

anot_colors <- list(
  "Tipo de VPH" = c("16/18"="limegreen",
                    "Otros"="navy"),
  Condicion = c("SL/NIC1"="lightblue",
                "NIC2"="violet",
                "NIC3"="red")
)

library(pheatmap)
pheatmap(mat[down_c2_nic2,] , annotation_col = anot, show_colnames = T,
         fontsize_col=7, annotation_colors= anot_colors,
         main="Heatmap miRNA DOWN from NIC2 vs SL/NIC1")

# NIC3
up_c2_nic3 <- DEG_c2_nic3 %>% filter(expressed == "UP") %>% rownames() # no hay genes a la alta
down_c2_nic3 <- DEG_c2_nic3 %>% filter(expressed == "DOWN") %>% rownames()

pheatmap(mat[up_c2_nic3,] , annotation_col = anot, show_colnames = T,
         fontsize_col=7, annotation_colors= anot_colors,
         main="Heatmap miRNA UP from NIC3 vs SL/NIC1")

pheatmap(mat[down_c2_nic3,] , annotation_col = anot, show_colnames = T,
         fontsize_col=7, annotation_colors= anot_colors,
         main="Heatmap miRNA DOWN from NIC3 vs SL/NIC1")

# 9. Interseccion ---------------------------------------------------------

miRNA_ints <- intersect(rownames(DEG_c1), rownames(DEG_c2_nic3))
length(miRNa_ints)

a <- DEG_c1[miRNA_ints,]
b <- DEG_c2_nic3[miRNA_ints,]
  
tab1 <- data.frame(
  "FC_NIC2+" = a$FC,
  "pvalue_NIC2+" = a$padj,
  "FC_NIC3" = b$FC,
  "pvalue_NIC3" = b$padj,
  "CV" = apply(mat[miRNa_ints,], 1, CV),
  "Expressed" = a$expressed
)

# Guardar Tabla 1
write.xlsx(tab1, "Tabla1.xlsx", rowNames=T)


# 9. PCA sobre muestras ---------------------------------------------------------------------

plotPCA(rld, intgroup = "condition1")
plotPCA(rld, intgroup = "condition2") # las muestras no se diferencian

library(ade4)
library(factoextra)
pca <- dudi.pca(t(data_filt[rownames(DEG_c1),]), scannf = FALSE, nf = 2)
summary(pca2)

# Variabilidad retenida
fviz_eig(pca)
pve <- 100*pca$eig/sum(pca$eig)
cumsum(pve) # 57.2 % las 2 primeras componentes

# Circulo de correlaciones
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             cex = 1
)

# Graficos PCA
set.seed(12)
cl <- kmeans(t(data_filt[up_c1,]), 4)
s.class(pca$li, as.factor(cl$cluster), col=c("#00AFBB","red"))
s.class(pca$li, as.factor(coldata2$condition1), col=c("#00AFBB","red"))

tab2 <- table(as.factor(coldata2$condition1), as.factor(cl$cluster))


pca2 <- dudi.pca(t(data_filt[rownames(DEG_c2_nic3),]), scannf = FALSE, nf = 2)
summary(pca2)

set.seed(1234)
cl2 <- kmeans(t(data_filt[up_c2_nic3,]), 3)
s.class(pca2$li, as.factor(cl2$cluster), col=c("#00AFBB","orange","red"))
s.class(pca2$li, as.factor(coldata2$condition2), col=c("#00AFBB","orange","red"))

tab3 <- table(as.factor(coldata2$condition2), as.factor(cl2$cluster))
tab3


# 10. PCA sobre genes ---------------------------------------------------------

pca2 <- dudi.pca(data_filt[rownames(DEG_c1),], scannf = FALSE, nf = 2)
summary(pca2)

set.seed(9)
cl3 <- kmeans(data_filt[rownames(DEG_c1),], 2)
s.class(pca2$li, as.factor(cl3$cluster), col=c(1,2))

# NOTA: hay muy pocos DEGs para hacerlo



# 11. LLamado de Genes Blanco ---------------------------------------------

#https://www.bioconductor.org/packages/release/bioc/html/multiMiR.html

targets_genes <- get_multimir(org = "Homo Sapiens", # organismo
                              mirna = rownames(DEG_c1), # lista de miRNA
                              table = "validated") # interaciones validadas

head(targets_genes@data)
dim(targets_genes@data) # 19615 x 10
colnames(targets_genes@data)
columns(targets_genes)
# database: base de datos de la interaccion
# mature_mirna_id: codigo de miRNA
# "target_symbol": simbolo del gen
# target_entrez
# target_ensembl

# Lista de Genes asociados a los miRNA
targets <- targets_genes@data

list_miRNA_targets <- list()

for(mi in miRNa_ints){
  list_miRNA_targets[[mi]] <- targets %>% filter(mature_mirna_id == mi) %>% 
    dplyr::select(target_entrez) %>% unique()
}
# reemplazar target_entrez por target_symbol o target_ensembl

# Para acceder a los genes de de un miRNA, list_miRNA_targets$ y seleccionar el miRNA
list_miRNA_targets$`hsa-miR-96-5p` # genes asociados a hsa-miR-96-5p

# Todos los Genes asociados a todos los miRNA
all_genes_target <- unique(unlist(list_miRNA_targets, use.names=F))
length(all_genes_target) # 8718


##con los genes target ponerlos en lista para cada microRNA


