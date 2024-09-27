## ----setup, echo=FALSE, cache=FALSE-------------------------------------------
library(knitr) ## kable()
library(kableExtra) ## kable_styling(), save_kable()
library(here) ## here()
library(usethis) ## use_directory()

knitr::opts_chunk$set(
  collapse=TRUE,
  comment="",
  fig.align="center",
  fig.wide=TRUE,
  cache=FALSE
)

## this option avoid use_directory() being verbose
options(usethis.quiet=TRUE)

## create these paths at build time if they do not exist
use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))

## fetch the package root directory
path2pkg <- here()

## ----eval=FALSE---------------------------------------------------------------
#  devtools::build_vignettes()

## ----eval=FALSE---------------------------------------------------------------
#  devtools::document()

## ----message=FALSE------------------------------------------------------------
library(SummarizedExperiment)

se <- readRDS(file.path(system.file("extdata", package="IEOproject"), "GSE181674.rds"))
se



## -----------------------------------------------------------------------------
head(rowData(se))

colnames(rowData(se))

dim(rowData(se))

## -----------------------------------------------------------------------------
head(colData(se), n=3)

colnames(colData(se))

dim(colData(se))

## -----------------------------------------------------------------------------
ncol(se)
length(unique(se$geo_accession))
table(lengths(split(colnames(se), se$geo_accession)))

## ----message=FALSE------------------------------------------------------------
library(edgeR)
dge <- DGEList(counts=assays(se)$counts, genes=rowData(se), samples=colData(se))
dim(dge)
dim(dge)==dim(se)

## -----------------------------------------------------------------------------
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.25)
assays(se)$logCPM[1:3, 1:6] #Visualization

## -----------------------------------------------------------------------------
table(se$organism_ch1)

table(se$`tissue:ch1`)

## -----------------------------------------------------------------------------
table(se$extract_protocol_ch1)

table(se$data_processing.1)

## -----------------------------------------------------------------------------
table(se$`disease state:ch1`)

## -----------------------------------------------------------------------------
se$disease_state <- as.factor(se$`disease state:ch1`)
levels(se$disease_state) <- c("ND", "T1D", "AAb")

dge$samples$disease_state <- as.factor(dge$samples$disease.state.ch1)
levels(dge$samples$disease_state) <- c("ND", "T1D", "AAb")

se$bioreplicate <- se$`biorep:ch1`

#Visualization
se$disease_state
dge$samples$disease_state

## ----pheno, echo=FALSE, message=FALSE-----------------------------------------
tmpdf <- data.frame("Identifer"=colnames(se),
                    "Sample Title" = se$title,
                    "Disease State"=se$`disease state:ch1`,
                    "Bioreplicate"=se$bioreplicate,
                    check.names=FALSE)
ktab<-kable(tmpdf, caption="Phenotypic variables", format = "html", booktabs = TRUE, row.names = FALSE)

kable_styling(ktab, position="center")

## ----frecuency, echo=FALSE, message=FALSE-------------------------------------

freq_table <- as.data.frame(table(se$disease_state))
colnames(freq_table) <- c("Disease state", "Frequency")
kable_styling(kable(freq_table, caption = "Frequency table for disease state", bootstrap_options = "striped", full_width = FALSE, format = "html", booktabs = TRUE, row.names = FALSE))



## ----libsizes, echo=FALSE, out.width="90%", fig.align="center", fig.cap="Library sizes in increasing order."----
par(mar=c(7, 5, 2, 2))
ord <- order(dge$sample$lib.size/1e6)
bp <- barplot(dge$sample$lib.size[ord]/1e6, las = 1, ylab = "Millions of reads",
              xlab = "", col = se$disease_state[ord], las = 2)
axis(1, at = bp, labels = se$title[ord], las = 2, cex.axis=0.65)
legend("topleft", c("Diabetic", "Control", "AAb"), fill = unique(se$disease_state[ord]), inset = 0.01, cex = 0.65)
mtext("Samples", side = 1, line = 4, cex.lab = 1.2)

## ----distRawExp, echo=FALSE, fig.height=5, fig.width=5, out.width="600px", fig.cap="Non-parametric density distribution of expression profiles per sample.", message=FALSE----
library(geneplotter)
par(mar=c(4, 5, 1, 1))
lst <- as.list(as.data.frame(assays(se)$logCPM))
multidensity(lst, xlab=expression(log[2] * "CPM"), legend=NULL, main="", las=1)
legend("topright", legend = se$title, fill = 1:length(colnames(se)), cex=0.5)

## ----boxplot, echo=FALSE, fig.height=5, fig.width=5, out.width="600px", fig.cap="Boxplot of expression profiles per sample", message=FALSE----

boxplot(assays(se)$logCPM, col="gray", xlab="Samples",
          ylab=expression(log[2] * "CPM"), cex.axis=0.5, cex.lab=1.0, las=2, names=se$title)


## ----exprdist, echo=FALSE, out.width="600px",fig.align="center",  fig.cap="Distribution of average expression level per gene."----
avgexp <- rowMeans(assays(se)$logCPM)
h<-hist(avgexp, xlab=expression(log[2] * "CPM"), main="", las=1, col = "gray")
abline(v=1, col="red", lwd=2)
axis(1, at = seq(-8, max(avgexp)), labels = TRUE)
axis(1, at = h$mids, labels = TRUE)

## -----------------------------------------------------------------------------
mask <- filterByExpr(dge, group=se$disease_state)
se.filt <- se[mask, ]
dim(se.filt)
dge.filt <- dge[mask, ]
dim(dge.filt)

## ----filtering1, echo=FALSE, out.width="600px", fig.align="center", fig.cap="Distribution of average expression level per gene with filtered genes highlighted in red 1"----
par(mar=c(4, 5, 1, 1))
h<-hist(avgexp, xlab=expression(log[2] * "CPM"), main="", las=1, col = "gray")
x <- cut(rowMeans(assays(se.filt)$logCPM), breaks=h$breaks)
lines(h$mids, table(x), type="h", lwd=10, lend=1, col="red")
legend("topright", c("All genes", "Filtered genes"), fill=c("grey", "red"))
axis(1, at = seq(-8, max(avgexp)), labels = TRUE)
axis(1, at = h$mids, labels = TRUE)

## -----------------------------------------------------------------------------
mask <- rowMeans(assays(se)$logCPM) > 1
se.filt <- se[mask, ]
dge.filt <- dge[mask, ]
dim(dge.filt)
par(mar=c(4, 5, 1, 1))

## ----filtering2, echo=FALSE, out.width="600px", fig.align="center",  fig.cap="Distribution of average expression level per gene with filtered genes highlighted in red 2"----
par(mar=c(4, 5, 1, 1))
h<-hist(avgexp, xlab=expression(log[2] * "CPM"), main="", las=1, col = "gray")
x <- cut(rowMeans(assays(se.filt)$logCPM), breaks=h$breaks)
lines(h$mids, table(x), type="h", lwd=10, lend=1, col="red")
legend("topright", c("All genes", "Filtered genes"), fill=c("grey", "red"))
axis(1, at = seq(-8, max(avgexp)), labels = TRUE)
axis(1, at = h$mids, labels = TRUE)

## -----------------------------------------------------------------------------
dge.filt <- calcNormFactors(dge.filt)

## -----------------------------------------------------------------------------
assays(se.filt)$logCPM <- cpm(dge.filt, log=TRUE, normalized.lib.sizes=TRUE)

## ----maPlots, echo=FALSE, message=FALSE, fig.height=18, fig.width=10, out.width="100%", dpi=100, echo=FALSE, fig.cap="MA-plots of expression values before and after filtering and normalisation."----
# Non normalized data
dge$samples$group <- dge$samples$disease_state
# Define the comparison groups
group1 <- "AAb"
group2 <- "ND"
group3 <- "T1D"

# Create subsets of data for each pairwise comparison
subset1 <- dge[, dge$samples$group %in% c(group1, group2)]
subset2 <- dge[, dge$samples$group %in% c(group1, group3)]
subset3 <- dge[, dge$samples$group %in% c(group2, group3)]

# In normalized and filtered data
dge.filt$samples$disease_state <- as.factor(dge.filt$samples$disease.state.ch1)
levels(dge.filt$samples$disease_state) <- c("AAb", "ND", "T1D")
dge.filt$samples$group <- dge.filt$samples$disease_state

# Create subsets of data for each pairwise comparison
subset1.filt <- dge.filt[, dge.filt$samples$group %in% c(group1, group2)]
subset2.filt <- dge.filt[, dge.filt$samples$group %in% c(group1, group3)]
subset3.filt <- dge.filt[, dge.filt$samples$group %in% c(group2, group3)]

# Create MA plots for each comparison, suppressing warnings
par(mfrow=c(3,2))

# AAb and ND
suppressWarnings({
  plotSmear(subset1, lowess=TRUE, las=1, cex.lab=1.5, cex.axis=1.2)
  abline(h=0, col="blue", lwd=2)
  title(main="AAb and ND")

  plotSmear(subset1.filt, lowess=TRUE, las=1, cex.lab=1.5, cex.axis=1.2)
  abline(h=0, col="blue", lwd=2)
  title(main="Normalized/Filtered AAb and ND")
})

# AAb and T1D
suppressWarnings({
  plotSmear(subset2, lowess=TRUE, las=1, cex.lab=1.5, cex.axis=1.2)
  abline(h=0, col="blue", lwd=2)
  title(main="AAb and T1D")

  plotSmear(subset2.filt, lowess=TRUE, las=1, cex.lab=1.5, cex.axis=1.2)
  abline(h=0, col="blue", lwd=2)
  title(main="Normalized/Filtered AAb and T1D")
})

# ND and T1D
suppressWarnings({
  plotSmear(subset3, lowess=TRUE, las=1, cex.lab=1.5, cex.axis=1.2)
  abline(h=0, col="blue", lwd=2)
  title(main="ND and T1D")

  plotSmear(subset3.filt, lowess=TRUE, las=1, cex.lab=1.5, cex.axis=1.2)
  abline(h=0, col="blue", lwd=2)
  title(main="Normalized/Filtered ND and T1D")
})



## ----maPlots2, fig.height=18, fig.width=10, dpi=100,out.width="100%", echo=FALSE, fig.cap="MA-plots of filtered and normalized expression values for each sample AAb vs T1D."----

par(mfrow=c(6, 3), mar=c(4, 5, 3, 1))
for (i in 1:ncol(se.filt)) {
  
  A <- rowMeans(assays(se.filt)$logCPM)
  M <- assays(se.filt)$logCPM[, i] - A
  sample_title <- se.filt$title[i]  
  smoothScatter(A, M, main=sample_title, las=1, cex.axis=1.2,
                cex.lab=1.5, cex.main=2)
  abline(h=0, col="blue", lwd=2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col="red", lwd=2)
}


## ----mdsPlot, fig.height=5, fig.width=14, dpi=100, out.width="100%", echo=FALSE, fig.cap="Multidimensional scaling plot (MDS) and hierarchical clustering by gene expression of the samples. Labels correspond to DiseaseState_sampleID and colors indicate disease state."----
par(mfrow = c(1, 2))

colnames(dge.filt)<-dge.filt$samples$title
plotMDS(dge.filt, col = c("red", "blue", "darkgreen")[dge.filt$samples$disease_state], cex = 0.5)
legend("topright", c("ND", "T1D", "AAb"), fill = c("red", "blue", "darkgreen"), inset = 0)
title(main = "MDS colored by disease state")


logCPM <- cpm(dge.filt, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
disease_colors <- c("red", "blue", "darkgreen")  
disease_state_colors <- disease_colors[as.numeric(factor(se.filt$disease_state))]
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
batch <- disease_state_colors
names(batch) <- se$title
outcome <- se$title
names(outcome) <- se$title
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)

plot(sampleDendrogram, main = "Hierarchical clustering colored by disease state", cex = 0.6)
legend("right", c("ND", "T1D", "AAb"), fill=c("red", "blue", "darkgreen"), inset=0)


## -----------------------------------------------------------------------------

sex <- c("female","male","male","male","female","female","male","female","female","male","female","female","male","male","female","male","female","male")
cause_death <- c("Anoxia", "Head Trauma", "Head Trauma", "Anoxia", "Head Trauma", "Head Trauma", "Head Trauma", "Anoxia", "Head Trauma", "Head Trauma", "Head Trauma", "Cerebral Edema", "Anoxia", "Head Trauma", "Anoxia", "Head Trauma", "Cerebral Edema", "Anoxia")
age <- c(23.0, 17.65, 22.10, 24.2, 14.0, 20.0, 19.0, 9, 14.3, 23.3, 25.07, 5.0, 13.0, 19.0, 14.0, 24.9, 12.5, 23.1)
bmi <- c(18.60, 20.80, 28.90, 34.00, 25.60, 25.00, 24.80, 31.30, 51.40, 19.60, 
         23.50, 15.90, 16.60, 17.40, 24.30, 24.50, 28.40, 28.50)

colnames(dge.filt)<-dge.filt$samples$title

colData(se.filt)$sex <- as.factor(sex)
colData(se.filt)$cause_death <- as.factor(cause_death)

colData(se.filt)$age <- as.numeric(age)
colData(se.filt)$bmi <- as.numeric(bmi)

dge.filt$samples$sex <- as.factor(sex)
dge.filt$samples$cause_death <- as.factor(cause_death)
dge.filt$samples$age <- as.numeric(age)
dge.filt$samples$bmi <- as.numeric(bmi)

## ----mdspheno, fig.height=10, fig.width=15, dpi=100, out.width="100%", echo=FALSE, fig.cap="Multidimensional scaling plot (MDS) of the samples. Labels correspond to DiseaseState_sampleID and colors indicate the different potential confounding variables"----
par(mfrow=c(2, 2))
#By sex
plotMDS(dge.filt, col=c("red", "blue")[as.factor(dge.filt$samples$sex)])
legend("topright", c("female", "male"), fill=c("red", "blue"), inset = 0)
title(main = "By sex")

#By cause of death
plotMDS(dge.filt, col=c("red", "blue", "darkgreen")[as.factor(dge.filt$samples$cause_death)])
legend("topright", c("Anoxia", "Cerebral Edema", "Head Trauma"), fill=c("red", "blue", "darkgreen"), inset = 0)
title(main = "By cause of death")

#By age
unique_ages <- unique(dge.filt$samples$age)
sorted_unique_ages <- sort(unique_ages)
color_range <- colorRampPalette(c("lightblue", "darkblue"))
age_colors <- color_range(length(sorted_unique_ages))
age_color_map <- age_colors[match(dge.filt$samples$age, sorted_unique_ages)]

# Plot MDS with Age color gradient
plotMDS(dge.filt, col = age_color_map)
title(main = "By Age")
# Legend for BMI color gradient
legend("topright", legend = sorted_unique_ages, fill = age_colors, cex = 0.6)

#By BMI
unique_bmi <- unique(dge.filt$samples$bmi)
sorted_unique_bmi <- sort(unique_bmi)
bmi_colors <- color_range(length(sorted_unique_bmi))
bmi_color_map <- bmi_colors[match(dge.filt$samples$bmi, sorted_unique_bmi)]
# Plot MDS with BMI color gradient
plotMDS(dge.filt, col = bmi_color_map)
title(main = "By BMI")
# Legend for BMI color gradient
legend("topright", legend = sorted_unique_bmi, fill = bmi_colors, cex = 0.6)


## ----clustering, fig.height=10, fig.width=15, dpi=100, out.width="100%", echo=FALSE, fig.cap="Hierarchical clustering of the samples. Labels correspond to DiseaseState_sampleID and colors indicate the different potential confounding variables."----
par(mfrow=c(2, 2))

#Sex
sex_colors <- c("red", "blue")  
sex_state_colors <- sex_colors[as.numeric(factor(dge.filt$samples$sex))]
sampleDendrogram_sex <- hclust(d, method = "complete")  
sampleDendrogram_sex <- as.dendrogram(sampleDendrogram_sex)
batch <- sex_state_colors
names(batch) <- se$title
sampleDendrogram_sex <- dendrapply(sampleDendrogram_sex,
                                   function(x, batch, labels) {
                                     if (is.leaf(x)) {
                                       attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")]))
                                       attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                     }
                                     x
                                   }, batch, outcome)
plot(sampleDendrogram_sex, main = "Colored by sex", cex = 0.6)
legend("topright", c("female", "male"), fill = c("red", "blue"), cex=0.6)

#Death_cause
death_cause_colors <- c("red", "blue", "darkgreen")  
death_cause_state_colors <- death_cause_colors[as.numeric(factor(dge.filt$samples$cause_death))]
sampleDendrogram_death_cause <- hclust(d, method = "complete")  
sampleDendrogram_death_cause <- as.dendrogram(sampleDendrogram_death_cause)
batch <- death_cause_state_colors
names(batch) <- se$title
sampleDendrogram_death_cause <- dendrapply(sampleDendrogram_death_cause,
                                           function(x, batch, labels) {
                                             if (is.leaf(x)) {
                                               attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")]))
                                               attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                             }
                                             x
                                           }, batch, outcome)
plot(sampleDendrogram_death_cause, main = "Colored by cause of death", cex = 0.6)
legend("topright", c("Anoxia", "Cerebral Edema", "Head Trauma"), fill = c("red", "blue", "darkgreen"), cex=0.6)

#Age
# Age
logCPM <- cpm(dge.filt, log=TRUE, prior.count=3)
d <- as.dist(1 - cor(logCPM, method = "spearman"))
sampleClustering <- hclust(d)

sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
outcome <- se$title
names(outcome) <- se$title
batch <- age_color_map
names(batch)<- se$title
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main = "Colored by age", cex = 0.7)
legend("topright", legend = sorted_unique_ages, fill = age_colors, cex = 0.7)

#BMI
logCPM <- cpm(dge.filt, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)

batch <- bmi_color_map
names(batch) <- se$title
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
outcome <- se$title
names(outcome) <- se$title
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main = "Colored by BMI", cex = 0.7)
legend("topright", legend = sorted_unique_bmi, fill = bmi_colors, cex=0.6)




## -----------------------------------------------------------------------------
# Subset for T1D vs AAb
#T1D vs AAb
se.filt_T1D_AAb <- se.filt[, se.filt$disease_state %in% c("T1D", "AAb")]
se.filt_T1D_AAb$disease_state <- droplevels(se.filt_T1D_AAb$disease_state)

#AAb vs ND
se.filt_AAb_ND <- se.filt[, se.filt$disease_state %in% c("ND", "AAb")]
se.filt_AAb_ND$disease_state <- droplevels(se.filt_AAb_ND$disease_state)

#ND vs T1D
se.filt_ND_T1D <- se.filt[, se.filt$disease_state %in% c("T1D", "ND")]
se.filt_ND_T1D$disease_state <- droplevels(se.filt_ND_T1D$disease_state)

## -----------------------------------------------------------------------------
# Subset for T1D vs AAb
mod_T1D_AAb <- model.matrix(~ disease_state, colData(se.filt_T1D_AAb))
mod0_T1D_AAb <- model.matrix(~ 1, colData(se.filt_T1D_AAb))

# Subset for AAb vs ND
mod_AAb_ND <- model.matrix(~ disease_state, colData(se.filt_AAb_ND))
mod0_AAb_ND <- model.matrix(~ 1, colData(se.filt_AAb_ND))

# Subset for ND vs T1D
mod_ND_T1D <- model.matrix(~ disease_state, colData(se.filt_ND_T1D))
mod0_ND_T1D <- model.matrix(~ 1, colData(se.filt_ND_T1D))


## ----message=FALSE------------------------------------------------------------
library(sva)

# T1D vs AAb
pv_T1D_AAb <- f.pvalue(assays(se.filt_T1D_AAb)$logCPM, mod_T1D_AAb, mod0_T1D_AAb)
sum(p.adjust(pv_T1D_AAb, method="fdr") < 0.05)
sum(p.adjust(pv_T1D_AAb, method="fdr") < 0.1)

# AAb vs ND
pv_AAb_ND <- f.pvalue(assays(se.filt_AAb_ND)$logCPM, mod_AAb_ND, mod0_AAb_ND)
sum(p.adjust(pv_AAb_ND, method="fdr") < 0.05)
sum(p.adjust(pv_AAb_ND, method="fdr") < 0.1)

# ND vs T1D
pv_ND_T1D <- f.pvalue(assays(se.filt_ND_T1D)$logCPM, mod_ND_T1D, mod0_ND_T1D)
sum(p.adjust(pv_ND_T1D, method="fdr") < 0.05)
sum(p.adjust(pv_ND_T1D, method="fdr") < 0.1)


## ----pdist, echo=FALSE, out.width="600px", fig.cap="Distribution of raw p-values for an F-test on every gene between T1D vs AAb, AAb vs ND and ND vs T1D samples."----
par(mfrow=c(1,3))

hist(pv_T1D_AAb, main="T1D vs AAb", xlab="p-value", ylab="Frequency", las=1)
hist(pv_AAb_ND, main="AAb vs ND", xlab="p-value", ylab="Frequency", las=1)
hist(pv_ND_T1D, main="ND vs T1D", xlab="p-value", ylab="Frequency", las=1)


## -----------------------------------------------------------------------------
# Subset for T1D vs AAb
mod_T1D_AAb <- model.matrix(~ disease_state + age, colData(se.filt_T1D_AAb))
mod0_T1D_AAb <- model.matrix(~ age, colData(se.filt_T1D_AAb))

# Subset for AAb vs ND
mod_AAb_ND <- model.matrix(~ disease_state + age, colData(se.filt_AAb_ND))
mod0_AAb_ND <- model.matrix(~ age, colData(se.filt_AAb_ND))

# Subset for ND vs T1D
mod_ND_T1D <- model.matrix(~ disease_state + age, colData(se.filt_ND_T1D))
mod0_ND_T1D <- model.matrix(~ age, colData(se.filt_ND_T1D))

## -----------------------------------------------------------------------------
library(sva)
sv1 <- sva(assays(se.filt_T1D_AAb)$logCPM, mod_T1D_AAb, mod0_T1D_AAb)
names(sv1)

sv2 <- sva(assays(se.filt_AAb_ND)$logCPM, mod_AAb_ND, mod0_AAb_ND)
names(sv2)

sv3 <- sva(assays(se.filt_ND_T1D)$logCPM, mod_ND_T1D, mod0_ND_T1D)
names(sv3)

#For SV1
modSV1 <- cbind(mod_T1D_AAb, sv1$sv)
mod0SV1 <- cbind(mod0_T1D_AAb, sv1$sv)
pvalsSV1 <- f.pvalue(assays(se.filt_T1D_AAb)$logCPM, modSV1, mod0SV1)
sum(p.adjust(pvalsSV1, method="fdr") < 0.05)

#For SV2
modSV2 <- cbind(mod_AAb_ND, sv2$sv)
mod0SV2 <- cbind(mod0_AAb_ND, sv2$sv)
pvalsSV2 <- f.pvalue(assays(se.filt_AAb_ND)$logCPM, modSV2, mod0SV2)
sum(p.adjust(pvalsSV2, method="BH") < 0.05)

#For SV3
modSV3 <- cbind(mod_ND_T1D, sv3$sv)
mod0SV3 <- cbind(mod0_ND_T1D, sv3$sv)
pvalsSV3 <- f.pvalue(assays(se.filt_ND_T1D)$logCPM, modSV3, mod0SV3)
sum(p.adjust(pvalsSV3, method="BH") < 0.05)

## ----pdist2, echo=FALSE, out.width="600px", fig.cap="Distribution of raw p-values for an F-test on every gene between T1D vs AAb, AAb vs ND and ND vs T1D samples."----
par(mfrow=c(1,3))

hist(pvalsSV1, main="T1D vs AAb (SV adjusted)", xlab="p-value", ylab="Frequency", las=1)
hist(pvalsSV2, main="AAb vs ND (SV adjusted)", xlab="p-value", ylab="Frequency", las=1)
hist(pvalsSV3, main="ND vs T1D (SV adjusted)", xlab="p-value", ylab="Frequency", las=1)


## ----message=FALSE, echo=FALSE------------------------------------------------
# Contar cuántos elementos son TRUE en mask_T1D_AAb
mask_ND_T1D <- p.adjust(pv_ND_T1D, method="fdr") < 0.1
DEgenesEGs_ND_T1D <- names(pv_ND_T1D)[mask_ND_T1D]
DEgenesSyms_ND_T1D <- mcols(se.filt_ND_T1D)[DEgenesEGs_ND_T1D, "symbol"]
DEgenesPvalue_ND_T1D <- pv_ND_T1D[mask_ND_T1D]
DEgenesDesc_ND_T1D <- mcols(se.filt_ND_T1D)[DEgenesEGs_ND_T1D, "description"]
DEgenesDesc_ND_T1D <- sub(" \\[.+\\]", "", DEgenesDesc_ND_T1D)
DEgenesTab_ND_T1D <- data.frame(EntrezID=DEgenesEGs_ND_T1D,
                                Symbol=DEgenesSyms_ND_T1D,
                                Description=DEgenesDesc_ND_T1D,
                                "P value"=DEgenesPvalue_ND_T1D,
                                stringsAsFactors=FALSE, check.names=FALSE)
DEgenesTab_ND_T1D <- DEgenesTab_ND_T1D[order(DEgenesTab_ND_T1D[["P value"]]), ] ## order by p-value
rownames(DEgenesTab_ND_T1D) <- 1:nrow(DEgenesTab_ND_T1D)

## ----echo=FALSE---------------------------------------------------------------
## Generate full table in a CSV file and store it in the 'doc' directory
## twice, once in 'doc' to enable quick lookup during vignette editing
## and building with 'devtools::build_vignettes()' and a second time in
## 'inst/doc' to make these files available at install.
fnameCSV_ND_T1D <- "DEgenesND_T1D.csv"
fpathCSV_ND_T1D <- file.path(path2pkg, "doc", fnameCSV_ND_T1D)
write.csv(DEgenesTab_ND_T1D, fpathCSV_ND_T1D, row.names=FALSE)
fpathCSV_ND_T1D <- file.path(path2pkg, "inst", "doc", fnameCSV_ND_T1D)
write.csv(DEgenesTab_ND_T1D, fpathCSV_ND_T1D, row.names=FALSE)

## Generate full table in HTML and store it into the 'doc' directory
## twice, just as we did with the CSV file. Note that because the
## table caption is not translated from Markdown, but directly copied
## into HTML, we need to avoid using the '<' symbol, as in FDR < 10%,
## and put its HTML code instead (&lt;)

ktab_ND_T1D <- kable(DEgenesTab_ND_T1D, "html", escape=FALSE, row.names=TRUE,
                      caption=sprintf("Differentially expressed genes between ND and T1D. DE genes with FDR &lt; 10%% (CSV <a href=\"%s\" download>file</a>).",
                                      fnameCSV_ND_T1D))
ktab_ND_T1D <- kable_styling(ktab_ND_T1D,
                              bootstrap_options=c("stripped", "hover", "responsive"),
                              fixed_thead=TRUE)
fnameHTML_ND_T1D <- "DEgenesND_T1D.html"
fpathHTML_ND_T1D <- file.path(path2pkg, "doc", fnameHTML_ND_T1D)
save_kable(ktab_ND_T1D, file=fpathHTML_ND_T1D, self_contained=TRUE)
fpathHTML_ND_T1D <- file.path(path2pkg, "inst", "doc", fnameHTML_ND_T1D)
save_kable(ktab_ND_T1D, file=fpathHTML_ND_T1D, self_contained=TRUE)

## ----tableDE, message=FALSE, echo=FALSE---------------------------------------
ktab <- kable(DEgenesTab_ND_T1D[1:10, ], format="html", escape=FALSE, row.names=TRUE,
              caption=sprintf("Top-10 differentially expressed genes with lowest p-value between ND and T1D with FDR &lt; 10%%. To see the full list of DE genes, follow this <a href=\"%s\" target=\"_blank\">link</a> or download this CSV <a href=\"%s\" download>file</a>.",
                              fnameHTML_ND_T1D, fnameCSV_ND_T1D), booktabs=TRUE, ref.label="DEgenesND_T1D")
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
# T1D vs AAb
se.filt_T1D_AAb$disease_state <- relevel(se.filt_T1D_AAb$disease_state, ref="T1D")
mod1 <- model.matrix(~ disease_state + age, data=colData(se.filt_T1D_AAb))
mod0.1 <- model.matrix(~age, data=colData(se.filt_T1D_AAb))
sv1 <- sva(assays(se.filt_T1D_AAb)$logCPM, mod=mod1, mod0=mod0.1)
mod1 <- cbind(mod1, sv1$sv)
colnames(mod1) <- c(colnames(mod1)[1:3], paste0("SV", 1:sv1$n))
head(mod1, 3)
head(rowData(se.filt_AAb_ND))
#Fit the model
fit1 <- lmFit(assays(se.filt_T1D_AAb)$logCPM, mod1)
#Calculate moderate t statistics
fit1 <- eBayes(fit1, trend=TRUE)
#Add gene metadata
genesmd <- data.frame(chr=as.character(seqnames(rowRanges(se.filt_T1D_AAb))),
                       symbol=rowData(se.filt_T1D_AAb)[, 5],
                       stringsAsFactors=FALSE)

fit1$genes <- genesmd
tt1 <- topTable(fit1, coef=2, n=Inf)
DEgenes1 <- rownames(tt1)[tt1$adj.P.Val < 0.1]
length(DEgenes1)

## -----------------------------------------------------------------------------
# AAb vs ND
se.filt_AAb_ND$disease_state <- relevel(se.filt_AAb_ND$disease_state, ref = "AAb")
mod2 <- model.matrix(~ disease_state + age, data = colData(se.filt_AAb_ND))
mod0.2 <- model.matrix(~ age, data = colData(se.filt_AAb_ND))
sv2 <- sva(assays(se.filt_AAb_ND)$logCPM, mod = mod2, mod0 = mod0.2)
mod2 <- cbind(mod2, sv2$sv)
colnames(mod2) <- c(colnames(mod2)[1:3], paste0("SV", 1:sv2$n))
head(mod2, 3)

# Fit the model
fit2 <- lmFit(assays(se.filt_AAb_ND)$logCPM, mod2)
# Calculate moderated t statistics
fit2 <- eBayes(fit2, trend = TRUE)
# Add gene metadata
genesmd2 <- data.frame(chr = as.character(seqnames(rowRanges(se.filt_AAb_ND))),
                       symbol = rowData(se.filt_AAb_ND)[, 5],
                       stringsAsFactors = FALSE)
fit2$genes <- genesmd2
tt2 <- topTable(fit2, coef = 2, n = Inf)
DEgenes2 <- rownames(tt2)[tt2$adj.P.Val < 0.1]
length(DEgenes2)

## -----------------------------------------------------------------------------
# ND vs T1D
se.filt_ND_T1D$disease_state <- relevel(se.filt_ND_T1D$disease_state, ref = "ND")
mod3 <- model.matrix(~ disease_state + age, data = colData(se.filt_ND_T1D))
mod0.3 <- model.matrix(~ age, data = colData(se.filt_ND_T1D))
sv3 <- sva(assays(se.filt_ND_T1D)$logCPM, mod = mod3, mod0 = mod0.3)
mod3 <- cbind(mod3, sv3$sv)
colnames(mod3) <- c(colnames(mod3)[1:3], paste0("SV", 1:sv3$n))
head(mod3, 3)

# Fit the model
fit3 <- lmFit(assays(se.filt_ND_T1D)$logCPM, mod3)
# Calculate moderated t statistics
fit3 <- eBayes(fit3, trend = TRUE)
# Add gene metadata
genesmd3 <- data.frame(chr = as.character(seqnames(rowRanges(se.filt_ND_T1D))),
                       symbol = rowData(se.filt_ND_T1D)[, 5],
                       stringsAsFactors = FALSE)
fit3$genes <- genesmd3
tt3 <- topTable(fit3, coef = 2, n = Inf)
DEgenes3 <- rownames(tt3)[tt3$adj.P.Val < 0.1]
length(DEgenes3)

## ----pdist3, echo=FALSE, out.width="600px", fig.cap="Distribution of raw p-values and qq plots of the moderated t statistics for the test on every gene between T1D vs AAb, AAb vs ND and ND vs T1D samples with limma-trend adjusting for known and unknown variables."----

# Set up the layout for six plots in a 2x3 grid
par(mfrow = c(3, 2), mar = c(5, 5, 2, 2))

# Histogram of raw p-values for T1D vs AAb
hist(tt1$P.Value, xlab = "Raw P-values", main = "T1D vs AAb", las = 1, col = "lightblue", border = "black")

# QQ plot for T1D vs AAb
qqt(fit1$t[, 2], df = fit1$df.prior + fit1$df.residual, main = "QQ Plot T1D vs AAb", pch = ".", cex = 3)
abline(0, 1, lwd = 2)

# Histogram of raw p-values for AAb vs ND
hist(tt2$P.Value, xlab = "Raw P-values", main = "AAb vs ND", las = 1, col = "lightgreen", border = "black")

# QQ plot for AAb vs ND
qqt(fit2$t[, 2], df = fit2$df.prior + fit2$df.residual, main = "QQ Plot AAb vs ND", pch = ".", cex = 3)
abline(0, 1, lwd = 2)

# Histogram of raw p-values for ND vs T1D
hist(tt3$P.Value, xlab = "Raw P-values", main = "ND vs T1D", las = 1, col = "lightcoral", border = "black")

# QQ plot for ND vs T1D
qqt(fit3$t[, 2], df = fit3$df.prior + fit3$df.residual, main = "QQ Plot ND vs T1D", pch = ".", cex = 3)
abline(0, 1, lwd = 2)



## ----volcano,message=FALSE, echo=FALSE, message=FALSE, fig.cap="Volcano Plots showing DE genes, with upregulated and downregulated genes in the right and left respectively.", fig.width = 14, fig.height = 7----
library(gridExtra)
library(ggplot2)
library(limma)

# Define thresholds
logFC_threshold <- 1
adj.P.Val_threshold <- 0.1

# Function to create volcano plot
create_volcano_plot <- function(data, title) {
  data$threshold <- as.factor(
    ifelse(data$adj.P.Val < adj.P.Val_threshold & abs(data$logFC) >= logFC_threshold,
           ifelse(data$logFC > logFC_threshold, "Upregulated", "Downregulated"),
           "Not Significant")
  )
  data <- data[order(data$adj.P.Val), ]
  top10_genes <- head(data, 10)

  volcano_plot <- ggplot(data, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
    geom_point(alpha = 0.7, size = 1.5) +
    theme_minimal() +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    labs(title = title,
         x = "Log Fold Change",
         y = "-Log10 Adjusted P-Value") +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(adj.P.Val_threshold), col = "black", linetype = "dashed") +
    geom_text(data = top10_genes, aes(label = symbol), vjust = 1, size = 4)

  return(volcano_plot)
}

plot1 <- create_volcano_plot(tt1, "AAb vs T1D volcano plot")
plot2 <- create_volcano_plot(tt2, "AAb vs ND volcano plot")
plot3 <- create_volcano_plot(tt3, "ND vs T1D volcano plot")

# Arrange the plots in one row
grid.arrange(plot1, plot2, plot3, nrow = 2)

## ----GO, message=FALSE--------------------------------------------------------
library(org.Hs.eg.db)
library(GOstats)

# Define geneUniverse as all genes
geneUniverse <- rownames(se)

# ND vs T1D (tt3)

# Build parameters object
params3 <- new("GOHyperGParams", geneIds=DEgenes3,    
                 universeGeneIds=geneUniverse,
                 annotation="org.Hs.eg.db", ontology="BP",
                 pvalueCutoff=0.05, testDirection="over")

#Run functional enrichment
hgOver3 <- hyperGTest(params3)
hgOver3

# Results
resHgOver3<- summary(hgOver3)
dim(resHgOver3)
head(resHgOver3)

## ----pdistGO, echo=FALSE, out.width="600px", fig.cap="Distribution of p-values of the GO enrichment analysis."----

hist(pvalues(hgOver3), xlab="P-values", main="ND vs T1D", las=1, col= "lightgreen", border="black")


## ----GO cond, message=FALSE---------------------------------------------------
# ND vs T1D (tt3)
conditional(params3) <- TRUE
hgOverCond3 <- hyperGTest(params3)
hgOverCond3

#into a df
resHgOverCond3 <- summary(hgOverCond3)
dim(resHgOverCond3)
head(resHgOverCond3)

#Filter with a minimum number of count and size
mask <- resHgOverCond3$Size >= 3 & resHgOverCond3$Size <= 300 & + resHgOverCond3$Count >= 3
resHgOverCond3 <- resHgOverCond3[mask, ]
dim(resHgOverCond3)
ord <- order(resHgOverCond3$OddsRatio, decreasing=TRUE)
resHgOverCond3 <- resHgOverCond3[ord, ]

#Extratc genes of each GO term
geneIDs <- geneIdsByCategory(hgOverCond3)[resHgOverCond3$GOBPID]
geneSYMs <- sapply(geneIDs, function(id) mapIds(org.Hs.eg.db, keys = id, column = "SYMBOL", keytype = "ENTREZID"))
geneSYMs <- sapply(geneSYMs, paste, collapse=", ")
goresults <- cbind(resHgOverCond3, Genes=geneSYMs)
rownames(goresults) <- 1:nrow(goresults)
head(goresults)

ktab <- kable(goresults, "html", caption="GO results.") %>% kable_styling(bootstrap_options = c("striped", "hover", "responsive"), fixed_thead = TRUE)
fnameHTML_ND_D <- "ND_T1D_go.html"
fpathHTML_ND_D <- file.path(path2pkg, "doc", fnameHTML_ND_D)
save_kable(ktab, file = fpathHTML_ND_D, self_contained = TRUE)


## ----pdistGO2, echo=FALSE, out.width="600px", fig.cap="Distribution of p-values of the GO enrichment analysis after conditionaltest."----

# Set up the layout for six plots in a 2x3 grid
par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))

# Histogram of p-values before filtering for ND vs T1D
hist(pvalues(hgOverCond3), xlab = "P-values", main = "Conditional test before filtering (ND vs T1D)", las = 1, col = "lightgreen", border = "black", cex.main=0.8)

# Histogram of p-values after filtering for ND vs T1D
pvals_ndt1d <- pvalues(hgOverCond3)
hist(pvals_ndt1d[resHgOverCond3$GOBPID], xlab = "P-values", xlim = c(0, 1),
     main = "Conditional test after filtering (ND vs T1D)", las = 1, col = "lightgreen", border = "black", cex.main=0.8)

## ----GOdotplot, echo=FALSE, message=FALSE, out.width="600px", fig.cap="Dotplot of the top 25 GO enriched terms in T1D vs ND comparison."----
# ND vs T1D (tt3)
# Sort the data frame based on P-value and select the top
library(dplyr)
top_terms3 <- resHgOverCond3 %>%
  arrange(Pvalue) %>%
  head(25)

# Create the dot plot
ggplot(top_terms3, aes(x = OddsRatio, y = reorder(Term, OddsRatio), size = Count, color = Pvalue)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "darkblue", high = "red") +
  labs(x = "Odds Ratio", y = "Term", size = "Count", color = "P-value") +
  ggtitle("Top 25 GO Terms with Lowest P-values
                        ND vs T1D") +
  theme_minimal()


## ----GSEAdata, message=FALSE--------------------------------------------------
library(GSVAdata)

# C2 Curated Gene Sets
data(c2BroadSets)
c2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
               grep("^REACTOME", names(c2BroadSets)), grep("^BIOCARTA", names(c2BroadSets)))]

gsc<-GeneSetCollection(c2BroadSets)

## ----GSEA Analysis, message=FALSE---------------------------------------------
library(fgsea)

gsets <- geneIds(gsc)

#AAb vs T1D
stats1 <- tt1$t
names(stats1) <- rownames(tt1)
fgseares1 <- fgsea(gsets, stats1, minSize=5, maxSize=300)

#AAb vs ND
stats2 <- tt2$t
names(stats2) <- rownames(tt2)
fgseares2 <- fgsea(gsets, stats2, minSize=5, maxSize=300)

#T1D vs ND
stats3 <- tt3$t
names(stats3) <- rownames(tt3)
fgseares3 <- fgsea(gsets, stats3, minSize=5, maxSize=300)


## ----GSEAfilter1,echo=FALSE, message=FALSE------------------------------------
#AAb vs T1D
gsea1 <- fgseares1[fgseares1$padj < 0.01, ]

gsea1$log2err <- NULL
colnames(gsea1) <- c("Pathway", "PValue", "PadjValue", "ES", "NES","Size", "LeadingEdge")

ktab1 <- kable(gsea1[order(gsea1$PadjValue), ], caption="Enriched Pathways between AAb and T1D samples", format = "html", booktabs = TRUE, row.names = FALSE) 
kable_styling(ktab1, position="center")

## ----GSEAdotplot1, echo=FALSE, message=FALSE, out.width="600px", fig.cap="Dotplot of the top 25 GSEA enriched pathways between AAb vs T1D."----
# Select the top 25 pathways
gsea1_top <- head(gsea1, 25)

# Create the dotplot
ggplot(gsea1_top, aes(x = NES, y = reorder(Pathway, NES), size = Size, color = PValue)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "darkblue", high = "red") +
  labs(x = "NES", y = "Pathway", size = "Size", color = "P-value") +
  ggtitle("AAb vs T1D") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) 

## ----GSEAfilter2,echo=FALSE,message=FALSE-------------------------------------
#AAb vs ND
gsea2 <- fgseares2[fgseares2$padj < 0.01, ]

gsea2$log2err <- NULL
colnames(gsea2) <- c("Pathway", "PValue", "PadjValue", "ES", "NES","Size", "LeadingEdge")

ktab2 <- kable(gsea2[order(gsea2$PadjValue), ], caption="Enriched Pathways between AAb and ND samples", format = "html", booktabs = TRUE, row.names = FALSE)
kable_styling(ktab2, position="center")

## ----GSEAdotplot2, echo=FALSE, message=FALSE, out.width="600px", fig.cap="Dotplot of the top 25 GSEA enriched pathways between AAb and ND."----
# Select the top 25 pathways
gsea2_top <- head(gsea2, 25)

# Create the dotplot
ggplot(gsea2_top, aes(x = NES, y = reorder(Pathway, NES), size = Size, color = PValue)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "darkblue", high = "red") +
  labs(x = "NES", y = "Pathway", size = "Size", color = "P-value") +
  ggtitle("AAb vs ND") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) 

## ----GSEAfilter3,echo=FALSE, message=FALSE------------------------------------
#T1D vs ND
gsea3 <- fgseares3[fgseares3$padj < 0.01, ]

gsea3$log2err <- NULL
colnames(gsea3) <- c("Pathway", "PValue", "PadjValue", "ES", "NES","Size", "LeadingEdge")

ktab3 <- kable(gsea3[order(gsea3$PadjValue), ], caption="Enriched Pathways between T1D and ND samples", format = "html", booktabs = TRUE, row.names = FALSE)
kable_styling(ktab3, position="center")


## ----GSEAdotplot3, echo=FALSE, message=FALSE, out.width="600px", fig.cap="Dotplot of the top 25 GSEA enriched pathways between T1D and ND."----
# Select the top 25 pathways
gsea3_top <- head(gsea3, 25)

# Create the dotplot
ggplot(gsea3_top, aes(x = NES, y = reorder(Pathway, NES), size = Size, color = PValue)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "darkblue", high = "red") +
  labs(x = "NES", y = "Pathway", size = "Size", color = "P-value") +
  ggtitle("T1D vs ND") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))  

## -----------------------------------------------------------------------------
sessionInfo()

