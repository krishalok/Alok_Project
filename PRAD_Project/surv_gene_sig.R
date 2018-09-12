# Lia Harrington 2017 - GBM Gene Signatures
# scripts/surv_gene_sig.R
#
# Usage:
# Run in command line:
#
#       Rscript scripts/surv_gene_sig.R
#
# Output:
# Produces survival curves and analyses for cell TCGA GBM data

# Set Dependencies --------------------------------------------------------

library(survival)
library(gplots)

# Prepare Expression Data -------------------------------------------------

data_file = file.path("data", "data_RNA_Seq_v2_mRNA_median_Zscores.txt")
data_file = file.path("data", "genomicMatrix")

expression <- read.table(data_file, sep = "\t",header = T, row.names = 1)[,-1]
expression <- read.table(data_file, sep = "\t",header = T, row.names = 1)
#expression<-subset(expression,expression_T$Gene_Symbol %in% row.names(expression))
#transpose dataset so each patient a row
expression <- data.frame(t(expression))

sig_file<-file.path("data", "jha et al.txt")
expression_T<-read.table(sig_file, sep = "\t", 
                        header = T)

# just grab genes in signature 
expression2 <- data.frame(expression$IL1RAPL2,expression$PHOX2A,expression$EIF2B2,
                          expression$FOXS1,expression$PIGW,expression$ALPL,
                          expression$GPR4,expression$FHOD1,expression$EVI2B,
                          expression$CCL27,expression$ALPK1,expression$TOMM34,
                          expression$B3GNT9,expression$RAX2,expression$S100A3,
                          expression$AZI1,expression$PDLIM4,expression$MEIS3,
                          expression$NCRNA00115,expression$YY2,expression$CCDC39,
                          expression$RACGAP1P,expression$FAM76A,expression$C18orf10,
                          expression$KLHL30,expression$PRPF39,expression$FBXW8,
                          expression$SLC15A3,expression$SLC25A12,expression$ALG9,
                          expression$TEX101,expression$EFCAB10,expression$TMEM175,
                          expression$SLC16A14,expression$CARNS1,expression$CMPK2,
                          expression$RNGTT,expression$TMEM26,expression$PIK3CD,
                          expression$ARL11,expression$LOC283663,expression$ANXA10,
                          expression$ZNF611,expression$PLAG1,expression$FAT4,
                          expression$PSPC1,expression$CHIC2,expression$GUCY2D,
                          expression$SLC27A3,expression$LPAR1,expression$TBX15,
                          expression$TOP3A,expression$FOXR1,expression$C5orf25,
                          expression$KLHDC8B,expression$BCL10,expression$RBMY2EP,
                          expression$OVCA2,expression$KIAA0355,expression$FUT1,
                          expression$BBS12,expression$ST7,expression$ARHGAP11A,
                          expression$MTERF,expression$THBS3,expression$CIITA,
                          expression$SLC35F4,expression$ABCB4,expression$GJD3,
                          expression$CD86,expression$EFNA4,expression$ANTXRL,
                          expression$ANKRD2,expression$FOLR3,expression$LRRC70,
                          expression$MAL,expression$SIKE1,expression$OPN3,expression$DAPP1,
                          expression$FAM115A,expression$SAMD11,expression$RAB11FIP2,
                          expression$KIF18B,expression$ZNF564,expression$ZNF404,
                          expression$NFATC4,expression$IL16,expression$CBR4,expression$MS4A7,
                          expression$ROPN1,expression$NEDD1,expression$TWIST2,expression$GDF7,
                          expression$SLC41A2,expression$ENTPD2,expression$ASPM,expression$APBB1,
                          expression$FTO,expression$LPA,expression$MRAS,expression$BTLA,
                          expression$FOXD3,expression$ICAM4,expression$LOC152217,expression$DCAF4L2,
                          expression$PRLHR,expression$UGT1A6,expression$UCK2,expression$SRPK3,
                          expression$IL1F8,expression$CHI3L1,expression$ORMDL1,expression$FAM55B,
                          expression$PYROXD1,expression$METTL13,expression$HOXB9,expression$ZNF507,
                          expression$HOXB4,expression$KLF1,expression$KIAA0415,expression$TLL1,
                          expression$BTBD16,expression$DUSP8,expression$AP3M2,expression$ALOX5AP,
                          expression$ACER1,expression$MYF6,expression$SPRED1,expression$FGD1,
                          expression$PXMP2,expression$KDELR3,expression$JRK,expression$DAPK2,
                          expression$RCAN1,expression$OR8A1,expression$RPL10L,expression$NCEH1,
                          expression$PTPLB,expression$CHRNB1,expression$CYP7A1,expression$RIPPLY2)
# normalize expression
expression2 <- data.frame(scale(expression2))

new_row_names <- c()
current_names <- row.names(expression)
for (r in 1:length(current_names)){
  name <- current_names[r]
  new_row_names = c(new_row_names, paste("TCGA-", substr(name, 6,7), "-", 
                                         substr(name, 9, 12), sep = ""))
}

expression2$PATIENT_ID <- new_row_names

# average over those with multiple measurements 
expression_ag <- aggregate(.~PATIENT_ID, FUN=mean, data=expression2)

# Prepare Clinical Data ---------------------------------------------------

data_file = file.path("data", "data_bcr_clinical_data_patient.txt")
clinical <- read.table(data_file, header = T, 
                       row.names = 1, sep = "\t")

small_clinical <- data.frame(clinical$PATIENT_ID, clinical$OS_MONTHS, 
                             clinical$OS_STATUS)

small_clinical <- small_clinical[small_clinical$clinical.OS_STATUS != "[Not Available]", ]

# create binary variable if dead or alive 
Dead = rep(NA, dim(small_clinical)[1])
for (i in 1:length(small_clinical$clinical.OS_STATUS)){
  if (small_clinical$clinical.OS_STATUS[i] == "DECEASED"){
    Dead[i] = 1
  }
  else if (small_clinical$clinical.OS_STATUS[i] == "LIVING"){
    Dead[i] = 0
  }
}

small_clinical$Dead = Dead
small_clinical$clinical.OS_MONTHS <- as.integer(small_clinical$clinical.OS_MONTHS)

# Merge Data --------------------------------------------------------------

merged_data <- merge(small_clinical, expression_ag, by.x = "clinical.PATIENT_ID", 
                     by.y = "PATIENT_ID")

head(merged_data)

# Survival Analysis -------------------------------------------------------

coxfit <- coxph(Surv(clinical.OS_MONTHS, Dead) ~., data=merged_data[, c(2, 4:146)])

cox.zph(coxfit)

coefs <- as.matrix(coxfit$coefficients)

risk_score <- as.matrix((merged_data[, 5:146] - colMeans(merged_data[, 5:146]))) %*% coefs 

merged_data$risk_score <- risk_score

group <- rep(NA, dim(merged_data)[1])

for (i in 1:dim(merged_data)[1]){
  if (risk_score[i] > 0){
    group[i] = "H"
  }
  else if (risk_score[i] < 0){
    group[i] = "L"
  }
}

merged_data$group <- group

my_surv <- survfit(Surv(clinical.OS_MONTHS, Dead) ~ group, data = merged_data)
par(mfrow=c(1,1))
plot(my_surv, conf.int=F, lty = 1, col = c("red", "black"), 
     xlab = "Months of Overall Survival", ylab = "Percent Survival")
legend("bottomleft", legend=c("H", "L"), col=c("red","black"), lty=1, 
       horiz=F, bty='n')
png('gene_signature_survival.png', res = 300)
dev.off()

risk_low <- merged_data$risk_score[merged_data$group == 'L']
risk_high <- merged_data$risk_score[merged_data$group == 'H']

mean(risk_low)
mean(risk_high)

t.test(risk_low, risk_high)

survdiff(Surv(clinical.OS_MONTHS, Dead) ~ group, data = merged_data)

get_95_ci <- function(values){
  n = length(values)
  s = sd(values)
  error <- qnorm(.975)*(s/sqrt(n))
  ci <- mean(values) + c(-1,1)*error
  return(ci)
}

get_95_ci(risk_low)
get_95_ci(risk_high)

coxfit2 <- coxph(Surv(clinical.OS_MONTHS, Dead) ~ group, data=merged_data)

# get hazard ratio
hz <- exp(1)^coxfit2$coefficients

percen_reduct <- (1-hz)*100
