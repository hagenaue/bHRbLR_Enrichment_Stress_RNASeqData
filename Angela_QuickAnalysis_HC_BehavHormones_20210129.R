

#Corticosterone:
#I tried this one with the hope that it would provide easy validation

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/Corticosterone_Log10")

#do we want to log transform it first?

hist(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT, breaks=12)
hist(log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT), breaks=12)
#The variable definitely has a much more normal distribution when log transformed... but does that make sense in terms of biological impact?

TempDFForDesign<-data.frame(Sample_MetaData_HC_NoOutliers_Ordered, VarOfInterest=log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT))

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)

#I should probably double-check whether cort differs by day of dissection:

boxplot(VarOfInterest~date_of_dissection, data=TempDFForDesign2)
#Nope, looks good. #phew


design <- model.matrix(~VarOfInterest+date_of_dissection, data=TempDFForDesign)

str(design)
#num [1:25, 1:4] 1 1 1 1 1 1 1 1 1 1 ...

write.csv(design, "design.csv")

design


str(HC_RNASeq_dge_TMM)

#This is annoying, but I'm not sure how to subset the dge object, so I think I may need to subset the counts and re-normalize each time I run an analysis on a subset of the data.
str(HC_RNASeq_Matrix_noLowHits)

Temp_dge<-DGEList(counts=HC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

Temp_HC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_HC_RNASeq_dge_TMM)


v <- voom(Temp_HC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) VarOfInterest date_of_dissection13/3/2019 date_of_dissection17/1/2019
# -1         638             0                        2460                           0
# 0         5369         17629                       12293                       17629
# 1        11622             0                        2876                           0

topTable(efit, coef=2)
#Everything FDR>0.22 for CORT


topTable(efit, coef=4)
#dissection 17/1/2019 is FDR>0.19

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_Dissect.txt")

#write.csv(HC_RNASeq_Annotation_noLowHits, "HC_RNASeq_Annotation_noLowHits.csv")

#Does it look valid?

#Not sig, but the top-ranked genes are super pretty:
#1 is a heat-shock protein Cryab
#2 is Plekhd1
#3 is Ucp3 (one of my favorites!)
#4 is Htra1
#8 is Dbp
#18 Olig1
#17 is a cytochrome p450 enzyme
#21 is Prodh

#13 is Arc
#Fos is 83 (validation)

#Example plots:


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Cryab"),]

pdf("Cryab_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT, ylab="Cryab Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Cryab_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT), ylab="Cryab Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()
#very pretty

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Ucp3"),]

pdf("Ucp3_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT, ylab="Ucp3 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Ucp3_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT), ylab="Ucp3 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

#big outlier

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Htra1"),]

pdf("Htra1_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT, ylab="Htra1 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Htra1_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT), ylab="Htra1 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

#Less  pretty but still convincing.


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Prodh1"),]

pdf("Prodh1_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT, ylab="Prodh1 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Prodh1_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT), ylab="Prodh1 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

#That relationship looks driven by treatment group.

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Arc"),]

pdf("Arc_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT, ylab="Arc Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Arc_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT), ylab="Arc Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

#OOh - that one is pretty.

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Fos"),]

pdf("Fos_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT, ylab="Fos Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Fos_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT), ylab="Fos Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
dev.off()


#############
#Corticosterone:

#Simplest model:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/JustCort_Log10")

#do we want to log transform it first?

hist(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT, breaks=12)
hist(log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT), breaks=12)
#The variable definitely has a much more normal distribution when log transformed... but does that make sense in terms of biological impact?

TempDFForDesign<-data.frame(Sample_MetaData_HC_NoOutliers_Ordered, VarOfInterest=log10(Sample_BehavHormonalData_HC_NoOutliers_Ordered$CORT))

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)


design <- model.matrix(~VarOfInterest, data=TempDFForDesign)

str(design)
#num [1:25, 1:4] 1 1 1 1 1 1 1 1 1 1 ...

write.csv(design, "design.csv")

design


str(HC_RNASeq_dge_TMM)

#This is annoying, but I'm not sure how to subset the dge object, so I think I may need to subset the counts and re-normalize each time I run an analysis on a subset of the data.
str(HC_RNASeq_Matrix_noLowHits)

Temp_dge<-DGEList(counts=HC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

Temp_HC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_HC_RNASeq_dge_TMM)


v <- voom(Temp_HC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) VarOfInterest
# -1         550             0
# 0         5635         17629
# 1        11444             0

topTable(efit, coef=2)
#Everything FDR>0.77 for CORT

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_Dissect.txt")

