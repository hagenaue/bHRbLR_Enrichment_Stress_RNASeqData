#Quickly examining the relationship between Angela's NACC RNA-Seq data and other variables besides enrichment/SD (hormones, behavior)
#01-29-2021

#Builds on code/workspace: 
#Angela_QuickAnalysis_NACC_20210129.R
#Angela_Workspace_NACC_20210129.RData

############


#I should probably double-check really quick which behaviors and hormones actually differ by treatment group in this small sample.
#We can output prettier/formal versions of this later

(colnames(Sample_BehavHormonalData_NACC_NoOutliers_Ordered))

boxplot(CORT~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#Yep, looks similar to the full sample - very pretty - elevated in SD, especially w/o EE

#This one is a cleaner output measure but has a dinky sample size:
boxplot(time_approaching_stimulus_animal_Video~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#But shows a clear effect of EE increasing social approach even in this small sample size.

boxplot(time_on_top_of_stimulus_animal_Video~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#Also a very pretty effect of EE

boxplot(time_social_interaction_EthoVision~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#The messier version - bigger sample size, but data is noisy.

boxplot(time_social_avoidance_Ethovision~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#also larger sample size but messy

Sample_BehavHormonalData_NACC_NoOutliers_Ordered$SUB_DefeatDay1
#all of the defeat related variables are half the sample size - so much smaller.
boxplot(SUB_DefeatDay1~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#no difference in this subgroup

Sample_BehavHormonalData_NACC_NoOutliers_Ordered$SUB_DefeatDay4
sum(is.na(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$SUB_DefeatDay4)==F)
#18
boxplot(SUB_DefeatDay4~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#there's a difference due to EE in this group in submissiveness like the larger sample, but the sample size is puny.
summary.lm(lm(SUB_DefeatDay4~Enrichment, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered))
#not sig in small sample size.


#No effect of treatment group in this small sample (or maybe not an interesting/useful effect)
boxplot(time_open_arms_EPM~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
boxplot(distance_open_field~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
boxplot(time_centre_open_field~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#Goes up w/SD like in main sample but maybe not the best measure.

#lower priority these didn't differ in the main sample in bLRs
boxplot(testosterone~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
boxplot(Oxytocin..pg.ml.~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
boxplot(IL.6..pg.ml.~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#Goes down in SD in this subsample
summary.lm(lm(IL.6..pg.ml.~Enrichment+Social_Defeat, data=Sample_BehavHormonalData_NACC_NoOutliers_Ordered))
#Not sig, at least using basic stats.


#So most likely to produce results that might be useful for interpreting in EE and SD results:
#1) Cort
#2) Other hormones, just because they were collected at the same time as the tissue and more likely to directly effect the brain. (but they don't differ by treatment group, so might be an awkward story)
#3) Video versions of time on top and time approaching
#4) EPM because it is a relatively good measure of anxiety? (but again doesn't differ by treatment group...)



#******************************


#Corticosterone:
#I tried this one with the hope that it would provide easy validation

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129/Corticosterone_Log10")

#do we want to log transform it first?

hist(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, breaks=12)
hist(log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), breaks=12)
#The variable definitely has a much more normal distribution when log transformed... but does that make sense in terms of biological impact?

TempDFForDesign<-data.frame(Sample_MetaData_NACC_NoOutliers_Ordered, VarOfInterest=log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT))

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)

#I should probably double-check whether cort differs by day of dissection:

boxplot(VarOfInterest~date_of_dissection, data=TempDFForDesign2)
#Nope, looks good. #phew


design <- model.matrix(~VarOfInterest+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=TempDFForDesign)

str(design)
#num [1:29, 1:5] 1 1 1 1 1 1 1 1 1 1 ...

# 4 variables fit to 29 datapoints is a little tight. We may want to reduce the design.
write.csv(design, "design.csv")

design

# (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# 68           1      4.432649     47.0                               0                           0
# 69           1      4.251395     67.4                               0                           0
# 70           1      4.249198    114.9                               1                           0
# 61           1      5.273464     61.3                               1                           0
# 62           1      5.424718     50.3                               1                           0
# 63           1      5.107549     92.2                               0                           0
# 64           1      5.028164    156.4                               0                           0
# 65           1      5.361917     58.0                               1                           0
# 41           1      4.829304     91.2                               0                           0
# 42           1      5.175512     40.6                               1                           1
# 43           1      5.455910     98.9                               0                           0
# 44           1      5.693111     68.8                               0                           0
# 45           1      6.026942     81.0                               0                           1
# 46           1      5.451940     97.8                               0                           1
# 47           1      5.813714     97.8                               0                           1
# 48           1      5.754272     54.2                               0                           0
# 74           1      5.742489     63.3                               0                           0
# 75           1      5.118265     50.7                               0                           0
# 76           1      5.184123     47.8                               0                           0
# 32           1      4.616686     42.7                               1                           0
# 33           1      5.064083     51.5                               0                           1
# 34           1      4.799065     56.6                               0                           0
# 35           1      5.358506     78.1                               0                           0
# 71           1      5.626853     76.6                               1                           0
# 72           1      5.639984     69.5                               1                           0
# 73           1      5.445137     55.7                               0                           0
# 37           1      5.232742     62.7                               1                           0
# 38           1      4.608098     27.7                               1                           1
# 40           1      4.507046     55.4                               0                           1
# attr(,"assign")
# [1] 0 1 2 3 4
# attr(,"contrasts")
# attr(,"contrasts")$date_of_RNA_extraction
# [1] "contr.treatment"
# 
# attr(,"contrasts")$date_of_dissection
# [1] "contr.treatment"

str(NACC_RNASeq_dge_TMM)

#This is annoying, but I'm not sure how to subset the dge object, so I think I may need to subset the counts and re-normalize each time I run an analysis on a subset of the data.
str(NACC_RNASeq_Matrix_noLowHits)

Temp_dge<-DGEList(counts=NACC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:17765, 1:29] 3 12 9824 2734 153 2089 222 381 1800 8 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:17765] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# .. .. .. ..$ : chr [1:29] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
# .. ..$ :'data.frame':	29 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:29] 29653288 53054264 32601268 30882918 28005917 ...
# .. .. ..$ norm.factors: num [1:29] 1 1 1 1 1 1 1 1 1 1 ...

Temp_NACC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_NACC_RNASeq_dge_TMM)


v <- voom(Temp_NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covConcExtractDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# -1        1130             0        0                               0                           0
# 0         4855         17765    17765                           17765                       17765
# 1        11780             0        0                               0                           0

#It doesn't seem like the covariates are doing much.


topTable(efit, coef=2)
#Some things are FDR<0.10 for CORT

topTable(efit, coef=3)
#Nothing is close to significant for RNAconc

topTable(efit, coef=4)
#there are things that matter for RNA extraction  (FDR<0.10)

topTable(efit, coef=5)
#there are things that matter for dissection (FDR<0.10)

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_RNAconc_RNAextract_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")

#Does it look valid?

#Top gene regulated by cort is Hif3a
#known to be regulated by cort (with a GR binding site)... but in the opposite direction (upregulation):
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3544487/
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6553161/
  
 #AABR07044388.2
  #This gene overlaps Cyp21a1 (a cytochrome P450 enzyme that is required to produce corticosterone): Chr 20:4,436,419 - 4,486,680 (+)
  https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=15004044
https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP21A2&keywords=Cyp21a1
  
#AABR07037536.1
  #On the x-chromosome, can't find much about it.
  
#Not sig, just looking for additional validation:
  
  #Zbtb16
  #Yep, it is a GR-response gene:
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892747/
  
  #Apold1
  #Yep, responds to acute stress too:
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5183579/
  
  #Apcdd1
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3507324/
  
  #Yep, these all seem legitimate. We're probably just underpowered.
  #... and maybe we shouldn't be log transforming???
  
  
  #How about the top SD affected genes?
  
  #RT1-N2: upregulated with cort, not sig
  #RT1-CE4: upregulated with cort, p=0.028
  #Slc26a8: upregulated with cort, p=0.01023
  #RT1-CE5: upregulated with cort, not sig
  #Tmtc1: downregulated with cort, p=0.00129

#Well that's cool. It fits the pattern of effects for SD.
  

#Example plots:


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Hif3a"),]

pdf("Hif3a_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="Hif3a Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Hif3a_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="Hif3a Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

#That relationship is also really driven by 3 animals. :(

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="AABR07044388.2"),]

pdf("AABR07044388_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="AABR07044388.2 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("AABR07044388_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="AABR07044388.2 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

#That is really pretty.

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="AABR07037536.1"),]

pdf("AABR07037536_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="AABR07037536.1 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("AABR07037536_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="AABR07037536.1 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

#Less  pretty but still convincing.


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Zbtb16"),]

pdf("Zbtb16_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="Zbtb16 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Zbtb16_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="Zbtb16 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

#That relationship is also really driven by the same 3 animals. :(


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Apold1"),]

pdf("Apold1_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="Apold1 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Apold1_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="Apold1 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

#OOh - that one is pretty.



GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="RT1-CE4"),]

pdf("RT1-CE4_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="RT1-CE4 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("RT1-CE4_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="RT1-CE4 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Slc26a8"),]

pdf("Slc26a8_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="Slc26a8 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Slc26a8_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="Slc26a8 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

summary.lm(lm(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)+as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)))

# Call:
#   lm(formula = GeneY ~ log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT) + 
#        as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.42729 -0.18675 -0.03475  0.13505  0.59849 
# 
# Coefficients:
#   Estimate
# (Intercept)                                                                     2.90305
# log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)                    0.38518
# as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)bLR EE + SD  -0.02761
# as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)bLR NIL      -0.75743
# as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)bLR NIL + SD -0.52887
# Std. Error t value
# (Intercept)                                                                       0.68544   4.235
# log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)                      0.13198   2.918
# as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)bLR EE + SD     0.16191  -0.171
# as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)bLR NIL         0.15366  -4.929
# as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)bLR NIL + SD    0.15944  -3.317
# Pr(>|t|)    
# (Intercept)                                                                     0.00029 ***
#   log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)                    0.00753 ** 
#   as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)bLR EE + SD   0.86603    
# as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)bLR NIL      4.97e-05 ***
#   as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)bLR NIL + SD  0.00289 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2908 on 24 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.6776,	Adjusted R-squared:  0.6238 
# F-statistic: 12.61 on 4 and 24 DF,  p-value: 1.153e-05


summary.lm(lm(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)+Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat))
#Just controlling for SD makes the CORT effect weaken and the SD effect go away.

summary.lm(lm(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)+Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat+Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment))
#Extremely related to enrichment, cort effect comes back, sd effect stil gone.

summary.lm(lm(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)+Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat+Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment+Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat+Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection))

# Call:
#   lm(formula = GeneY ~ log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT) + 
#        Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat + Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment + 
#        Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat + Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.40883 -0.19803 -0.01206  0.13922  0.68803 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)                                                          1.91792    0.65103   2.946
# log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)         0.44106    0.12807   3.444
# Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeatSD               0.07064    0.13338   0.530
# Sample_MetaData_NACC_NoOutliers_Ordered$EnrichmentEE                 0.63761    0.11112   5.738
# Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection17/1/2019  0.06250    0.14332   0.436
# Pr(>|t|)    
# (Intercept)                                                          0.00705 ** 
#   log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)         0.00212 ** 
#   Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeatSD               0.60123    
# Sample_MetaData_NACC_NoOutliers_Ordered$EnrichmentEE                6.51e-06 ***
#   Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection17/1/2019  0.66670    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.297 on 24 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.6637,	Adjusted R-squared:  0.6076 
# F-statistic: 11.84 on 4 and 24 DF,  p-value: 1.878e-05

#Relationship with SD seems dependent on CORT, relationship with EE is independent.


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Tmtc1"),]

pdf("Tmtc1_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="Tmtc1 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Tmtc1_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="Tmtc1 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

summary.lm(lm(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)+Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat+Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment+Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat+Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection))
Call:
  lm(formula = GeneY ~ log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT) + 
       Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat + Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment + 
       Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat + Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.26684 -0.10315 -0.03325  0.09792  0.28721 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)                                                          6.62137    0.32676  20.264
# log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)        -0.15671    0.06428  -2.438
# Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeatSD              -0.18101    0.06694  -2.704
# Sample_MetaData_NACC_NoOutliers_Ordered$EnrichmentEE                -0.01580    0.05578  -0.283
# Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection17/1/2019  0.19935    0.07194   2.771
# Pr(>|t|)    
# (Intercept)                                                           <2e-16 ***
#   log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)          0.0226 *  
#   Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeatSD                0.0124 *  
#   Sample_MetaData_NACC_NoOutliers_Ordered$EnrichmentEE                  0.7794    
# Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection17/1/2019   0.0106 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1491 on 24 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.5069,	Adjusted R-squared:  0.4248 
# F-statistic: 6.169 on 4 and 24 DF,  p-value: 0.001462

#Tmtc1 is the only one that seems convincing. Part of ER. To my knowledge not already known to respond to cort: https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=1564868


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Slc19a3"),]

pdf("Slc19a3_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="Slc19a3 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Slc19a3_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="Slc19a3 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

#That's reasonably pretty, even though it isn't significant in the dataset (but highly ranked and nominally significant)

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Htra1"),]

pdf("Htra1_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="Htra1 Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Htra1_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="Htra1 Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

#Very pretty.

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Fos"),]

pdf("Fos_VsCort_byGroup.pdf", width=6, height=6)
plot(GeneY~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, ylab="Fos Log2 Cpm", xlab="Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

pdf("Fos_Vslog10Cort_byGroup.pdf", width=6, height=6)
plot(GeneY~log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT), ylab="Fos Log2 Cpm", xlab="Log10 Corticosterone", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
dev.off()

#Also very pretty.



  #I'm going to try running it again without log transform, because I'm not sure if we should assume the biological impact of CORT follows a loglinear relationship.
  setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129/Corticosterone") 
  

TempDFForDesign<-data.frame(Sample_MetaData_NACC_NoOutliers_Ordered, VarOfInterest=Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)

#....

#
dt<-decideTests(efit)
summary(decideTests(efit))
# (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# -1        3623             0        0                               0                           0
# 0         1375         17765    17765                           17764                       17765
# 1        12767             0        0                               1                           0

#weaker


##Version using a reduced model (because it has a smaller sample size than the full study)

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129/Corticosterone_Log10_fewerDF")

TempDFForDesign<-data.frame(Sample_MetaData_NACC_NoOutliers_Ordered, VarOfInterest=log10(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT))

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)

design <- model.matrix(~VarOfInterest+date_of_dissection, data=TempDFForDesign)

str(design)
#num [1:29, 1:3] 1 1 1 1 1 1 1 1 1 1 ...

write.csv(design, "design.csv")

design

Temp_dge<-DGEList(counts=NACC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

Temp_NACC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_NACC_RNASeq_dge_TMM)


v <- voom(Temp_NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covConcExtractDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) VarOfInterest date_of_dissection17/1/2019
# -1        1197             0                           0
# 0         4770         17765                       17765
# 1        11798             0                           0

#Not much gained

topTable(efit, coef=2)
#Still only two things are FDR<0.10 for CORT - I should see if they are the same genes as before

topTable(efit, coef=3)
#there are quite a few things that matter for dissection (FDR<0.10)

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_Dissect.txt")
#Very similar to what we had before




#################

#Testosterone:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129/Testosterone")

#do we want to log transform it first?

hist(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$testosterone, breaks=12)
hist(log2(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$testosterone), breaks=12)
#The variable definitely has a much more normal distribution when log transformed... but does that make sense in terms of biological impact?

TempDFForDesign<-data.frame(Sample_MetaData_NACC_NoOutliers_Ordered, VarOfInterest=log2(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$testosterone))

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)

design <- model.matrix(~VarOfInterest+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=TempDFForDesign)

str(design)
#num [1:29, 1:5] 1 1 1 1 1 1 1 1 1 1 ...

# 4 variables fit to 29 datapoints is a little tight. We may want to reduce the design.

design

# (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# 68           1      9.422696     47.0                               0                           0
# 69           1      8.769507     67.4                               0                           0
# 70           1      9.017365    114.9                               1                           0
# 61           1     10.404077     61.3                               1                           0
# 62           1     11.419960     50.3                               1                           0
# 63           1     11.438792     92.2                               0                           0
# 64           1     11.021674    156.4                               0                           0
# 65           1     10.417853     58.0                               1                           0
# 41           1     11.692180     91.2                               0                           0
# 42           1     10.442943     40.6                               1                           1
# 43           1     10.052568     98.9                               0                           0
# 44           1     10.524542     68.8                               0                           0
# 45           1      9.688425     81.0                               0                           1
# 46           1     10.709084     97.8                               0                           1
# 47           1      9.733863     97.8                               0                           1
# 48           1     10.720244     54.2                               0                           0
# 74           1     10.399812     63.3                               0                           0
# 75           1     10.472691     50.7                               0                           0
# 76           1      9.987264     47.8                               0                           0
# 32           1      8.295539     42.7                               1                           0
# 33           1     10.433585     51.5                               0                           1
# 34           1      9.490851     56.6                               0                           0
# 35           1      9.842507     78.1                               0                           0
# 71           1      9.401306     76.6                               1                           0
# 72           1      9.519440     69.5                               1                           0
# 73           1     10.329796     55.7                               0                           0
# 37           1     10.982281     62.7                               1                           0
# 38           1      9.310840     27.7                               1                           1
# 40           1      8.694532     55.4                               0                           1
# attr(,"assign")
# [1] 0 1 2 3 4
# attr(,"contrasts")
# attr(,"contrasts")$date_of_RNA_extraction
# [1] "contr.treatment"
# 
# attr(,"contrasts")$date_of_dissection
# [1] "contr.treatment"


str(NACC_RNASeq_dge_TMM)

#This is annoying, but I'm not sure how to subset the dge object, so I think I may need to subset the counts and re-normalize each time I run an analysis on a subset of the data.
str(NACC_RNASeq_Matrix_noLowHits)

Temp_dge<-DGEList(counts=NACC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:17765, 1:29] 3 12 9824 2734 153 2089 222 381 1800 8 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:17765] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# .. .. .. ..$ : chr [1:29] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
# .. ..$ :'data.frame':	29 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:29] 29653288 53054264 32601268 30882918 28005917 ...
# .. .. ..$ norm.factors: num [1:29] 1 1 1 1 1 1 1 1 1 1 ...

Temp_NACC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_NACC_RNASeq_dge_TMM)


v <- voom(Temp_NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covConcExtractDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# -1         914             0        0                               0                           1
# 0         5223         17765    17765                           17765                       17764
# 1        11628             0        0                               0                           0

#It doesn't seem like the covariates are doing much.


topTable(efit, coef=2)
#Some things are FDR<0.10 for testosterone

topTable(efit, coef=3)
#Nothing is close to significant for RNAconc

topTable(efit, coef=4)
#Maybe there are things that matter for RNA extraction but still not close (FDR=0.28)

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_RNAconc_RNAextract_Dissect.txt")

write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")


#################


#IL6:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129/IL6")

#do we want to log transform it first?

hist(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$IL.6..pg.ml., breaks=12)
#Nope.

TempDFForDesign<-data.frame(Sample_MetaData_NACC_NoOutliers_Ordered, VarOfInterest=Sample_BehavHormonalData_NACC_NoOutliers_Ordered$IL.6..pg.ml.)

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)

#I should probably double-check whether cort differs by day of dissection:

boxplot(VarOfInterest~date_of_dissection, data=TempDFForDesign2)
#Nope, looks good. #phew


design <- model.matrix(~VarOfInterest+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=TempDFForDesign)

str(design)
# num [1:28, 1:5] 1 1 1 1 1 1 1 1 1 1 ...
# 4 variables fit to 28 datapoints is a little tight. We may want to reduce the design.

write.csv(design, "design.csv")

design


#This is annoying, but I'm not sure how to subset the dge object, so I think I may need to subset the counts and re-normalize each time I run an analysis on a subset of the data.

Temp_dge<-DGEList(counts=NACC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:17765, 1:28] 3 12 9824 2734 153 2089 222 381 1800 8 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:17765] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# .. .. .. ..$ : chr [1:28] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
# .. ..$ :'data.frame':	28 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:28] 29653288 53054264 32601268 30882918 28005917 ...
# .. .. ..$ norm.factors: num [1:28] 1 1 1 1 1 1 1 1 1 1 ...

Temp_NACC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_NACC_RNASeq_dge_TMM)


v <- voom(Temp_NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covConcExtractDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# -1        3216             0        0                               0                           0
# 0         1892         17765    17765                           17765                       17765
# 1        12657             0        0                               0                           0
#It doesn't seem like the covariates are doing much.


topTable(efit, coef=2)
#Nothing is even close to sig

topTable(efit, coef=3)
#Nothing is close to significant for RNAconc

topTable(efit, coef=4)
#Nothing is close to significant RNA extraction

topTable(efit, coef=5)
#Maybe for dissection?

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_RNAconc_RNAextract_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")

#Does it look valid?
#pretty questionable/noisy. Ccl4 is 4th on the list (which is chemokine responsive) and Cd1d1 is 5th (which is antigen presenting)


############

#Oxytocin:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129/Oxytocin")

#do we want to log transform it first?

hist(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$Oxytocin..pg.ml., breaks=12)
hist(log2(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$Oxytocin..pg.ml.), breaks=12)
#Toss-up

TempDFForDesign<-data.frame(Sample_MetaData_NACC_NoOutliers_Ordered, VarOfInterest=Sample_BehavHormonalData_NACC_NoOutliers_Ordered$Oxytocin..pg.ml.)

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)

#I should probably double-check whether cort differs by day of dissection:

boxplot(VarOfInterest~date_of_dissection, data=TempDFForDesign2)
#Maybe

summary.lm(lm(VarOfInterest~date_of_dissection, data=TempDFForDesign2))
#Nope

design <- model.matrix(~VarOfInterest+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=TempDFForDesign)

str(design)
# num [1:29, 1:5] 1 1 1 1 1 1 1 1 1 1 ...
# 4 variables fit to 28 datapoints is a little tight. We may want to reduce the design.

write.csv(design, "design.csv")

design


#This is annoying, but I'm not sure how to subset the dge object, so I think I may need to subset the counts and re-normalize each time I run an analysis on a subset of the data.

Temp_dge<-DGEList(counts=NACC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:17765, 1:29] 3 12 9824 2734 153 2089 222 381 1800 8 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:17765] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# .. .. .. ..$ : chr [1:29] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
# .. ..$ :'data.frame':	29 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:29] 29653288 53054264 32601268 30882918 28005917 ...
# .. .. ..$ norm.factors: num [1:29] 1 1 1 1 1 1 1 1 1 1 ...

Temp_NACC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_NACC_RNASeq_dge_TMM)


v <- voom(Temp_NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covConcExtractDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# -1        3455             0        0                               0                           0
# 0         1682         17765    17765                           17765                       17765
# 1        12628             0        0                               0                           0
#It doesn't seem like the covariates are doing much.


topTable(efit, coef=2)
#Nothing is even close to sig (all FDR>0.9999, ouch)

topTable(efit, coef=3)
#Nothing is close to significant for RNAconc

topTable(efit, coef=4)
#Nothing is close to significant RNA extraction

topTable(efit, coef=5)
#Maybe for dissection?

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_RNAconc_RNAextract_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")

################

#How about behavior?



setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129/EPM")

#do we want to log transform it first?

hist(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_open_arms_EPM, breaks=15)
hist(log2(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_open_arms_EPM+1), breaks=15)
#Yeah, it's pretty skewed - transformation is probably a good idea. Otherwise the results will be largely driven by 4 animals.
#except the EPM data contains quite a few zeros. I can add one. Not sure if this is the best strategy. Meh.

TempDFForDesign<-data.frame(Sample_MetaData_NACC_NoOutliers_Ordered, VarOfInterest=log2(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_open_arms_EPM+1))

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)
#'data.frame':	46 obs. of  18 variables:

design <- model.matrix(~VarOfInterest+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=TempDFForDesign)

str(design)
# num [1:46, 1:6] 1 1 1 1 1 1 1 1 1 1 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:46] "68" "69" "70" "61" ...
# ..$ : chr [1:6] "(Intercept)" "VarOfInterest" "RNA_conc" "date_of_RNA_extraction28.2.2019" ...
# - attr(*, "assign")= int [1:6] 0 1 2 3 4 4
# - attr(*, "contrasts")=List of 2
# ..$ date_of_RNA_extraction: chr "contr.treatment"
# ..$ date_of_dissection    : chr "contr.treatment"

write.csv(design, "design.csv")

design


#This is annoying, but I'm not sure how to subset the dge object, so I think I may need to subset the counts and re-normalize each time I run an analysis on a subset of the data.

Temp_dge<-DGEList(counts=NACC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:17765, 1:46] 3 12 9824 2734 153 2089 222 381 1800 8 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:17765] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# .. .. .. ..$ : chr [1:46] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
# .. ..$ :'data.frame':	46 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:46] 29653288 53054264 32601268 30882918 28005917 ...
# .. .. ..$ norm.factors: num [1:46] 1 1 1 1 1 1 1 1 1 1 ...

Temp_NACC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_NACC_RNASeq_dge_TMM)


v <- voom(Temp_NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covConcExtractDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

#     (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# -1        3732             0        0                               0                           0
# 0         1161         17765    17765                           17765                       17765
# 1        12872             0        0                               0                           0
# date_of_dissection21/2/2019
# -1                        1132
# 0                        15512
# 1                         1121

#the 12/2/2019 date of dissection looks much more different than the other two

topTable(efit, coef=2)
#Nothing is even close to sig (all FDR>0.35, ouch)

topTable(efit, coef=3)
#Many things less than FDR<0.1 for RNAconc

topTable(efit, coef=4)
#There are some things that are <0.2 for Extraction

topTable(efit, coef=5)
#Lots of things with FDR<0.15

topTable(efit, coef=6)
#... and that one dissection date matters a lot

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_RNAconc_RNAextract_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")


#######################

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129/SocialApproach")

#do we want to log transform it first?

hist(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_approaching_stimulus_animal_Video, breaks=15)
#skewed, but not in a way that log transformation is going to help - basically bimodal.

TempDFForDesign<-data.frame(Sample_MetaData_NACC_NoOutliers_Ordered, VarOfInterest=Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_approaching_stimulus_animal_Video)

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)
#'data.frame':	29 obs. of  18 variables:

design <- model.matrix(~VarOfInterest+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=TempDFForDesign)

str(design)
# num [1:29, 1:5] 1 1 1 1 1 1 1 1 1 1 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:29] "68" "69" "70" "61" ...
# ..$ : chr [1:5] "(Intercept)" "VarOfInterest" "RNA_conc" "date_of_RNA_extraction28.2.2019" ...
# - attr(*, "assign")= int [1:5] 0 1 2 3 4
# - attr(*, "contrasts")=List of 2
# ..$ date_of_RNA_extraction: chr "contr.treatment"
# ..$ date_of_dissection    : chr "contr.treatment"

write.csv(design, "design.csv")

design


#This is annoying, but I'm not sure how to subset the dge object, so I think I may need to subset the counts and re-normalize each time I run an analysis on a subset of the data.

Temp_dge<-DGEList(counts=NACC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

Temp_NACC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_NACC_RNASeq_dge_TMM)


v <- voom(Temp_NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covConcExtractDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# -1        3484             0        0                               0                           0
# 0         1530         17765    17765                           17765                       17765
# 1        12751             0        0                               0                           0

topTable(efit, coef=2)
#Nothing is even close to sig (all FDR>0.75, ouch)

topTable(efit, coef=3)
#all FDR<0.65 for RNAconc

topTable(efit, coef=4)
#There are some things that are <0.25 for Extraction

topTable(efit, coef=5)
#Lots of things with FDR<0.15

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_RNAconc_RNAextract_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")

#############

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129/TimeOnTop")

#do we want to log transform it first?

hist(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_on_top_of_stimulus_animal_Video, breaks=15)
#super skewed, but not in a way that log transformation is going to help - super floor effect.

TempDFForDesign<-data.frame(Sample_MetaData_NACC_NoOutliers_Ordered, VarOfInterest=Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_on_top_of_stimulus_animal_Video)

TempDFForDesign2<-TempDFForDesign[is.na(TempDFForDesign$VarOfInterest)==F, ]
str(TempDFForDesign2)
#'data.frame':	29 obs. of  18 variables:

design <- model.matrix(~VarOfInterest+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=TempDFForDesign)

str(design)
# num [1:29, 1:5] 1 1 1 1 1 1 1 1 1 1 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:29] "68" "69" "70" "61" ...
# ..$ : chr [1:5] "(Intercept)" "VarOfInterest" "RNA_conc" "date_of_RNA_extraction28.2.2019" ...
# - attr(*, "assign")= int [1:5] 0 1 2 3 4
# - attr(*, "contrasts")=List of 2
# ..$ date_of_RNA_extraction: chr "contr.treatment"
# ..$ date_of_dissection    : chr "contr.treatment"

write.csv(design, "design.csv")

design


#This is annoying, but I'm not sure how to subset the dge object, so I think I may need to subset the counts and re-normalize each time I run an analysis on a subset of the data.

Temp_dge<-DGEList(counts=NACC_RNASeq_Matrix_noLowHits[,is.na(TempDFForDesign$VarOfInterest)==F])
str(Temp_dge)

Temp_NACC_RNASeq_dge_TMM<-calcNormFactors(Temp_dge, method="TMM")
str(Temp_NACC_RNASeq_dge_TMM)


v <- voom(Temp_NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("VarOfInterest_covConcExtractDissect_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) VarOfInterest RNA_conc date_of_RNA_extraction28.2.2019 date_of_dissection17/1/2019
# -1        3450             0        0                               0                           0
# 0         1529         17765    17765                           17765                       17765
# 1        12786             0        0                               0                           0

topTable(efit, coef=2)
#There are three genes with FDR=0.1

topTable(efit, coef=3)
#all FDR>0.29 for RNAconc

topTable(efit, coef=4)
#There are some things that are <0.20 for Extraction

topTable(efit, coef=5)
#Lots of things with FDR<0.20

write.fit(efit, adjust="BH", file="Limma_results_Model_Var_RNAconc_RNAextract_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")

