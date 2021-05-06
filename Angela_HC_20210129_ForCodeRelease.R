#Angela's bLR social defeat and environmental enrichment RNA-Seq data from the HC
#Megan Hagenauer

###########################

#Rstudio (v.1.0.153, R v. 3.4.1)

#################


#Overview of analysis versions:

## Fan ran a simple version of the differential expression analysis within Ahub (DESeq2) that didn't control for confounds and didn't use a full EE*SD interaction model: August 2019
## My first quick analysis of Fan's preprocessed output to check for potential confounds: September 5, 2019
# Updated Jan 29 2021 for a full analysis (completed after the larger NACC dataset)

########################

#The initial preprocessing was performed using the MBNI Analysis Hub (https://ahub.mbni.org) by Dr. Fan Meng in late August 2019 (the BAM files are labeled 20190830).

#Preprocessing methodological details (from Fan):
# - Alignment algorithm (w/version & any non-default parameters): STAR 2.7.0f + Encode parameter with unique hit only
# - Genome assembly: Rnor6 
# - Algorithm for read summarization (w/version & any non-default parameters): featureCounts v1.6.4, (not 100% sure)
# - Gene annotation source (w/version): Ensembl 96

##############################################

#Libraries loaded:

library(edgeR)
library(limma)
library(affy)
library(fgsea)
library(BrainInABlender)
library(plyr)

sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.15.7
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] plyr_1.8.4                 BrainInABlender_0.0.0.9000 affy_1.54.0                Biobase_2.36.2             BiocGenerics_0.22.0       
# [6] edgeR_3.18.1               limma_3.32.5               fgsea_1.2.1                Rcpp_1.0.4.6              
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.6        zlibbioc_1.22.0       munsell_0.4.3         BiocParallel_1.10.1   colorspace_1.3-2      lattice_0.20-35       rlang_0.3.1          
# [8] fastmatch_1.1-0       tools_3.4.1           grid_3.4.1            data.table_1.10.4     gtable_0.2.0          lazyeval_0.2.0        preprocessCore_1.38.1
# [15] tibble_2.1.1          crayon_1.3.4          affyio_1.46.0         gridExtra_2.3         ggplot2_2.2.1         BiocInstaller_1.26.0  compiler_3.4.1       
# [22] pillar_1.3.1          scales_0.4.1          locfit_1.5-9.1        pkgconfig_2.0.2      


##############################################

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples")

#Reading in Fan's gene-level summary output for the hippocampus data, with 6 outliers removed based on the PC1/PC2 plot:

HC_RNASeq<-read.delim("Hippocampus_remove6samples_gene_featureCounts_counts.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)

str(HC_RNASeq)
# 'data.frame':	32883 obs. of  25 variables:
#   $ Sample_126043: int  0 0 0 0 0 0 1 0 0 0 ...
# $ Sample_126044: int  0 0 0 0 0 0 1 0 0 0 ...
# $ Sample_126045: int  0 0 0 0 0 0 0 0 0 0 ...
# $ Sample_126046: int  0 0 0 0 0 0 1 0 0 0 ...
# $ Sample_126048: int  0 0 0 0 0 0 4 0 0 0 ...
# $ Sample_126049: int  0 0 0 0 0 0 0 0 0 0 ...
# $ Sample_126050: int  0 0 0 0 0 0 1 0 0 0 ...

#This data is just basic counts, so we will need to run some transformations to make graphs.

#Double-checking whether a log transform has already been performed:
max(HC_RNASeq)
[1] 773896
#Nope.


HC_RNASeq_Matrix<-as.matrix(HC_RNASeq)
str(HC_RNASeq_Matrix)
# int [1:32883, 1:25] 0 0 0 0 0 0 1 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:25] "Sample_126043" "Sample_126044" "Sample_126045" "Sample_126046" ...

HC_RNASeq_Annotation<-row.names(HC_RNASeq)


#Counts per sample (library size):
HC_RNASeq_CountsPerSample<-apply(HC_RNASeq_Matrix, 2, sum)
head(HC_RNASeq_CountsPerSample)

write.csv(HC_RNASeq_CountsPerSample, "HC_RNASeq_CountsPerSample.csv")

pdf("Histogram_CountsPerSample.pdf", width=5, height=5)
hist(HC_RNASeq_CountsPerSample, breaks=10, xlab="Library Size: Counts per Sample", main="")
dev.off()

summary(HC_RNASeq_CountsPerSample)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1811379  9937370 19628096 22789287 34609573 46551330 

#Wow - some of those counts are really low.
#1,811,379 - less than 2 million counts!!!
#9,937,370 - and a quarter of the sample has less than 10 million. :(
#19,628,096 - the median sample is still less than 20 million. :(


#Which samples have low total counts?

HC_RNASeq_CountsPerSample[HC_RNASeq_CountsPerSample<20000000]
# Sample_126043 Sample_126044 Sample_126045 Sample_126046 Sample_126048 Sample_126049 Sample_126050 Sample_126053 Sample_126054 
# 19881617      17151268      11540367       9937370       8309574       6513551       2665616      19628096      18722011 
# Sample_126055 Sample_126056 Sample_126058 Sample_126059 Sample_126065 
# 17471307      16450483       9485002       6351533       1811379 

HC_RNASeq_CountsPerSample[HC_RNASeq_CountsPerSample<10000000]
# Sample_126046 Sample_126048 Sample_126049 Sample_126050 Sample_126058 Sample_126059 Sample_126065 
# 9937370       8309574       6513551       2665616       9485002       6351533       1811379 

#Sample_126065 is the sample with the lowest counts (less than 2 million!)


#Running some basic filtering:

table(rowSums(HC_RNASeq_Matrix==0))
# 0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20 
# 13838   821   477   380   317   274   220   229   217   209   213   216   206   244   207   223   259   292   298   339   410 
# 21    22    23    24    25 
# 438   621   925  1671  9339 
#There are 13838 rows with 0's for all subjects!  Get rid of them.
#Also: I should align my analysis with what Fan did.  He says that he filtered out any genes with average counts <1. I think he means actual counts (not CPM), because he didn't calculate CPM (or something similar) for running DE, he just used DESeq2.


#Filtering out rows of data with average counts<1
keep.exprs<-(apply(HC_RNASeq_Matrix, 1, mean)>1)
table(keep.exprs)
# keep.exprs
# FALSE  TRUE 
# 15254 17629 

HC_RNASeq_Matrix_noLowHits<-HC_RNASeq_Matrix[keep.exprs,]
dim(HC_RNASeq_Matrix_noLowHits)
#[1] 17629    25
HC_RNASeq_Annotation_noLowHits<-HC_RNASeq_Annotation[keep.exprs]

#TMM normalization to correct for differences in estimated RNA production levels (recommended in limma manual) (note: I did not due this for my first round of analysis - this is new as of 01-29-2021)
head(HC_RNASeq_Matrix_noLowHits)
# Sample_126043 Sample_126044 Sample_126045 Sample_126046 Sample_126048 Sample_126049 Sample_126050 Sample_126051
# ENSRNOG00000061316             1             1             0             1             4             0             1             7
# ENSRNOG00000029897            24            12             5             4             5             2             0            29
# ENSRNOG00000014303          4891          4248          3745          3159          2727          2067           806          7687

HC_RNASeq_Log2_Filtered <- cpm(HC_RNASeq_Matrix_noLowHits, log=TRUE)

head(HC_RNASeq_Log2_Filtered)
# Sample_126043 Sample_126044 Sample_126045 Sample_126046 Sample_126048 Sample_126049 Sample_126050 Sample_126051
# ENSRNOG00000061316    -4.0286298    -3.8514251     -6.510185     -3.163483    -1.0221672     -6.510185     -1.372776   -1.97346712
# ENSRNOG00000029897     0.2847483    -0.4927351     -1.170513     -1.273968    -0.7066825     -1.652688     -6.510185    0.02921798
# ENSRNOG00000014303     7.9427124     7.9524895      8.342284      8.312542     8.3584713      8.310041      8.240331    8.06390689

#As comparison:

#no log:
HC_cpm_Filtered <- cpm(HC_RNASeq_Matrix_noLowHits)
str(HC_cpm_Filtered)

#no log, no filter:
HC_cpm <- cpm(HC_RNASeq_Matrix)

#log, no filter:
HC_cpm_Log2 <- cpm(HC_RNASeq_Matrix, log=TRUE)

#Let's compare the distributions before and after filtering:

library(RColorBrewer)

png("HC_logCPM_Distribution_BeforeAndAfterFiltering.png")
nsamples <- ncol(HC_RNASeq_Matrix)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(HC_cpm_Log2[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(HC_cpm_Log2[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
plot(density(HC_RNASeq_Log2_Filtered [,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(HC_RNASeq_Log2_Filtered [,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
dev.off()

pdf("HC_boxplot_logCPM_afterFiltering.pdf", width=10, height=5)
boxplot(HC_RNASeq_Log2_Filtered , las=2, col=col, main="")
title(main="Unnormalised data",ylab="Log-cpm")
dev.off()

#The boxes for general distribution of cpm data across all genes look very similar across all samples - I don't think any additional normalization is necessary (e.g., quantile normalization)


################################

#TMM normalization: Scaling to correct for estimated relative RNA production levels

#I didn't originally include this step, but my analysis of the NACC data revealed that library size was correlated with the principal components of variation in the data, despite cpm normalization. Also, the limma manual recommends it.

dge<-DGEList(counts=HC_RNASeq_Matrix_noLowHits)
str(dge)

# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:17629, 1:25] 1 24 4891 1351 73 1479 137 331 1357 4 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:17629] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# .. .. .. ..$ : chr [1:25] "Sample_126043" "Sample_126044" "Sample_126045" "Sample_126046" ...
# .. ..$ :'data.frame':	25 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:25] 19880273 17150068 11539545 9936653 8309006 ...
# .. .. ..$ norm.factors: num [1:25] 1 1 1 1 1 1 1 1 1 1 ...

str(dge[[2]])

dge[[2]]$lib.size
#This is the same as HC_RNASeq_CountsPerSample - useful
#Group is currently empty - do we need it?

head(dge[[1]])
# Sample_126043 Sample_126044 Sample_126045 Sample_126046 Sample_126048 Sample_126049 Sample_126050 Sample_126051
# ENSRNOG00000061316             1             1             0             1             4             0             1             7
# ENSRNOG00000029897            24            12             5             4             5             2             0            29
# ENSRNOG00000014303          4891          4248          3745          3159          2727          2067           806          7687

HC_RNASeq_dge_TMM<-calcNormFactors(dge, method="TMM")
str(HC_RNASeq_dge_TMM)

# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:17629, 1:25] 1 24 4891 1351 73 1479 137 331 1357 4 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:17629] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# .. .. .. ..$ : chr [1:25] "Sample_126043" "Sample_126044" "Sample_126045" "Sample_126046" ...
# .. ..$ :'data.frame':	25 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:25] 19880273 17150068 11539545 9936653 8309006 ...
# .. .. ..$ norm.factors: num [1:25] 0.9 0.931 1.006 0.996 1.01 ...

# 
head(HC_RNASeq_dge_TMM[[1]])
# Sample_126043 Sample_126044 Sample_126045 Sample_126046 Sample_126048 Sample_126049 Sample_126050 Sample_126051
# ENSRNOG00000061316             1             1             0             1             4             0             1             7
# ENSRNOG00000029897            24            12             5             4             5             2             0            29
# ENSRNOG00000014303          4891          4248          3745          3159          2727          2067           806          7687

#Interesting - the object seems to store the normalization factors but not the normalized data.
#Double-checked - that is what it is supposed to be. Interesting. At what point is the data actually scaled with those normalization factors?

#Let's see if it happens during the cpm function:

HC_RNASeq_Log2_Filtered_TMM <- cpm(HC_RNASeq_dge_TMM, log=TRUE)
head(HC_RNASeq_Log2_Filtered_TMM)

# Sample_126043 Sample_126044 Sample_126045 Sample_126046 Sample_126048 Sample_126049 Sample_126050 Sample_126051
# ENSRNOG00000061316    -3.9051882    -3.7666691     -6.525408     -3.159986    -1.0368771     -6.525408     -1.404262   -1.98119055
# ENSRNOG00000029897     0.4353179    -0.3916694     -1.178936     -1.268991    -0.7213901     -1.667689     -6.525408    0.02174725
# ENSRNOG00000014303     8.0946976     8.0553010      8.334032      8.318065     8.3437731      8.295048      8.208367    8.05651985

head(HC_RNASeq_Log2_Filtered)
# Sample_126043 Sample_126044 Sample_126045 Sample_126046 Sample_126048 Sample_126049 Sample_126050 Sample_126051
# ENSRNOG00000061316    -4.0286298    -3.8514251     -6.510185     -3.163483    -1.0221672     -6.510185     -1.372776   -1.97346712
# ENSRNOG00000029897     0.2847483    -0.4927351     -1.170513     -1.273968    -0.7066825     -1.652688     -6.510185    0.02921798
# ENSRNOG00000014303     7.9427124     7.9524895      8.342284      8.312542     8.3584713      8.310041      8.240331    8.06390689
#Yep - it looks like the normalization does influence the final cpm. Good.


pdf("HC_boxplot_logCPM_afterFilteringTMM.pdf", width=10, height=5)
boxplot(HC_RNASeq_Log2_Filtered_TMM , las=2, col=col, main="")
title(main="TMM normalised data",ylab="Log-cpm")
dev.off()

#Overwriting our previous object so we don't need to change the rest of the code pipeline:
HC_RNASeq_Log2_Filtered<-HC_RNASeq_Log2_Filtered_TMM


###################################

#Getting additional annotation:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress")

HC_Annotation<-read.delim("Annotation_Fan.txt", sep="\t", header=T, stringsAsFactors = F)
str(HC_Annotation)
# 'data.frame':	17335 obs. of  4 variables:
#   $ gene_id     : chr  "ENSRNOG00000037206" "ENSRNOG00000010452" "ENSRNOG00000011762" "ENSRNOG00000020388" ...
# $ gene_symbol : chr  "Ccdc77" "LOC100363502" "Elf1" "Inpp5f" ...
# $ gene_biotype: chr  "protein_coding" "protein_coding" "protein_coding" "protein_coding" ...
# $ description : chr  "coiled-coil domain containing 77 [Source:RGD Symbol;Acc:1310710]" "cytochrome c, somatic-like [Source:RGD Symbol;Acc:2322845]" "E74 like ETS transcription factor 1 [Source:RGD Symbol;Acc:620697]" "inositol polyphosphate-5-phosphatase F [Source:RGD Symbol;Acc:1305777]" ...

#For ease, let's put the annotation in the same order as the Count data:

gene_id<-HC_RNASeq_Annotation_noLowHits

ToJoin<-data.frame(gene_id, HC_RNASeq_Log2_Filtered, stringsAsFactors = F)
str(ToJoin)

HC_RNASeq_Log2_Annotated<-join(ToJoin, HC_Annotation, by="gene_id", type="left")
str(HC_RNASeq_Log2_Annotated)

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples")
write.csv(HC_RNASeq_Log2_Annotated, "HC_RNASeq_Log2_Annotated.csv")

##################################

#Meta-data

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress")

#New version of the analysis (from 01/2021)
setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/ReDone_OpenFieldData_20201008")
Sample_BehavHormonalData<-read.csv("HRLR_EE_Stress_AllBehavData_forR_withNewCORTOxytIL6_SI_OFSDScoresFixed_FixedFormatIDs_TimeOnTop_forHC2.csv", header=T, stringsAsFactors = F)

str(Sample_BehavHormonalData)
#'data.frame':	142 obs. of  58 variables:
colnames(Sample_BehavHormonalData)
# [1] "Rat_ID"                                                                                                         
# [2] "Litter"                                                                                                         
# [3] "From.Kathryn.or.not...Kathryn.may.not.have.culled.on.P1."                                                       
# [4] "Cage..pair.housed_.if.uneven.last.three.together."                                                              
# [5] "Generation"                                                                                                     
# [6] "Line"                                                                                                           
# [7] "Treatment_group..EC.Enriched.Cage.Control..EE.Environmental.Enrichment..NIL.Standard.Housing..SD.Social.Defeat."
# [8] "Enrichment"                                                                                                     
# [9] "Social_Defeat"                                                                                                  
# [10] "distance_open_field"                                                                                            
# [11] "time_centre_open_field"                                                                                         
# [12] "time_open_arms_EPM"                                                                                             
# [13] "time_social_avoidance_Ethovision"                                                                               
# [14] "time_social_interaction_EthoVision"                                                                             
# [15] "time_other_areas_in_SI_Test_EthoVision"                                                                         
# [16] "time_on_top_of_stimulus_animal_Ethovision"                                                                      
# [17] "time_on_top_of_stimulus_animal_Video"                                                                           
# [18] "time_approaching_stimulus_animal_Video"                                                                         
# [19] "time_other_behav_in_SI_Test_Video"                                                                              
# [20] "cranky_USVs"                                                                                                    
# [21] "happy_USVs"                                                                                                     
# [22] "SD.Dates"                                                                                                       
# [23] "D1.Aggressor.ID"                                                                                                
# [24] "SD.d1.score"                                                                                                    
# [25] "SD.d1.time.caged..min."                                                                                         
# [26] "D2.Aggressor.ID"                                                                                                
# [27] "SD.d2.score"                                                                                                    
# [28] "SD.d2.time.caged..min."                                                                                         
# [29] "D3.Aggressor.ID"                                                                                                
# [30] "SD.d3.score"                                                                                                    
# [31] "SD.d3.time.caged..min."                                                                                         
# [32] "D4.Aggressor.ID"                                                                                                
# [33] "SD.d4.score"                                                                                                    
# [34] "SD.d4.time.caged..min."                                                                                         
# [35] "SUB_DefeatDay1"                                                                                                 
# [36] "AGG_DefeatDay1"                                                                                                 
# [37] "OTHER_DefeatDay1"                                                                                               
# [38] "SUB_DefeatDay2"                                                                                                 
# [39] "AGG_DefeatDay2"                                                                                                 
# [40] "OTHER_DefeatDay2"                                                                                               
# [41] "SUB_DefeatDay3"                                                                                                 
# [42] "AGG_DefeatDay3"                                                                                                 
# [43] "OTHER_DefeatDay3"                                                                                               
# [44] "SUB_DefeatDay4"                                                                                                 
# [45] "AGG_DefeatDay4"                                                                                                 
# [46] "OTHER_DefeatDay4"                                                                                               
# [47] "date_of_sac"                                                                                                    
# [48] "CORT"                                                                                                           
# [49] "testosterone"                                                                                                   
# [50] "Oxytocin..pg.ml."                                                                                               
# [51] "IL.6..pg.ml."                                                                                                   
# [52] "HC_date_of_dissection"                                                                                          
# [53] "HC_hemisphere_dissected"                                                                                        
# [54] "HC_date_of_RNA_extraction"                                                                                      
# [55] "HC_RNA_conc"                                                                                                    
# [56] "HC_ratio_260_280"                                                                                               
# [57] "HC_brain_region"                                                                                                
# [58] "HC_Sequencing_core_sample_ID"  


#Let's filter both data-frames just down to the samples that were included in the RNA-Seq experiment:

table(Sample_MetaData$brain_region)
# hippocampus        Nacc 
# 31          49 

sum(Sample_MetaData$brain_region=="hippocampus")
#[1] 31

#Ah - there are 6 more samples included in the meta data because Fan removed 6 outlier samples (which is a lot in such a small dataset... ouch...)

Sample_MetaData<-Sample_MetaData[Sample_MetaData$brain_region=="hippocampus",]
str(Sample_MetaData)
#'data.frame':	31 obs. of  14 variables:
#$ Sequencing_core_sample_ID: int  126059 126060 127120 126062 126066 126067 126068 126069 126070 126048 ...

table(Sample_BehavHormonalData$HC_brain_region)
# HC 
# 111  31 

Sample_BehavHormonalData<-Sample_BehavHormonalData[Sample_BehavHormonalData$HC_brain_region=="HC",]
str(Sample_BehavHormonalData)
#'data.frame':	31 obs. of  58 variables:

#To make things easier, I need to create a MetaData file that is in the same order as the HC RNASeq data file.
colnames(HC_RNASeq_Log2_Filtered)
# [1] "Sample_126043" "Sample_126044" "Sample_126045" "Sample_126046" "Sample_126048" "Sample_126049" "Sample_126050" "Sample_126051" "Sample_126052"
# [10] "Sample_126053" "Sample_126054" "Sample_126055" "Sample_126056" "Sample_126057" "Sample_126058" "Sample_126059" "Sample_126065" "Sample_126066"
# [19] "Sample_126067" "Sample_126068" "Sample_126069" "Sample_126070" "Sample_127119" "Sample_127120" "Sample_127121"

#the column names are just in sample ID order.

Sample_MetaData$Sequencing_core_sample_ID
# [1] 126059 126060 127120 126062 126066 126067 126068 126069 126070 126048 126049 126050 126051 126052 126053 126054 126055 126043
# [19] 126044 126045 126046 127119 126040 126041 126042 126063 127121 126065 126056 126057 126058

Sample_MetaData$SampleID<-paste("Sample_", Sample_MetaData$Sequencing_core_sample_ID, sep="")
Sample_MetaData$SampleID
# [1] "Sample_126059" "Sample_126060" "Sample_127120" "Sample_126062" "Sample_126066" "Sample_126067" "Sample_126068"
# [8] "Sample_126069" "Sample_126070" "Sample_126048" "Sample_126049" "Sample_126050" "Sample_126051" "Sample_126052"
# [15] "Sample_126053" "Sample_126054" "Sample_126055" "Sample_126043" "Sample_126044" "Sample_126045" "Sample_126046"
# [22] "Sample_127119" "Sample_126040" "Sample_126041" "Sample_126042" "Sample_126063" "Sample_127121" "Sample_126065"
# [29] "Sample_126056" "Sample_126057" "Sample_126058"

Sample_MetaData_HC_NoOutliers<-Sample_MetaData[which(Sample_MetaData$SampleID%in%colnames(HC_RNASeq_Log2_Filtered)), ]
str(Sample_MetaData_HC_NoOutliers)
# 'data.frame':	25 obs. of  15 variables:
#   $ Sequencing_core_sample_ID: int  126059 127120 126066 126067 126068 126069 126070 126048 126049 126050 ...
# $ Rat_ID                   : chr  "KH L07 A1" "KH L07 B1" "KH L07 C1" "KH L07 C2" ...
# $ Generation               : chr  "F56" "F56" "F56" "F56" ...
# $ Line                     : chr  "bLR" "bLR" "bLR" "bLR" ...
# $ Treatment_group          : chr  "bLR EE" "bLR EE" "bLR EE + SD" "bLR EE + SD" ...
# $ Enrichment               : chr  "EE" "EE" "EE" "EE" ...
# $ SocialDefeat             : chr  "NIL" "NIL" "SD" "SD" ...
# $ brain_region             : chr  "hippocampus" "hippocampus" "hippocampus" "hippocampus" ...
# $ date_of_sac              : chr  "2-Nov-17" "2-Nov-17" "2-Nov-17" "2-Nov-17" ...
# $ date_of_dissection       : chr  "16/1/2019" "13/3/2019" "16/1/2019" "16/1/2019" ...
# $ hemisphere_dissected     : chr  "right" "right" "left" "left" ...
# $ date_of_RNA_extraction   : chr  "28.2.2019" "13.3.2019" "28.2.2019" "28.2.2019" ...
# $ RNA_conc                 : num  80.1 56.7 103.9 97.5 102.9 ...
# $ ratio_260_280            : num  2.06 2.01 2.02 2.06 2.07 2.01 2.01 2.05 2.13 2.16 ...
# $ SampleID                 : chr  "Sample_126059" "Sample_127120" "Sample_126066" "Sample_126067" ...

table(Sample_MetaData_HC_NoOutliers$brain_region)
# hippocampus  
# 25 
#Good.

#But it may still be in the wrong order:
Sample_MetaData_HC_NoOutliers$Sequencing_core_sample_ID
#Yep.

Sample_MetaData_HC_NoOutliers_Ordered<-Sample_MetaData_HC_NoOutliers[order(Sample_MetaData_HC_NoOutliers$Sequencing_core_sample_ID),]
Sample_MetaData_HC_NoOutliers_Ordered$Sequencing_core_sample_ID
# [1] 126043 126044 126045 126046 126048 126049 126050 126051 126052 126053 126054 126055 126056 126057 126058 126059 126065 126066
# [19] 126067 126068 126069 126070 127119 127120 127121

#Looks good now.

cbind(colnames(HC_RNASeq_Log2_Filtered), Sample_MetaData_HC_NoOutliers_Ordered$SampleID)
#Looks good now.
sum(colnames(HC_RNASeq_Log2_Filtered)==Sample_MetaData_HC_NoOutliers_Ordered$SampleID)
#[1] 25
#Yep, they're all in the same order.

cbind(names(HC_RNASeq_CountsPerSample), Sample_MetaData_HC_NoOutliers_Ordered$SampleID)
#Same order

Sample_MetaData_HC_NoOutliers_Ordered$CountsPerSample<-HC_RNASeq_CountsPerSample

Sample_BehavHormonalData$SampleID<-paste("Sample_", Sample_BehavHormonalData$HC_Sequencing_core_sample_ID, sep="")
Sample_BehavHormonalData$SampleID
# [1] "Sample_126059" "Sample_126060" "Sample_127120" "Sample_126062" "Sample_126066" "Sample_126067" "Sample_126068"
# [8] "Sample_126069" "Sample_126070" "Sample_126048" "Sample_126049" "Sample_126050" "Sample_126051" "Sample_126052"
# [15] "Sample_126053" "Sample_126054" "Sample_126055" "Sample_126043" "Sample_126044" "Sample_126045" "Sample_126046"
# [22] "Sample_127119" "Sample_126040" "Sample_126041" "Sample_126042" "Sample_126063" "Sample_127121" "Sample_126065"
# [29] "Sample_126056" "Sample_126057" "Sample_126058"

Sample_BehavHormonalData_HC_NoOutliers<-Sample_BehavHormonalData[which(Sample_BehavHormonalData$SampleID%in%colnames(HC_RNASeq_Log2_Filtered)), ]
str(Sample_BehavHormonalData_HC_NoOutliers)
#'data.frame':	25 obs. of  59 variables:

Sample_BehavHormonalData_HC_NoOutliers_Ordered<-Sample_BehavHormonalData_HC_NoOutliers[order(Sample_BehavHormonalData_HC_NoOutliers$HC_Sequencing_core_sample_ID),]
Sample_BehavHormonalData_HC_NoOutliers_Ordered$HC_Sequencing_core_sample_ID
# [1] 126043 126044 126045 126046 126048 126049 126050 126051 126052 126053 126054 126055 126056 126057 126058 126059 126065 126066
# [19] 126067 126068 126069 126070 127119 127120 127121

cbind(colnames(HC_RNASeq_Log2_Filtered), Sample_BehavHormonalData_HC_NoOutliers_Ordered$SampleID)
#Looks good now.
sum(colnames(HC_RNASeq_Log2_Filtered)==Sample_BehavHormonalData_HC_NoOutliers_Ordered$SampleID)
#[1] 25
#yep

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples")
write.csv(Sample_BehavHormonalData_HC_NoOutliers_Ordered, "Sample_BehavHormonalData_HC_NoOutliers_Ordered.csv")


#Getting the factors set-up properly:
str(Sample_MetaData_HC_NoOutliers_Ordered)
Sample_MetaData_HC_NoOutliers_Ordered$Generation<-as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Generation)
levels(Sample_MetaData_HC_NoOutliers_Ordered$Generation)
#[1] "F53" "F56"
table(Sample_MetaData_HC_NoOutliers_Ordered$Generation)
# F53 F56 
# 5  20 

#Let's make F56 the reference level
Sample_MetaData_HC_NoOutliers_Ordered$Generation<-relevel(Sample_MetaData_HC_NoOutliers_Ordered$Generation, ref="F56")
levels(Sample_MetaData_HC_NoOutliers_Ordered$Generation)
#[1] "F56" "F53"

Sample_MetaData_HC_NoOutliers_Ordered$Enrichment<-as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Enrichment)
levels(Sample_MetaData_HC_NoOutliers_Ordered$Enrichment)
#[1] "EE"  "NIL"
Sample_MetaData_HC_NoOutliers_Ordered$Enrichment<-relevel(Sample_MetaData_HC_NoOutliers_Ordered$Enrichment, ref="NIL")
levels(Sample_MetaData_HC_NoOutliers_Ordered$Enrichment)
#[1] "NIL" "EE" 

Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat<-as.factor(Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat)
levels(Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat)
#[1] "NIL" "SD" 


######################

#Double-checking whether any of the batch-related variables are unevenly distributed by group:

table(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_HC_NoOutliers_Ordered$Generation)
#               F56 F53
# bLR EE         2   3
# bLR EE + SD    5   2
# bLR NIL        5   0
# bLR NIL + SD   8   0

fisher.test(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_HC_NoOutliers_Ordered$Generation)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group and Sample_MetaData_HC_NoOutliers_Ordered$Generation
# p-value = 0.02351
# alternative hypothesis: two.sided

fisher.test(Sample_MetaData_HC_NoOutliers_Ordered$Enrichment, Sample_MetaData_HC_NoOutliers_Ordered$Generation)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_HC_NoOutliers_Ordered$Enrichment and Sample_MetaData_HC_NoOutliers_Ordered$Generation
# p-value = 0.01491
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 1.23054     Inf
# sample estimates:
# odds ratio 
# Inf 

#Generation is unevenly distributed. :(
#But maybe it won't cause problems? F53 and F56 are pretty close, and in the NAcc dataset it was really the F49 collections that looked weird.

table(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)

#                   13/3/2019 16/1/2019 17/1/2019
# bLR EE               1         4         0
# bLR EE + SD          1         4         2
# bLR NIL              1         4         0
# bLR NIL + SD         0         4         4

#Hmmm..dissection date may be unevenly distributed across groups but not in a way that can be easily included as a co-variate. Meh.

fisher.test(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group and Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection
# p-value = 0.3121
# alternative hypothesis: two.sided


table(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)
#               1.3.2019 13.3.2019 28.2.2019
# LR EE              2         1         2
# LR EE + SD         0         1         6
# LR Nil             2         1         2
# LR Nil + SD        2         0         6
#.... and RNA-extraction date too.


#Both of those variables would be hard to include as co-variates - the sample size simply isn't big enough. I guess we can check and see if they correlate with the PCs, although some of the subgroups are going to make it hard to tell. :(

#Although I think they may be partially redundant. Let's check:

table(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection, Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)
#             1.3.2019 13.3.2019 28.2.2019
# 13/3/2019        0         3         0
# 16/1/2019        6         0        10
# 17/1/2019        0         0         6
#Yep. There are basically just 4 batches.  Let's create a combined variable:

Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch<-paste(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection, Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)
table(Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch)

# 13/3/2019 13.3.2019  16/1/2019 1.3.2019 16/1/2019 28.2.2019 17/1/2019 28.2.2019 
# 3                   6                  10                   6 

table(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch)

#                         13/3/2019 13.3.2019 16/1/2019 1.3.2019 16/1/2019 28.2.2019 17/1/2019 28.2.2019
# bLR EE                         1                  2                   2                   0
# bLR EE + SD                    1                  0                   4                   2
# bLR NIL                        1                  2                   2                   0
# bLR NIL + SD                   0                  2                   2                   4

fisher.test(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group and Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch
# p-value = 0.2899
# alternative hypothesis: two.sided

#Not significantly related, but some subgroups are super small.

#Out of curiousity, I wonder what things looked like in the original data before Fan removed the outliers:

table(Sample_MetaData$date_of_dissection[Sample_MetaData$brain_region=="hippocampus"], Sample_MetaData$date_of_RNA_extraction[Sample_MetaData$brain_region=="hippocampus"])
#               1.3.2019 13.3.2019 28.2.2019
# 13/3/2019        0         3         0
# 16/1/2019        7         0        14
# 17/1/2019        1         0         6

#So the outliers weren't from a single dissection/extraction group - 4 were from 28/2 and 2 were from 1/3 extraction groups, most were from 16/1 dissection group (but a minority of samples -5- collected that day)
#So the outliers weren't simply a bad batch. Interesting.

table(Sample_MetaData$date_of_dissection[Sample_MetaData$brain_region=="hippocampus"], Sample_MetaData$Treatment_group[Sample_MetaData$brain_region=="hippocampus"])

#             bLR EE bLR EE + SD bLR NIL bLR NIL + SD
# 13/3/2019      1           1       1            0
# 16/1/2019      5           5       7            4
# 17/1/2019      1           2       0            4

table(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection, Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)
#                bLR EE bLR EE + SD bLR NIL bLR NIL + SD
# 13/3/2019      1           1       1            0
# 16/1/2019      4           4       4            4
# 17/1/2019      0           2       0            4

#The outlier samples also weren't from a single group - 1 from EE, 1 from EE+SD, and 3 from NIL


#Are these also related to generation?
table(Sample_MetaData_HC_NoOutliers_Ordered$Generation, Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch)
#           13/3/2019 13.3.2019 16/1/2019 1.3.2019 16/1/2019 28.2.2019 17/1/2019 28.2.2019
# F56                   2                  4                   8                   6
# F53                   1                  2                   2                   0
#Partially, but not overwhelmingly (unlike for the Nacc)


anova(lm(RNA_conc~Treatment_group, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: RNA_conc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Treatment_group  3  2235.4  745.13  1.1717 0.3442
# Residuals       21 13354.8  635.94  

summary.lm(lm(Sample_MetaData_HC_NoOutliers_Ordered$RNA_conc~Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat+Sample_MetaData_HC_NoOutliers_Ordered$Enrichment))
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -37.719 -21.130   1.121   9.430  59.381 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                            69.670      9.222   7.555 1.51e-07 ***
#   Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeatSD   16.749     10.064   1.664    0.110    
# Sample_MetaData_HC_NoOutliers_Ordered$EnrichmentEE      9.960      9.868   1.009    0.324    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 24.64 on 22 degrees of freedom
# Multiple R-squared:  0.1434,	Adjusted R-squared:  0.06551 
# F-statistic: 1.841 on 2 and 22 DF,  p-value: 0.1822

pdf("HC_Boxplot_RNAConc_vs_TreatmentGroup.pdf", width=5, height=5)
boxplot(RNA_conc~Treatment_group, data=Sample_MetaData_HC_NoOutliers_Ordered, cex.axis=0.5, ylab="RNA Concentration")
dev.off()


anova(lm(ratio_260_280~Treatment_group, data=Sample_MetaData_HC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: ratio_260_280
# Df   Sum Sq   Mean Sq F value Pr(>F)
# Treatment_group  3 0.014214 0.0047379  1.4673 0.2521
# Residuals       21 0.067810 0.0032291 

summary.lm(lm(ratio_260_280~Enrichment+SocialDefeat, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.114556 -0.034652 -0.004652  0.015348  0.155444 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     2.10456    0.02132  98.730   <2e-16 ***
#   EnrichmentEE   -0.03711    0.02281  -1.627    0.118    
# SocialDefeatSD -0.01990    0.02326  -0.856    0.401    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05695 on 22 degrees of freedom
# Multiple R-squared:  0.1302,	Adjusted R-squared:  0.05108 
# F-statistic: 1.646 on 2 and 22 DF,  p-value: 0.2157


anova(lm(CountsPerSample~Treatment_group, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value  Pr(>F)  
# Treatment_group  3 1.2579e+15 4.1929e+14  2.4515 0.09169 .
# Residuals       21 3.5917e+15 1.7103e+14                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pdf("HC_Boxplot_CountsPerSample_vs_TreatmentGroup.pdf", width=5, height=5)
boxplot(CountsPerSample~Treatment_group, data=Sample_MetaData_HC_NoOutliers_Ordered, cex.axis=0.5, ylab="Total Counts Per Sample")
dev.off()

summary.lm(lm(CountsPerSample~Enrichment+SocialDefeat, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -27910540  -9797257   1563066   8344724  30963128 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)   
# (Intercept)    15588202    5099205   3.057  0.00578 **
# EnrichmentEE   10659538    5456478   1.954  0.06359 . 
# SocialDefeatSD  3474179    5564538   0.624  0.53882   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 13620000 on 22 degrees of freedom
# Multiple R-squared:  0.1581,	Adjusted R-squared:  0.08155 
# F-statistic: 2.065 on 2 and 22 DF,  p-value: 0.1506

summary.lm(lm(CountsPerSample~Enrichment*SocialDefeat, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -31784960  -9158711   1013234   7666067  25538940 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                 21012390    5848608   3.593  0.00171 **
# EnrichmentEE                 -188840    8271181  -0.023  0.98200   
# SocialDefeatSD              -5340128    7455542  -0.716  0.48172   
# EnrichmentEE:SocialDefeatSD 18112916   10687579   1.695  0.10490   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 13080000 on 21 degrees of freedom
# Multiple R-squared:  0.2594,	Adjusted R-squared:  0.1536 
# F-statistic: 2.452 on 3 and 21 DF,  p-value: 0.09169

pdf("CountsPerSample_vsSeqCoreID_byEE.pdf", height=5, width=8)
plot(CountsPerSample~Sequencing_core_sample_ID, data=Sample_MetaData_HC_NoOutliers_Ordered, col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Enrichment))
dev.off()


#Does RNAconc correlate with purity?

pdf("RNAconc_vs_Purity.pdf", width=5, height=5)
plot(Sample_MetaData_HC_NoOutliers_Ordered$RNA_conc, Sample_BehavHormonalData_HC_NoOutliers_Ordered$HC_ratio_260_280)
dev.off()
#maybe? purity isn't varying much though - between 2.0 and 2.25
#noteworthy: concentrations in general are pretty low (between 40-140)
cor(Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280, Sample_MetaData_HC_NoOutliers_Ordered$RNA_conc) 
#[1] -0.3536363

pdf("CountsPerSample_vs_Purity.pdf", width=5, height=5)
plot(HC_RNASeq_CountsPerSample, Sample_BehavHormonalData_HC_NoOutliers_Ordered$HC_ratio_260_280)
dev.off()
cor(HC_RNASeq_CountsPerSample, Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280)
#[1] -0.4273725
#huh. More pure has fewer counts??? 

pdf("CountsPerSample_vs_Conc.pdf", width=5, height=5)
plot(HC_RNASeq_CountsPerSample~Sample_MetaData_HC_NoOutliers_Ordered$RNA_conc)
dev.off()

cor(HC_RNASeq_CountsPerSample, Sample_MetaData_HC_NoOutliers_Ordered$RNA_conc)
#[1] 0.166099


#What about generation?

anova(lm(RNA_conc~Generation, data=Sample_MetaData_HC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: RNA_conc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Generation  1   939.4  939.42  1.4748 0.2369
# Residuals  23 14650.8  636.99 

anova(lm(ratio_260_280~Generation, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: ratio_260_280
# Df   Sum Sq   Mean Sq F value Pr(>F)
# Generation  1 0.000784 0.0007840   0.222  0.642
# Residuals  23 0.081240 0.0035322 

anova(lm(CountsPerSample~Generation, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value Pr(>F)
# Generation  1 9.6438e+13 9.6438e+13  0.4667 0.5014
# Residuals  23 4.7531e+15 2.0666e+14   

#No relationship with generation.

anova(lm(RNA_conc~DissectionExtractionBatch, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: RNA_conc
# Df  Sum Sq Mean Sq F value Pr(>F)
# DissectionExtractionBatch  3   825.5  275.15  0.3914 0.7605
# Residuals                 21 14764.7  703.08 

pdf("HC_Boxplot_RNA_conc_vs_DissectionExtractionBatch.pdf", width=5, height=5)
boxplot(Sample_MetaData_HC_NoOutliers_Ordered$RNA_conc~Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch, cex.axis=0.5)
dev.off()

anova(lm(ratio_260_280~DissectionExtractionBatch, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: ratio_260_280
# Df   Sum Sq   Mean Sq F value  Pr(>F)   
# DissectionExtractionBatch  3 0.040757 0.0135858  6.9136 0.00205 **
#   Residuals                 21 0.041267 0.0019651                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("HC_Boxplot_ratio_260_280_vs_DissectionExtractionBatch.pdf", width=5, height=5)
boxplot(Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280~Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch, cex.axis=0.5)
dev.off()

#Purity definitely differs by dissection/extraction batch
#Although it is still within a very narrow range, so I don't know if it would functionally matter (or might just represent batches of nanodrop readouts)

anova(lm(ratio_260_280~date_of_RNA_extraction, data=Sample_MetaData_HC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: ratio_260_280
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# date_of_RNA_extraction  2 0.038997 0.0194985  9.9697 0.0008276 ***
#   Residuals              22 0.043027 0.0019558                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#And batches of extraction seem to drive it.

pdf("HC_Boxplot_ratio_260_280_vs_ExtractionBatch.pdf", width=5, height=5)
boxplot(Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280~Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction, cex.axis=0.5, ylab="RNA Quality: Ratio 260/280")
dev.off()

summary.lm(lm(ratio_260_280~date_of_RNA_extraction, data=Sample_MetaData_HC_NoOutliers_Ordered))
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.078333 -0.020000 -0.003125  0.020000  0.121667 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      2.13833    0.01805 118.438  < 2e-16 ***
#   date_of_RNA_extraction13.3.2019 -0.12833    0.03127  -4.104 0.000468 ***
#   date_of_RNA_extraction28.2.2019 -0.07521    0.02117  -3.552 0.001784 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04422 on 22 degrees of freedom
# Multiple R-squared:  0.4754,	Adjusted R-squared:  0.4


anova(lm(CountsPerSample~DissectionExtractionBatch, data=Sample_MetaData_HC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value Pr(>F)
# DissectionExtractionBatch  3 1.1683e+15 3.8942e+14  2.2215 0.1156
# Residuals                 21 3.6813e+15 1.7530e+14 


pdf("HC_Boxplot_CountsPerSample_vs_DissectionExtractionBatch.pdf", width=8, height=5)
boxplot(CountsPerSample~DissectionExtractionBatch, data=Sample_MetaData_HC_NoOutliers_Ordered, cex.axis=0.5)
dev.off()

pdf("HC_Boxplot_CountsPerSample_vs_DissectionBatch.pdf", width=5, height=5)
boxplot(CountsPerSample~date_of_dissection, data=Sample_MetaData_HC_NoOutliers_Ordered, ylab="Library Size: Counts per Sample", cex.axis=0.5)
dev.off()

anova(lm(CountsPerSample~date_of_dissection, data=Sample_MetaData_HC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value Pr(>F)  
# date_of_dissection  2 1.1630e+15 5.8148e+14    3.47  0.049 *
#   Residuals          22 3.6866e+15 1.6757e+14                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


anova(lm(CountsPerSample~date_of_dissection+Treatment_group, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value  Pr(>F)  
# date_of_dissection  2 1.1630e+15 5.8148e+14  4.2033 0.03080 *
#   Treatment_group     3 1.0581e+15 3.5271e+14  2.5496 0.08623 .
# Residuals          19 2.6284e+15 1.3834e+14                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary.lm(lm(CountsPerSample~date_of_dissection+Enrichment+SocialDefeat, data=Sample_MetaData_HC_NoOutliers_Ordered))

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -24133758  -8157276   -541451  10389342  19406936 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  32566862    8036534   4.052 0.000622 ***
#   date_of_dissection16/1/2019 -20485044    7853727  -2.608 0.016821 *  
#   date_of_dissection17/1/2019 -16032265    9688178  -1.655 0.113568    
# EnrichmentEE                  9478287    5035024   1.882 0.074403 .  
# SocialDefeatSD                4385032    5709180   0.768 0.451422    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 12310000 on 20 degrees of freedom
# Multiple R-squared:  0.3752,	Adjusted R-squared:  0.2502 
# F-statistic: 3.002 on 4 and 20 DF,  p-value: 0.0431

table(Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat, Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)
#       13/3/2019 16/1/2019 17/1/2019
# NIL         2         8         0
# SD          1         8         6

fisher.test(Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat, Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat and Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection
# p-value = 0.0645
# alternative hypothesis: two.sided

table(Sample_MetaData_HC_NoOutliers_Ordered$Enrichment, Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)

#       13/3/2019 16/1/2019 17/1/2019
# NIL         1         8         4
# EE          2         8         2

fisher.test(Sample_MetaData_HC_NoOutliers_Ordered$Enrichment, Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_HC_NoOutliers_Ordered$Enrichment and Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection
# p-value = 0.7195
# alternative hypothesis: two.sided

table(Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat, Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)
#         1.3.2019 13.3.2019 28.2.2019
# NIL        4         2         4
# SD         2         1        12

fisher.test(Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat, Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat and Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction
# p-value = 0.1453
# alternative hypothesis: two.sided



#In review: Relationships that might be important:

#Generation relates to enrichment (all F53 are enriched - but hopefully it won't matter because it is a similar generation to F56)
#Date of dissection relates to social defeat (trend), but hopefully in a way that won't matter (17/1 is all SD - but very similar date to 16/1)
# Counts per sample is slightly higher in EE (trend) - opposite of NAcc. Also relates to date of dissection.
# Purity doesn't vary by treatment group, but varies by extraction batch, counts per sample, and concentration. 
# RNA conc doesn't vary by treatment group (unlike NAcc)
# Sequencing ID suggests that three samples were re-run - other batch variable?


#Meh. This sample size probably can't handle including any of these as co-variates. Hopefully they won't matter...



###############################################
#PCA/MDS

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples")

#MDS plots

png("MDSplots_Group_DissectionBatch.png")
par(mfrow=c(2,2))
col.group <- Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group
#levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.batch <- as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)
levels(col.batch) <-  brewer.pal(nlevels(col.batch), "Set2")
col.batch <- as.character(col.batch)
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)))
title(main="A. Sample groups")
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)), dim=c(3,4))
title(main="B. Sample groups")
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)), col=col.batch, dim=c(1,2))
title(main="C. Dissection Batch")
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)), col=col.batch, dim=c(3,4))
title(main="D. Dissection Batch")
dev.off()

png("MDSplots_Group_ExtractionBatch.png")
par(mfrow=c(2,2))
col.group <- Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group
#levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.batch <- as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)
levels(col.batch) <-  brewer.pal(nlevels(col.batch), "Set2")
col.batch <- as.character(col.batch)
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)))
title(main="A. Sample groups")
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)), dim=c(3,4))
title(main="B. Sample groups")
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)), col=col.batch, dim=c(1,2))
title(main="C. Extraction Batch")
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)), col=col.batch, dim=c(3,4))
title(main="D. Extraction Batch")
dev.off()
#Extraction and/or Dissection Batch are related to many of the top sources of variation in the data.
#But so does group. So that's cool. 


png("MDSplots_Group_Generation.png")
par(mfrow=c(2,2))
col.group <- Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group
#levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
#col.group <- as.character(col.group)
col.batch <- Sample_MetaData_HC_NoOutliers_Ordered$Generation
#levels(col.batch) <-  brewer.pal(nlevels(col.batch), "Set2")
#col.batch <- as.character(col.batch)
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)))
title(main="A. Sample groups")
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)), dim=c(3,4))
title(main="B. Sample groups")
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Generation)), col=as.numeric(as.factor(col.batch)), dim=c(1,2))
title(main="C. Generation")
plotMDS(HC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Generation)), col=as.numeric(as.factor(col.batch)), dim=c(3,4))
title(main="D. Generation")
dev.off()
#Generation seems to be irrelevant. *Phew*


#############################

#Switching over to principal components analysis (PCA)  to see if we can generate some concise co-variates that can be used to control for technical variation.


#I'm used to working with microarray data. RNA-Seq data needs to be centered and scaled before doing this.
#but in which direction?  I'm inclined to center and scale it by gene, so that highly variable genes aren't driving the effect.
#On the other hand, normally the normalization is performed by sample in microarray. It seems like the overall distribution of the genes looks very similar across subjects, so perhaps that is less important.

HC_RNASeq_Log2_Filtered_Scaled<-t(scale(t(HC_RNASeq_Log2_Filtered), center=T, scale=T))


pcaNormFiltered<-prcomp(t(HC_RNASeq_Log2_Filtered_Scaled))
tmp<-pcaNormFiltered$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")

PC1<-pcaNormFiltered$x[,1]
PC2<-pcaNormFiltered$x[,2]

PC3<-pcaNormFiltered$x[,3]
PC4<-pcaNormFiltered$x[,4]

tmp<-data.frame(gene_id=row.names(pcaNormFiltered$rotation),pcaNormFiltered$rotation[,c(1:4)])

HC_Eigenvectors_Annotated<-join(tmp, HC_RNASeq_Log2_Annotated[, c(1, 27:29)], by="gene_id", type="left")

write.csv(HC_Eigenvectors_Annotated, "Eigenvectors_PCA_1_4.csv")

#Output a scree plot for the PCA:
png("PCA Scree Plot1.png")
plot(summary(pcaNormFiltered)$importance[2,]~(c(1:25)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 explains around 30% of the variation in the dataset, all the others are less than 15%. Pretty normal.

png("PCA Scree Plot2.png")
plot(summary(pcaNormFiltered)$importance[3,]~(c(1:25)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

#Quickly coloring plots using group as factor:
levels(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("HC_PC1vsPC2_vsGroup.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group))
legend(min(PC2), max(PC1)-5, c("LR EE", "LR EE + SD", "LR Nil", "LR Nil + SD"), text.col=c(1, 2, 3, 4), pch=19, col=c(1, 2, 3, 4))
dev.off()
#PC1 isn't driven by outliers or group status
#Looks like Fan's graph - good. 


levels(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Generation))

png("HC_PC1vsPC2_vsGeneration.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$Generation))
legend(min(PC2), max(PC1)-5, c("F53", "F56"), text.col=c(1, 2), pch=19, col=c(1, 2))
dev.off()
#Not generation

levels(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_sac))
#[1] "2-Nov-17"  "24-Oct-17" "28-Oct-17" "9-Feb-17" 

png("HC_PC1vsPC2_vsDateOfSacrifice.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_sac))
legend(min(PC2), max(PC1)-5, c("2-Nov-17", "24-Oct-17" ,"28-Oct-17", "9-Feb-17" ), text.col=c(1, 2, 3, 4), pch=19, col=c(1, 2, 3, 4))
dev.off()
#Not date of sacrifice.

levels(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction))
#[1] "1.3.2019"  "13.3.2019" "28.2.2019"

png("HC_PC1vsPC2_vsDateOfExtraction.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction))
legend(min(PC2), max(PC1)-5, c("1.3.2019",  "13.3.2019", "28.2.2019"), text.col=c(1, 2, 3), pch=19, col=c(1, 2, 3))
dev.off()
#Date of Extraction is part of it, but not the full story.

png("HC_PC1_vsDateOfExtraction.png")
boxplot(PC1~Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)
dev.off()
#The 13-3 batch is particularly related to PC1

png("HC_PC2_vsDateOfExtraction.png")
boxplot(PC2~Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)
dev.off()

png("HC_PC3_vsDateOfExtraction.png")
boxplot(PC3~Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)
dev.off()

png("HC_PC4_vsDateOfExtraction.png")
boxplot(PC4~Sample_MetaData_HC_NoOutliers_Ordered$date_of_RNA_extraction)
dev.off()
#13-3 drives PC4


png("HC_PC1vsRNAconc.png")
plot(PC1~Sample_MetaData_HC_NoOutliers_Ordered$RNA_conc)
dev.off()

png("HC_PC2vsRNAconc.png")
plot(PC2~Sample_MetaData_HC_NoOutliers_Ordered$RNA_conc)
dev.off()
#Not RNA Concentration.

png("HC_PC3vsRNAconc.png")
plot(PC3~Sample_MetaData_HC_NoOutliers_Ordered$RNA_conc)
dev.off()
#maybe?

png("HC_PC1vsPurity.png")
plot(PC1~Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280)
dev.off()

png("HC_PC2vsPurity.png")
plot(PC2~Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280)
dev.off()
#Maybe?
cor(PC2,Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280)
#[1] -0.2891265

png("HC_PC3vsPurity.png")
plot(PC3~Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280)
dev.off()
#Maybe? outlier

png("HC_PC4vsPurity.png")
plot(PC4~Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280)
dev.off()
#Maybe?
cor(PC4,Sample_MetaData_HC_NoOutliers_Ordered$ratio_260_280)
#[1] 0.1334794

png("HC_PC1vsCountsPerSample.png")
plot(PC1~Sample_MetaData_HC_NoOutliers_Ordered$CountsPerSample)
dev.off()

png("HC_PC2vsCountsPerSample.png")
plot(PC2~Sample_MetaData_HC_NoOutliers_Ordered$CountsPerSample)
dev.off()
#Looks potentially convincing. Not sure what to do with that.
cor(PC2,Sample_MetaData_HC_NoOutliers_Ordered$CountsPerSample)
#[1] 0.6902021
#Ouch.

png("HC_PC3vsCountsPerSample.png")
plot(PC3~Sample_MetaData_HC_NoOutliers_Ordered$CountsPerSample)
dev.off()
#Also maybe convincing, but has big outlier.

png("HC_PC4vsCountsPerSample.png")
plot(PC4~Sample_MetaData_HC_NoOutliers_Ordered$CountsPerSample)
dev.off()


levels(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection))
#[1] "13/3/2019" "16/1/2019" "17/1/2019"

table(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)
# 13/3/2019 16/1/2019 17/1/2019 
# 3        16         6
table(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection, Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)
#              bLR EE bLR EE + SD bLR NIL bLR NIL + SD
# 13/3/2019      1           1       1            0
# 16/1/2019      4           4       4            4
# 17/1/2019      0           2       0            4

Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection<-as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)
Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection<-relevel(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection, ref="16/1/2019")
levels(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)
#[1] "16/1/2019" "13/3/2019" "17/1/2019"

png("HC_PC1vsPC2_vsDateOfDissection.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection))
legend(min(PC2), max(PC1)-5, c("13/3/2019",  "16/1/2019", "17/1/2019"), text.col=c(1, 2, 3,4), pch=19, col=c(1, 2, 3,4))
dev.off()
#That's part of the story too, but not the full story. In particular, the three 13/3/2019 samples look different from the rest, but 17/1 looks different too
#Dissection is related to both pc1 and pc2

png("HC_PC1_vsDateOfDissection.png")
boxplot(PC1~Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)
dev.off()
anova(lm(PC1~Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection))

# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection  2  42638 21318.8  6.3005 0.006866 **
#   Residuals                                                22  74441  3383.7                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Dissection batch definitely matters... and is unevenly distributed across groups.

anova(lm(PC1~Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection+Sample_MetaData_HC_NoOutliers_Ordered$Enrichment+Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat))
# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection  2  42638 21318.8  6.2968 0.007568 **
# Sample_MetaData_HC_NoOutliers_Ordered$Enrichment          1    537   537.3  0.1587 0.694583   
# Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat        1   6191  6190.6  1.8285 0.191398   
# Residuals                                                20  67713  3385.6                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Dissection batch still matters after controlling for treatment group. :(


png("HC_PC2_vsDateOfDissection.png")
boxplot(PC2~Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)
dev.off()

png("HC_PC3_vsDateOfDissection.png")
boxplot(PC3~Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)
dev.off()
#PC3 is an outlier.

png("HC_PC4_vsDateOfDissection.png")
boxplot(PC4~Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection)
dev.off()
#The 13/3 group again drives PC4. Sigh.

anova(lm(PC4~Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection+Sample_MetaData_HC_NoOutliers_Ordered$Enrichment+Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat))

# Analysis of Variance Table
# 
# Response: PC4
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection  2 19334.4  9667.2 21.6506 9.913e-06 ***
# Sample_MetaData_HC_NoOutliers_Ordered$Enrichment          1   173.7   173.7  0.3891    0.5398    
# Sample_MetaData_HC_NoOutliers_Ordered$SocialDefeat        1  2188.9  2188.9  4.9023    0.0386 *  
#   Residuals                                                20  8930.2   446.5                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#... and that is true after controlling for treatment group. :(

levels(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch) )
#[1]"13/3/2019 13.3.2019" "16/1/2019 1.3.2019"  "16/1/2019 28.2.2019" "17/1/2019 28.2.2019"

png("HC_PC1vsPC2_vsDateOfExtractionDissection.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch))
legend(min(PC2), max(PC1)-5, c("13/3/2019 13.3.2019" , "16/1/2019 1.3.2019",   "16/1/2019 28.2.2019",  "17/1/2019 28.2.2019" ), text.col=c(1, 2, 3,4, 5), pch=19, col=c(1, 2, 3,4, 5))
dev.off()

#Seems more like dissection batch than extraction batch (not much difference in the 16/1 samples)

png("HC_PC1_vsExtractionDissection.png")
boxplot(PC1~Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch)
dev.off()

anova(lm(PC1~Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch))

# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch  3  42909 14302.9  4.0496 0.02033 *
#   Residuals                                                       21  74169  3531.9                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Not a crazy strong relationship, but this is a super unbalanced design - we probably shouldn't even be using ANOVA...
#From the graph, looks like it is mostly driven by dissection, not extraction.

png("HC_PC2_vsExtractionDissection.png")
boxplot(PC2~Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch)
dev.off()

png("HC_PC3_vsExtractionDissection.png")
boxplot(PC3~Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch)
dev.off()
#PC3 is an outlier.

png("HC_PC4_vsExtractionDissection.png")
boxplot(PC4~Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch)
dev.off()
#The 13/3 group again drives PC4. Sigh.

anova(lm(PC4~Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch))

# Analysis of Variance Table
# 
# Response: PC4
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Sample_MetaData_HC_NoOutliers_Ordered$DissectionExtractionBatch  3  19336  6445.2  11.987 8.683e-05 ***
#   Residuals                                                       21  11292   537.7                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Definitely matters for PC4.


levels(as.factor(Sample_MetaData_HC_NoOutliers_Ordered$hemisphere_dissected))
#[1] "left"  "right"

png("HC_PC1vsPC2_vsHemisphere.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$hemisphere_dissected))
legend(min(PC2), max(PC1)-5, c("left", "right"), text.col=c(1, 2), pch=19, col=c(1, 2))
dev.off()
#Not hemisphere.

png("HC_PC1_vsGroup.png")
boxplot(PC1~Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)
dev.off()
#That's cool - social defeat looks like it is related to PC1.

png("HC_PC2_vsGroup.png")
boxplot(PC2~Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)
dev.off()

png("HC_PC3_vsGroup.png")
boxplot(PC3~Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)
dev.off()
#and EE looks maybe related to PC3

png("HC_PC4_vsGroup.png")
boxplot(PC4~Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)
dev.off()


colnames(Sample_MetaData_HC_NoOutliers_Ordered)


#I should probably take a look at the eigenvectors and see if there are any hints. :(
#Nope, nothing obvious.

write.csv(cor(cbind(PC1, PC2, PC3, PC4, as.matrix(Sample_MetaData_HC_NoOutliers_Ordered[, c(13:14, 16)]), as.matrix(Sample_BehavHormonalData_HC_NoOutliers_Ordered[,c(10:21,35:46, 48:51)])), use="pairwise.complete.obs"), "PCA_vsNumericVar_CorMatrix.csv")
#Note:  many of these variables are really really not distributed normally - let's try running a spearman instead.

write.csv(cor(cbind(PC1, PC2, PC3, PC4, as.matrix(Sample_MetaData_HC_NoOutliers_Ordered[, c(13:14, 16)]), as.matrix(Sample_BehavHormonalData_HC_NoOutliers_Ordered[,c(10:21,35:46, 48:51)])), use="pairwise.complete.obs", method="spearman"), "PCA_vsNumericVar_CorMatrix_Spearman.csv")

#I'm not sure how many of these relationships are actually likely to be derived from a reasonable sample size.

sum(is.na(Sample_BehavHormonalData_HC_NoOutliers_Ordered$Oxytocin..pg.ml.))
#[1] 0
#Oh- there's apparently actually quite a bit of overlap with the subjects with hormone measurements. That's exciting.
pdf("PC1_vs_Oxytocin_byDissection.pdf", width=5, height=5)
plot(PC1~Sample_BehavHormonalData_HC_NoOutliers_Ordered$Oxytocin..pg.ml., col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection))
dev.off()
#could be real - esp. after controlling for dissection. Weird.

sum(is.na(Sample_BehavHormonalData_HC_NoOutliers_Ordered$time_approaching_stimulus_animal_Video))
#[1] 0
#Awesome!

pdf("PC1_vs_ApproachVideo_byDissection.pdf", width=5, height=5)
plot(PC1~Sample_BehavHormonalData_HC_NoOutliers_Ordered$time_approaching_stimulus_animal_Video, col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection))
dev.off()
#Ouch - that variable is completely bimodal.

pdf("PC1_vs_EPM_byDissection.pdf", width=5, height=5)
plot(PC1~Sample_BehavHormonalData_HC_NoOutliers_Ordered$time_open_arms_EPM, col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection))
dev.off()
#Meh- seems like a batch effect.

pdf("PC2_vs_CenterOfOpenField_byDissection.pdf", width=5, height=5)
plot(PC2~Sample_BehavHormonalData_HC_NoOutliers_Ordered$time_centre_open_field, col=as.factor(Sample_MetaData_HC_NoOutliers_Ordered$date_of_dissection))
dev.off()
#Pretty unconvincing.

#So I think we are left with a resounding "Maybe" regarding the question of whether we need to control for dissection batch. 13.3.2019 sems to matter a lot, but it is only 3 subjects, spread out over 3 groups. But 17/1 seems to also matter... and dissection batch correlates with sd (trend), so it is potentially an important confound
#Counts per sample also strongly related to PC2, and may be elevated in EE (trend)

##############################

#Limma-voom analysis:

#

#Trying the design with the most co-variates first, although I suspect that it is probably way more df than this sample can handle:
#I did not include the interaction term

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/JustEESD_CovDissectionCountsPerSample")

design <- model.matrix(~Enrichment+SocialDefeat+date_of_dissection+CountsPerSample, data=Sample_MetaData_HC_NoOutliers_Ordered)

design
# (Intercept) EnrichmentEE SocialDefeatSD date_of_dissection13/3/2019 date_of_dissection17/1/2019 CountsPerSample
# 18           1            0              0                           0                           0        19881617
# 19           1            0              0                           0                           0        17151268
# 20           1            0              0                           0                           0        11540367

v <- voom(HC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_Cov_RNAConcExtractGen_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD date_of_dissection13/3/2019 date_of_dissection17/1/2019 CountsPerSample
# -1        3617            1              0                        1245                           0               0
# 0         1068        17627          17629                       15383                       17629           17629
# 1        12944            1              0                        1001                           0               0

#Oh good, it looks like I may be able to get rid of CountsPerSample as co-variate... and maybe also dissection 17/1 vs. 16/1
#I'm guessing that CountsPerSample was partially an effect of that weird 13/3 dissection day (which had high CountsPerSample)

toptable(efit, coef=2)
#Only 1 gene with FDR<0.1 for EE. I wonder why decideTests says two (?) - ah, the topTable summary isn't ordered by pval/FDR, so it is probably just no in the top10 by this definition
#ENSRNOG00000015354
#Aox1: one of the genes we found to be regulated by testosterone in the NACC

toptable(efit, coef=3)
#Nothing even close for SD - probably because everything from one of the dates of dissection couldn't be used to produce an estimate. :(

toptable(efit, coef=5)
#Nothing even close to sig (FDR>0.33)

toptable(efit, coef=6)
#3 genes are FDR<0.10 for CountsPerSample
#ENSRNOG00000023943: Gbp1
#ENSRNOG00000039568: Mgat4e
#ENSRNOG00000048545: LOC100910483

write.fit(efit, adjust="BH", file="Limma_results_Model_onlyEESD_Dissect_CountsPerSample.txt")

write.csv(HC_RNASeq_Annotation_noLowHits, "HC_RNASeq_Annotation_noLowHits.csv")

#################

#Alright, let's toss out Counts Per Sample as a co-variate, since it seems like a weak predictor of variation in the dataset

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/JustEESD_CovDissection")

design <- model.matrix(~Enrichment+SocialDefeat+date_of_dissection, data=Sample_MetaData_HC_NoOutliers_Ordered)

design
# (Intercept) EnrichmentEE SocialDefeatSD date_of_dissection13/3/2019 date_of_dissection17/1/2019
# 18           1            0              0                           0                           0
# 19           1            0              0                           0                           0
# 20           1            0              0                           0                           0

v <- voom(HC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_Cov_RNAConcExtractGen_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

#     (Intercept) EnrichmentEE SocialDefeatSD date_of_dissection13/3/2019 date_of_dissection17/1/2019
# -1        3702            2              0                        2335                           0
# 0          851        17623          17629                       12675                       17629
# 1        13076            4              0                        2619                           0

toptable(efit, coef=2)
#9 genes with FDR<0.1 for EE. I wonder why decideTests says two (?)
#ENSRNOG00000061359: this is a huge effect but it is a pseudogene. Psme1-ps1 hmmm...says logFC almost -5, seems suspicious. 
#ENSRNOG00000000281: Prodh1
#ENSRNOG00000015354: Aox1: one of the genes we found to be regulated by testosterone in the NACC
#ENSRNOG00000010597: Slc5a7
#... I'll look at the rest later.

toptable(efit, coef=3)
#Nothing even close for SD (FDR>0.68) - probably because everything from one of the dates of dissection couldn't be used to produce an estimate. :(

toptable(efit, coef=5)
#Nothing close for the 17/1 dissection day. Throw out? It's a little awkward to justify but the sample is soooo dinky meh.

write.fit(efit, adjust="BH", file="Limma_results_Model_onlyEESD_Dissect.txt")

write.csv(HC_RNASeq_Annotation_noLowHits, "HC_RNASeq_Annotation_noLowHits.csv")

##############

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/EEbySD_CovDissection")

design <- model.matrix(~Enrichment*SocialDefeat+date_of_dissection, data=Sample_MetaData_HC_NoOutliers_Ordered)

design
#   (Intercept) EnrichmentEE SocialDefeatSD date_of_dissection13/3/2019 date_of_dissection17/1/2019 EnrichmentEE:SocialDefeatSD
# 18           1            0              0                           0                           0                           0
# 19           1            0              0                           0                           0                           0
# 20           1            0              0                           0                           0                           0

v <- voom(HC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_Cov_RNAConcExtractGen_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD date_of_dissection13/3/2019 date_of_dissection17/1/2019 EnrichmentEE:SocialDefeatSD
# -1        3589            1              0                        2209                           0                           0
# 0         1032        17628          17629                       13008                       17629                       17629
# 1        13008            0              0                        2412                           0                           0

#yeah, we are underpowered for this.

toptable(efit, coef=2)
#Only 1 gene with FDR<0.1 for EE. 
#ENSRNOG00000015354
#Aox1: one of the genes we found to be regulated by testosterone in the NACC

toptable(efit, coef=3)
#Nothing even close for SD (FDR>0.99) - probably because everything from one of the dates of dissection couldn't be used to produce an estimate. :(

toptable(efit, coef=6)
#Nothing even close for SD*EE (FDR>0.99) 

write.fit(efit, adjust="BH", file="Limma_results_Model_EEbySD_Dissect.txt")


##############################

#Simplest model

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/EEbySD")

design <- model.matrix(~Enrichment*SocialDefeat, data=Sample_MetaData_HC_NoOutliers_Ordered)

design
#       (Intercept) EnrichmentEE SocialDefeatSD EnrichmentEE:SocialDefeatSD
# 18           1            0              0                           0
# 19           1            0              0                           0
# 20           1            0              0                           0

v <- voom(HC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_Cov_RNAConcExtractGen_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD EnrichmentEE:SocialDefeatSD
# -1        3644            1              0                           0
# 0         1041        17628          17629                       17629
# 1        12944            0              0                           0


toptable(efit, coef=2)
#I'm not sure what the one gene is related to EE
toptable(efit, coef=3)
#Nothing even close for SD (FDR>0.36) - probably because everything from one of the dates of dissection couldn't be used to produce an estimate. :(

toptable(efit, coef=4)
#Nothing even close for SD*EE (FDR>0.99) 

write.fit(efit, adjust="BH", file="Limma_results_Model_EEbySD.txt")

###################

#Simplest model (even more so)

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/JustEESD")

design <- model.matrix(~Enrichment+SocialDefeat, data=Sample_MetaData_HC_NoOutliers_Ordered)

design
#     (Intercept) EnrichmentEE SocialDefeatSD
# 18           1            0              0
# 19           1            0              0
# 20           1            0              0
# 21           1            0              0

v <- voom(HC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_Cov_RNAConcExtractGen_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD
# -1        3741            2              0
# 0          859        17623          17629
# 1        13029            4              0

toptable(efit, coef=2)
#Lots of genes related to EE with FDR<0.10

toptable(efit, coef=3)
#Closest for SD is FDR=0.15

write.fit(efit, adjust="BH", file="Limma_results_Model_JustEESD.txt")

#####################

head(HC_RNASeq_Log2_Annotated$gene_symbol)

GroupingVariable<-Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group
table(GroupingVariable)
levels(GroupingVariable)
GroupingVariable_Factor<-factor(GroupingVariable, levels=c("bLR NIL", "bLR EE", "bLR NIL + SD", "bLR EE + SD"))
levels(GroupingVariable_Factor)

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Gpc1"),]

pdf("Gpc1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Gpc1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GroupingVariable<-Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group
table(GroupingVariable)
levels(GroupingVariable)
GroupingVariable_Factor<-factor(GroupingVariable, levels=c("bLR NIL", "bLR EE", "bLR NIL + SD", "bLR EE + SD"))
levels(GroupingVariable_Factor)

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Sdc4"),]

pdf("Sdc4_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Sdc4 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

#I also outputted a lot of the top genes for the NACC since there seems to be a lot of overlap:

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhb6"),]

pdf("Pcdhb6_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhb6 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="AC120486.2"),]

pdf("AC120486_2_Prettier_vs_Group_v3.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="AC120486.2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhga2"),]

pdf("Pcdhga2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhga2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhb5"),]

pdf("Pcdhb5_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhb5 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhb8"),]

pdf("Pcdhb8_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhb8 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Dele1"),]

pdf("Dele1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Dele1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhb7"),]

pdf("Pcdhb7_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhb7 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="RT1-N2"),]

pdf("RT1-N2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="RT1-N2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="RT1-CE4"),]

pdf("RT1-CE4_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="RT1-CE4 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="RT1-CE5"),]

pdf("RT1-CE5_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="RT1-CE5 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Slc26a8"),]

pdf("Slc26a8_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Slc26a8 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Tmtc1"),]

pdf("Tmtc1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Tmtc1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Scn11a"),]

pdf("Scn11a_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Scn11a Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Dhrs2"),]

pdf("Dhrs2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Dhrs2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="AABR07027137.1"),]

pdf("AABR07027137.1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="AABR07027137.1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Megf8"),]

pdf("Megf8_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Megf8 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Scn11a"),]

pdf("Scn11a_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Scn11a Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Pcdh17"),]

pdf("Pcdh17_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdh17 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Abca12"),]

pdf("Abca12_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Abca12 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


#Top HRLR DE genes that also show effects of SD and EE 

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Tmem144"),]

pdf("Tmem144_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Tmem144 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhga1"),]

pdf("Pcdhga1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhga1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Etv4"),]

pdf("Etv4_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Etv4 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()
#I wonder why that doesn't seem to show any effects?


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Uhrf1"),]

pdf("Uhrf1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Uhrf1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Tdg"),]

pdf("Tdg_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Tdg Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Tubg1"),]

pdf("Tubg1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Tubg1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Slc19a3"),]

pdf("Slc19a3_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Slc19a3 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Prss55"),]

pdf("Prss55_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Prss55 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()



#For NIDA-related stuff:

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Oprm1"),]

pdf("Oprm1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Oprm1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()



GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Hba-a2"),][1,]

pdf("Hba-a2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Hba-a2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Hba-a2"),][2,]

pdf("Hba-a2v2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Hba-a2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Hbb"),]

pdf("Hbb_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Hbb Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Htra1"),]

pdf("Htra1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Htra1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Fos"),]

pdf("Fos_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Fos Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


#figure out which of the top DE results are also DE for HR/LR:
setwd("~/Documents/Microarray Gen/HRLR/ThesisMetaAnalysisOutput")

MetaAnalysesResults_AllAges<-read.csv("MetaAnalysesResults_AllAges.csv", header=T, stringsAsFactors = F)

str(MetaAnalysesResults_AllAges)
#Stupid date-gene problem. I'm going to ignore it for now - for a final version of these results we're going to have to fix that. Again.
colnames(MetaAnalysesResults_AllAges)[2]<-"gene_symbol"

HC_EffectSummaryVsMetaAnalysis<-join(HC_EffectSummary, MetaAnalysesResults_AllAges, by="gene_symbol", type="left")

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples")

write.csv(HC_EffectSummaryVsMetaAnalysis, "HC_EffectSummaryVsMetaAnalysis.csv")


#The top gene for Environmental enrichment effects is actually a gene that shows huge HR vs. LR differences - very cool.  
#There are some other cool findings here too.

################################################

#Plot some of the top results to confirm direction of effect, etc.
#I should also do a mediocre attempt to make sure that the effects aren't due to batch confounds.

head(HC_RNASeq_Log2_Annotated)
colnames(HC_RNASeq_Log2_Annotated)

head(HC_RNASeq_Log2_Annotated$gene_symbol)
sum(which(HC_RNASeq_Log2_Annotated$gene_symbol=="Pcmt1"))
#[1] 4

sum(which(HC_RNASeq_Log2_Annotated$gene_symbol=="Sgk1"))
#[1] 97

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Sgk1"),]
GroupingVariable<-Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group

pdf("Sgk1_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Sgk1 Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()
#Huh. Weird plot.

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Olig2"),]

pdf("Olig2_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Olig2 Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Mobp"),]

pdf("Mobp_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Mobp Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Rarres2"),]

pdf("Rarres2_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Rarres2 Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Slc39a12"),]

pdf("Slc39a12_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Slc39a12 Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()



GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Pdyn"),]

pdf("Pdyn_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Pdyn Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Fgfr1"),]

pdf("Fgfr1_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Fgfr1 Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()

table(Sample_MetaData_HC_NoOutliers_Ordered$Treatment_group)
# LR EE  LR EE + SD      LR Nil LR Nil + SD 
# 5           7           5           8 
#So we are best powered to detect an effect of EE in SD animals.


GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Th"),]

pdf("Th_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Th Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Ghr"),]

pdf("Ghr_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Ghr Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Prlhr"),]

pdf("Prlhr_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Prlhr Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()



GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Hba-a2"),][1,]

pdf("Hba-a2_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Hba-a2 Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()

GeneY<-HC_RNASeq_Log2_Filtered[which(HC_RNASeq_Log2_Annotated$gene_symbol=="Hba-a2"),][2,]

pdf("Hba-a2_v2_vs_Group.pdf", width=5, height=7)
boxplot(GeneY~GroupingVariable, ylab="Hba-a2 Cpm", col="grey")
stripchart(GeneY~GroupingVariable, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = 'red')
dev.off()

##########################

library(fgsea)

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/FGSEA")
HC_EESD_Combined<-read.csv("HC_Limma_results_Model_EESDCombined_Dissect_DateGenesFixed.csv", header=TRUE, stringsAsFactors = FALSE)
dim(HC_EESD_Combined)
#[1] 17629    55

colnames(HC_EESD_Combined)

sum(is.na(HC_EESD_Combined$gene_symbol))
#[1] 321

HC_EESD_Combined_noNA<-HC_EESD_Combined[is.na(HC_EESD_Combined$gene_symbol)==FALSE,]
dim(HC_EESD_Combined_noNA)
#[1] 17308    55

sum(duplicated(HC_EESD_Combined_noNA$gene_symbol))
#[1] 156

#MainEffectModel:
HC_EE_Betas_forGSEA<-tapply(X=HC_EESD_Combined_noNA$Coef.EnrichmentEE, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_EE_Betas_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

length(HC_EE_Betas_forGSEA)
#[1] 17152

HC_EE_Betas_forGSEARanked<-HC_EE_Betas_forGSEA[order(HC_EE_Betas_forGSEA)]
head(HC_EE_Betas_forGSEARanked)
# AABR07028009.1 AABR07044362.6 AABR07025819.1       Slc25a54          Otop3        Dnajc5b 
# -4.966         -2.802         -2.610         -1.946         -1.889         -1.839 

HC_SD_Betas_forGSEA<-tapply(X=HC_EESD_Combined_noNA$Coef.SocialDefeatSD, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_SD_Betas_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

HC_SD_Betas_forGSEARanked<-HC_SD_Betas_forGSEA[order(HC_SD_Betas_forGSEA)]
head(HC_SD_Betas_forGSEARanked)
# AABR07044362.6           Gbx2     AC120486.2           Ebf3          Foxb1           Pax7 
# -2.341         -2.097         -2.040         -1.862         -1.760         -1.723 

#InteractionModel:
HC_EEinteract_Betas_forGSEA<-tapply(X=HC_EESD_Combined_noNA$Coef.EnrichmentEE.1, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_EEinteract_Betas_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

HC_EEinteract_Betas_forGSEARanked<-HC_EEinteract_Betas_forGSEA[order(HC_EEinteract_Betas_forGSEA)]
head(HC_EEinteract_Betas_forGSEARanked)
# AABR07044362.6 AABR07028009.1 AABR07044302.1 AABR07044362.1          Khdc1 AABR07071891.2 
# -4.755         -4.705         -2.747         -2.434         -2.378         -2.364 

HC_SDinteract_Betas_forGSEA<-tapply(X=HC_EESD_Combined_noNA$Coef.SocialDefeatSD.1, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_SDinteract_Betas_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

HC_SDinteract_Betas_forGSEARanked<-HC_SDinteract_Betas_forGSEA[order(HC_SDinteract_Betas_forGSEA)]
head(HC_SDinteract_Betas_forGSEARanked)
# AABR07044362.6           Gldn     AC120486.2 AABR07044362.7 AABR07044302.1         Cyp4f5 
# -4.376         -4.262         -3.192         -2.826         -2.496         -2.282 


HC_EEbySD_Betas_forGSEA<-tapply(X=HC_EESD_Combined_noNA$Coef.EnrichmentEE.SocialDefeatSD, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_EEbySD_Betas_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

HC_EEbySD_Betas_forGSEARanked<-HC_EEbySD_Betas_forGSEA[order(HC_EEbySD_Betas_forGSEA)]
head(HC_EEbySD_Betas_forGSEARanked)
# Fam111a    Cnmd  Prss51   Erp27   Dmbx1 Il13ra2 
# -4.081  -3.963  -3.914  -3.800  -3.683  -3.642 

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("c5withBrainCellTypesFunction.all.v7.3.symbols_RatOrtholog.gmt.txt")
str(GMT_ForRats)

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/FGSEA")

temp1<-fgsea(GMT_ForRats, HC_EE_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "HC_EE_MainEffectsModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, HC_SD_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "HC_SD_MainEffectsModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, HC_EEinteract_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "HC_EE_InteractionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, HC_SDinteract_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "HC_SD_InteractionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction.csv")

rm(temp1)

#There is a pretty clear negative correlation between the extremeness of the log2FC and A (average log2 expression) that I think might be biasing the results towards pathways with a large number of low-expressed genes.
#I'm going to try running a version using t-statistics to rank genes instead of log2FC:

HC_EE_Tstat_forGSEA<-tapply(X=HC_EESD_Combined_noNA$t.EnrichmentEE, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_EE_Tstat_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

length(HC_EE_Tstat_forGSEA)
#[1] 17152

HC_EE_Tstat_forGSEARanked<-HC_EE_Tstat_forGSEA[order(HC_EE_Tstat_forGSEA)]
head(HC_EE_Tstat_forGSEARanked)
# AABR07028009.1          Otop2          Cryl1         Retsat         Plcxd2          Otop3 
# -11.72          -6.33          -5.20          -4.95          -4.80          -4.74 


HC_SD_Tstat_forGSEA<-tapply(X=HC_EESD_Combined_noNA$t.SocialDefeatSD, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_SD_Tstat_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

HC_SD_Tstat_forGSEARanked<-HC_SD_Tstat_forGSEA[order(HC_SD_Tstat_forGSEA)]
head(HC_SD_Tstat_forGSEARanked)
# Sgk1        Glra2 LOC100363502    LOC498154        Tfcp2       Pex11a 
# -4.34        -3.99        -3.91        -3.77        -3.74        -3.65 

HC_EEinteract_Tstat_forGSEA<-tapply(X=HC_EESD_Combined_noNA$t.EnrichmentEE.1, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_EEinteract_Tstat_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

HC_EEinteract_Tstat_forGSEARanked<-HC_EEinteract_Tstat_forGSEA[order(HC_EEinteract_Tstat_forGSEA)]
head(HC_EEinteract_Tstat_forGSEARanked)
# AABR07028009.1          Dele1        Pcdhga2        Pcdhga4         Pcdhb5         Pcdhb6 
# -6.93          -5.05          -4.91          -4.63          -4.34          -3.87 

HC_SDinteract_Tstat_forGSEA<-tapply(X=HC_EESD_Combined_noNA$t.SocialDefeatSD.1, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_SDinteract_Tstat_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

HC_SDinteract_Tstat_forGSEARanked<-HC_SDinteract_Tstat_forGSEA[order(HC_SDinteract_Tstat_forGSEA)]
head(HC_SDinteract_Tstat_forGSEARanked)
# Sgk1           Gldn         Trim69 AABR07044362.7 AABR07004404.1           Pdk4 
# -4.23          -4.02          -3.71          -3.63          -3.59          -3.57 

HC_EEbySD_Tstat_forGSEA<-tapply(X=HC_EESD_Combined_noNA$t.EnrichmentEE.SocialDefeatSD, INDEX=HC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(HC_EEbySD_Tstat_forGSEA)<-names(table(HC_EESD_Combined_noNA$gene_symbol))

HC_EEbySD_Tstat_forGSEARanked<-HC_EEbySD_Tstat_forGSEA[order(HC_EEbySD_Tstat_forGSEA)]
head(HC_EEbySD_Tstat_forGSEARanked)
# Hoxc4 LOC360919    Ccdc77   Il13ra2    Hmgn5b    Stk32a 
# -4.02     -3.99     -3.97     -3.86     -3.78     -3.72 

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/FGSEA")

temp1<-fgsea(GMT_ForRats, HC_EE_Tstat_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "HC_EE_MainEffectsModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction_Tstat.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, HC_SD_Tstat_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "HC_SD_MainEffectsModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction_Tstat.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, HC_EEinteract_Tstat_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "HC_EE_InteractionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction_Tstat.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, HC_SDinteract_Tstat_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "HC_SD_InteractionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction_Tstat.csv")

rm(temp1)

#################

#Pulling out gene sets containing the top genes for EE:

GeneSetsWithEEGenesAndInteractionModel_HC<-lapply(GMT_ForRats, function(y) sum(y%in%c("Otop2", "AABR07028009.1", "Prodh1", "Slc5a7", "Aox1", "Pld1", "Cryl1", "Arsb", "Retsat", "H3f3c", "Pid1", "Plcxd2")))
str(GeneSetsWithEEGenesAndInteractionModel_HC)
#Almost all of the top genes for EE from the NACC are also top genes in the interaction model in the HC, but not FDR<0.10
setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/FGSEA")

write.csv(simplify2array(GeneSetsWithEEGenesAndInteractionModel_HC), "GeneSetsWithEEGenesAndInteractionModel_HC.csv")

EE_Otop2<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Otop2"))))
EE_AABR07028009<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("AABR07028009.1"))))
EE_Prodh1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Prodh1"))))
EE_Slc5a7<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Slc5a7"))))
EE_Aox1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Aox1"))))
EE_Pld1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pld1"))))
EE_Cryl1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Cryl1"))))
EE_Arsb<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Arsb"))))
EE_Retsat<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Retsat"))))
EE_H3f3c<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("H3f3c"))))
EE_Pid1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pid1"))))
EE_Plcxd2<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Plcxd2"))))

write.csv(data.frame(EE_Otop2,EE_AABR07028009,EE_Prodh1,EE_Slc5a7,EE_Aox1,EE_Pld1,EE_Cryl1,EE_Arsb,EE_Retsat,EE_H3f3c,EE_Pid1,EE_Plcxd2),"GeneSetsWithEEGenes_wInteractionModel_ByGene_HC.csv")



