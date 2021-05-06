#Angela's bLR social defeat and environmental enrichment RNA-Seq data from the NACC
#Megan Hagenauer
#9/2019-4/2021

#################

#Rstudio (v.1.0.153, R v. 3.4.1)

#######################

#Overview of analysis versions:

## Fan ran a simple version of the differential expression analysis within Ahub (DESeq2) that didn't control for confounds and didn't use a full EE*SD interaction model: August 2019
## My first quick analysis of Fan's preprocessed output to check for potential confounds: September 8, 2019
## Updated September 10, 2020 for full analysis. Note: At this point, the open field and USV data within the metadata dataset are both a work in progress and any analyses including them will need to be re-done.
## Updated again September 28, 2020 to add in additional normalization recommended by limma (TMM)
## Updated again Jan 29 2021 with additional metadata (MBNI measured RNAConc and 260/280) and updated behavioral/hormonal data (for correlations with gene expression). Also note that RNA Purity (260/280) has been misnamed as quality/integrity throughout this analysis.

############################

#The initial preprocessing was performed using the MBNI Analysis Hub (https://ahub.mbni.org) by Dr. Fan Meng in late August 2019 (the BAM files are labeled 20190830).

#Preprocessing methodological details (from Fan):
# - Alignment algorithm (w/version & any non-default parameters): STAR 2.7.0f + Encode parameter with unique hit only
# - Genome assembly: Rnor6 
# - Algorithm for read summarization (w/version & any non-default parameters): featureCounts v1.6.4, (not 100% sure)
# - Gene annotation source (w/version): Ensembl 96

##############################

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

###########################

#Initial basic preprocessing and quality control:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers")

#Reading in Fan's gene-level summary output for the NucleusAccumbens data, with 2 outliers removed based on the PC1/PC2 plot:

NACC_RNASeq<-read.delim("NAcc_remove2samples_gene_featureCounts_counts.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)

str(NACC_RNASeq)
# 'data.frame':	32883 obs. of  47 variables:
#   $ Sample_125991: int  0 0 0 0 0 0 3 0 0 0 ...
# $ Sample_125992: int  0 0 0 0 0 4 2 0 0 2 ...
# $ Sample_125993: int  0 0 0 0 1 0 1 0 0 0 ...
# $ Sample_125994: int  0 0 0 0 0 0 6 0 0 0 ...
# $ Sample_125995: int  0 0 0 0 0 0 1 0 0 0 ...
# $ Sample_125996: int  0 0 0 0 0 0 3 0 0 0 ...

#This data is just basic counts, so we will need to run some transformations to make graphs & run analyses.

#Double-checking whether a log transform has already been performed:
max(NACC_RNASeq)
[1] 15900751
#Nope.


NACC_RNASeq_Matrix<-as.matrix(NACC_RNASeq)
str(NACC_RNASeq_Matrix)
# int [1:32883, 1:47] 0 0 0 0 0 0 3 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:47] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...

NACC_RNASeq_Annotation<-row.names(NACC_RNASeq)

#Counts per sample (library size):
NACC_RNASeq_CountsPerSample<-apply(NACC_RNASeq_Matrix, 2, sum)
head(NACC_RNASeq_CountsPerSample)
# Sample_125991 Sample_125992 Sample_125993 Sample_125994 Sample_125995 Sample_125996 
# 29654770      53057732      32602858      30884626      28007677      30692689 

write.csv(NACC_RNASeq_CountsPerSample, "NACC_RNASeq_CountsPerSample.csv")

pdf("Histogram_CountsPerSample.pdf", width=5, height=5)
hist(NACC_RNASeq_CountsPerSample, breaks=20, xlab="Library Size: Counts per Sample", main="")
dev.off()

summary(NACC_RNASeq_CountsPerSample)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1568502 22667425 29589175 28183629 33080397 53057732 
#Wow - that min is pretty low - 1.5 million! In general the samples tend to fall between 22-33 million.

#Which samples have the really low total counts?

NACC_RNASeq_CountsPerSample[NACC_RNASeq_CountsPerSample<20000000]
# Sample_126010 Sample_126011 Sample_126012 Sample_126015 Sample_126016 Sample_126018 
# 17114093      18055144      12796932       6370726       5458191       1568502 

#These samples all have less than 10 million total counts:
#Sample_126015 Sample_126016 Sample_126018 


#Running some basic filtering:

table(rowSums(NACC_RNASeq_Matrix==0))
# 0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25 
# 13674   970   444   279   266   167   170   163   140   115   109   113   116   113   103   125   100    91   115   101    98   105   110    96   118   108 
# 26    27    28    29    30    31    32    33    34    35    36    37    38    39    40    41    42    43    44    45    46    47 
# 110   123   124   113   134   157   141   155   159   162   194   210   194   258   284   314   376   515   568   867  1465  8151 
#There are 13674 rows with 0's for all subjects!  Get rid of them.
#Also: I should align my analysis with what Fan did.  He says that he filtered out any genes with average counts <1. I think he means actual counts (not CPM), because he didn't calculate CPM (or something similar) for running DE, he just used DESeq2.


#Filtering out rows of data with average counts<1
keep.exprs<-(apply(NACC_RNASeq_Matrix, 1, mean)>1)
table(keep.exprs)
# keep.exprs
# FALSE  TRUE 
# 15118 17765 
#That looks like the # of genes included in Fan's analysis - good.

NACC_RNASeq_Matrix_noLowHits<-NACC_RNASeq_Matrix[keep.exprs,]
dim(NACC_RNASeq_Matrix_noLowHits)
#[1] 17765    47
NACC_RNASeq_Annotation_noLowHits<-NACC_RNASeq_Annotation[keep.exprs]

#TMM normalization to correct for differences in estimated RNA production levels (recommended in limma manual) (note: I did not due this for my first round of analysis - this is new as of 09-28-20)
head(NACC_RNASeq_Matrix_noLowHits)

# Sample_125991 Sample_125992 Sample_125993 Sample_125994 Sample_125995 Sample_125996
# ENSRNOG00000061316             3             2             1             6             1             3
# ENSRNOG00000029897            12            22            17            28            25            23
# ENSRNOG00000014303          9824         18820         12184         10702          9234          9007
# ENSRNOG00000014330          2734          6542          4065          3268          3077          3161

NACC_RNASeq_Log2_Filtered <- cpm(NACC_RNASeq_Matrix_noLowHits, log=TRUE)

head(NACC_RNASeq_Log2_Filtered)
# Sample_125991 Sample_125992 Sample_125993 Sample_125994 Sample_125995 Sample_125996
# ENSRNOG00000061316     -3.183898     -4.424512    -4.6603770    -2.2993610    -4.4875358    -3.2294293
# ENSRNOG00000029897     -1.273874     -1.239428    -0.9150568    -0.1273352    -0.1495376    -0.3991739
# ENSRNOG00000014303      8.372013      8.470618     8.5458783     8.4368926     8.3651190     8.1971618

#As comparison:

NACC_cpm_Filtered <- cpm(NACC_RNASeq_Matrix_noLowHits)
str(NACC_cpm_Filtered)
# 
# # num [1:17765, 1:47] 0.101 0.405 331.295 92.199 5.16 ...
# # - attr(*, "dimnames")=List of 2
# # ..$ : chr [1:17765] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# # ..$ : chr [1:47] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...

NACC_cpm <- cpm(NACC_RNASeq_Matrix)
NACC_cpm_Log2 <- cpm(NACC_RNASeq_Matrix, log=TRUE)
# 
# str(NACC_cpm)
# # num [1:32883, 1:47] 0 0 0 0 0 ...
# # - attr(*, "dimnames")=List of 2
# # ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# # ..$ : chr [1:47] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
# 
# str(NACC_cpm_Log2)
# # num [1:32883, 1:47] -6.82 -6.82 -6.82 -6.82 -6.82 ...
# # - attr(*, "dimnames")=List of 2
# # ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# # ..$ : chr [1:47] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
# 
# max(NACC_cpm_Log2)
# #[1] 19.49466


#Let's compare the distributions before and after filtering:

library(RColorBrewer)

png("NACC_logCPM_Distribution_BeforeAndAfterFiltering.png")
nsamples <- ncol(NACC_RNASeq_Matrix)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(NACC_cpm_Log2[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(NACC_cpm_Log2[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
plot(density(NACC_RNASeq_Log2_Filtered [,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(NACC_RNASeq_Log2_Filtered [,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
dev.off()

pdf("NACC_boxplot_logCPM_afterFiltering.pdf", width=10, height=5)
boxplot(NACC_RNASeq_Log2_Filtered , las=2, col=col, main="")
title(main="Unnormalised data",ylab="Log-cpm")
dev.off()
#The boxes for general distribution of cpm data across all genes look very similar across all samples 
#But there is one sample with a strikingly low level of expression after normalization  (!) - 126017
#e.g., 
#"Sample_126017":
summary(NACC_RNASeq_Log2_Filtered[,27])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -6.817  -1.798   1.685   0.745   3.618  19.495 

#For comparison
summary(NACC_RNASeq_Log2_Filtered[,26])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -6.817  -0.431   3.484   2.280   5.465  14.360 
summary(NACC_RNASeq_Log2_Filtered[,28])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -6.8167 -0.6294  3.5217  1.9921  5.4803 14.3093 


#We will probably just need to remove it as an outlier:

names(NACC_RNASeq_CountsPerSample[27])
#[1] "Sample_126017"

NACC_RNASeq_Log2_Filtered<-NACC_RNASeq_Log2_Filtered[,-27]

NACC_RNASeq_Matrix_noLowHits<-NACC_RNASeq_Matrix_noLowHits[,-27]

NACC_RNASeq_CountsPerSample<-NACC_RNASeq_CountsPerSample[-27]


################################

#TMM normalization: Scaling to correct for estimated relative RNA production levels

#I didn't originally include this step, but my analysis revealed that library size (NACC_RNASeq_CountsPerSample) was correlated with the principal components of variation in the data, despite cpm normalization. Also, the limma manual recommends it.

dge<-DGEList(counts=NACC_RNASeq_Matrix_noLowHits)
str(dge)

# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:17765, 1:47] 3 12 9824 2734 153 2089 222 381 1800 8 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:17765] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# .. .. .. ..$ : chr [1:47] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
# .. ..$ :'data.frame':	47 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:47] 29653288 53054264 32601268 30882918 28005917 ...
# .. .. ..$ norm.factors: num [1:47] 1 1 1 1 1 1 1 1 1 1 ...

str(dge[[2]])

dge[[2]]$lib.size
#This is the same as NACC_RNASeq_CountsPerSample - useful
#Group is currently empty - do we need it?

head(dge[[1]])
# Sample_125991 Sample_125992 Sample_125993 Sample_125994 Sample_125995 Sample_125996
# ENSRNOG00000061316             3             2             1             6             1             3
# ENSRNOG00000029897            12            22            17            28            25            23
# ENSRNOG00000014303          9824         18820         12184         10702          9234          9007

NACC_RNASeq_dge_TMM<-calcNormFactors(dge, method="TMM")
str(NACC_RNASeq_dge_TMM)

# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:17765, 1:47] 3 12 9824 2734 153 2089 222 381 1800 8 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:17765] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# .. .. .. ..$ : chr [1:47] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
# .. ..$ :'data.frame':	47 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:47] 29653288 53054264 32601268 30882918 28005917 ...
# .. .. ..$ norm.factors: num [1:47] 1.01 1.11 1.05 1.02 1 ...
# 
head(NACC_RNASeq_dge_TMM[[1]])
# Sample_125991 Sample_125992 Sample_125993 Sample_125994 Sample_125995 Sample_125996
# ENSRNOG00000061316             3             2             1             6             1             3
# ENSRNOG00000029897            12            22            17            28            25            23
# ENSRNOG00000014303          9824         18820         12184         10702          9234          9007

#Interesting - the object seems to store the normalization factors but not the normalized data.
#Double-checked - that is what it is supposed to be. Interesting. At what point is the data actually scaled with those normalization factors?

#Let's see if it happens during the cpm function:

NACC_RNASeq_Log2_Filtered_TMM <- cpm(NACC_RNASeq_dge_TMM, log=TRUE)
head(NACC_RNASeq_Log2_Filtered_TMM)

# Sample_125991 Sample_125992 Sample_125993 Sample_125994 Sample_125995 Sample_125996
# ENSRNOG00000061316     -3.149538     -4.500036    -4.6810195    -2.2863866     -4.468329    -3.3346221
# ENSRNOG00000029897     -1.236785     -1.329346    -0.9389092    -0.1135963     -0.123816    -0.5121697
# ENSRNOG00000014303      8.410089      8.378914     8.5217673     8.4508488      8.391178     8.0828886


head(NACC_RNASeq_Log2_Filtered)
# Sample_125991 Sample_125992 Sample_125993 Sample_125994 Sample_125995 Sample_125996
# ENSRNOG00000061316     -3.183898     -4.424512    -4.6603770    -2.2993610    -4.4875358    -3.2294293
# ENSRNOG00000029897     -1.273874     -1.239428    -0.9150568    -0.1273352    -0.1495376    -0.3991739
# ENSRNOG00000014303      8.372013      8.470618     8.5458783     8.4368926     8.3651190     8.1971618

#Yep - it looks like the normalization does influence the final cpm. Good.


pdf("NACC_boxplot_logCPM_afterFilteringTMM.pdf", width=10, height=5)
boxplot(NACC_RNASeq_Log2_Filtered_TMM , las=2, col=col, main="")
title(main="TMM normalised data, outlier removed",ylab="Log-cpm")
dev.off()

#Overwriting our previous object so we don't need to change the rest of the code pipeline:
NACC_RNASeq_Log2_Filtered<-NACC_RNASeq_Log2_Filtered_TMM

###################################

#Getting additional annotation:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress")

NACC_Annotation<-read.delim("Annotation_Fan.txt", sep="\t", header=T, stringsAsFactors = F)
str(NACC_Annotation)
# 'data.frame':	17335 obs. of  4 variables:
#   $ gene_id     : chr  "ENSRNOG00000037206" "ENSRNOG00000010452" "ENSRNOG00000011762" "ENSRNOG00000020388" ...
# $ gene_symbol : chr  "Ccdc77" "LOC100363502" "Elf1" "Inpp5f" ...
# $ gene_biotype: chr  "protein_coding" "protein_coding" "protein_coding" "protein_coding" ...
# $ description : chr  "coiled-coil domain containing 77 [Source:RGD Symbol;Acc:1310710]" "cytochrome c, somatic-like [Source:RGD Symbol;Acc:2322845]" "E74 like ETS transcription factor 1 [Source:RGD Symbol;Acc:620697]" "inositol polyphosphate-5-phosphatase F [Source:RGD Symbol;Acc:1305777]" ...

#For ease, let's put the annotation in the same order as the Count data:

gene_id<-NACC_RNASeq_Annotation_noLowHits

ToJoin<-data.frame(gene_id, NACC_RNASeq_Log2_Filtered, stringsAsFactors = F)
str(ToJoin)

NACC_RNASeq_Log2_Annotated<-join(ToJoin, NACC_Annotation, by="gene_id", type="left")
str(NACC_RNASeq_Log2_Annotated)

write.csv(NACC_RNASeq_Log2_Annotated, "NACC_RNASeq_Log2_Annotated.csv")

##################################

#Meta-data:


#New version of the analysis (from 01/2021)
setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/ReDone_OpenFieldData_20201008")
Sample_BehavHormonalData<-read.csv("HRLR_EE_Stress_AllBehavData_forR_withNewCORTOxytIL6_SI_OFSDScoresFixed_FixedFormatIDs_TimeOnTop_forNACC2.csv", header=T, stringsAsFactors = F)

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
# [52] "date_of_dissection"                                                                                             
# [53] "hemisphere_dissected"                                                                                           
# [54] "date_of_RNA_extraction"                                                                                         
# [55] "RNA_conc"                                                                                                       
# [56] "ratio_260_280"                                                                                                  
# [57] "brain_region"                                                                                                   
# [58] "Sequencing_core_sample_ID"  

#Let's filter both data-frames just down to the samples that were included in the RNA-Seq experiment:

table(Sample_MetaData$brain_region)
# Nacc 
# 49

sum(Sample_MetaData$brain_region=="Nacc")
#[1] 49

Sample_MetaData<-Sample_MetaData[Sample_MetaData$brain_region=="Nacc",]
str(Sample_MetaData)
#'data.frame':	49 obs. of  14 variables:

table(Sample_BehavHormonalData$NACC_brain_region)
# Nacc 
# 93   49 
#Hmmm.... looks like I only previously annotated this data-frame with information about which subjects were included in the NAcc experiment. I'll have to fix that before running the hippocampal analysis.

Sample_BehavHormonalData<-Sample_BehavHormonalData[Sample_BehavHormonalData$NACC_brain_region=="Nacc",]
str(Sample_BehavHormonalData)
#'data.frame':	49 obs. of  58 variables:


#To make things easier, I should put the metadata in the same order as the NACC RNASeq data file:
colnames(NACC_RNASeq_Log2_Filtered)
# [1] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" "Sample_125995" "Sample_125996" "Sample_125997" "Sample_125998" "Sample_125999"
# [10] "Sample_126000" "Sample_126001" "Sample_126002" "Sample_126003" "Sample_126004" "Sample_126005" "Sample_126006" "Sample_126007" "Sample_126008"
# [19] "Sample_126009" "Sample_126010" "Sample_126011" "Sample_126012" "Sample_126013" "Sample_126014" "Sample_126015" "Sample_126016" "Sample_126018"
# [28] "Sample_126019" "Sample_126021" "Sample_126022" "Sample_126023" "Sample_126024" "Sample_126025" "Sample_126026" "Sample_126027" "Sample_126028"
# [37] "Sample_126029" "Sample_126030" "Sample_126031" "Sample_126032" "Sample_126033" "Sample_126034" "Sample_126035" "Sample_126036" "Sample_126037"
# [46] "Sample_126039"

#the column names are just in sample ID order.

Sample_MetaData$Sequencing_core_sample_ID
# [1] 126010 126011 126012 126013 126017 126018 126019 126020 126021 125999 126000 126001 126002 126003
# [15] 126004 126005 126006 126028 126029 126030 126031 126032 126033 126034 126035 126036 126037 126038
# [29] 126039 125994 125995 125996 125997 125998 126022 126023 125991 125992 125993 126014 126015 126016
# [43] 126007 126008 126009 126024 126025 126026 126027

Sample_MetaData$SampleID<-paste("Sample_", Sample_MetaData$Sequencing_core_sample_ID, sep="")
Sample_MetaData$SampleID
# [1] "Sample_126010" "Sample_126011" "Sample_126012" "Sample_126013" "Sample_126017" "Sample_126018"
# [7] "Sample_126019" "Sample_126020" "Sample_126021" "Sample_125999" "Sample_126000" "Sample_126001"
# [13] "Sample_126002" "Sample_126003" "Sample_126004" "Sample_126005" "Sample_126006" "Sample_126028"
# [19] "Sample_126029" "Sample_126030" "Sample_126031" "Sample_126032" "Sample_126033" "Sample_126034"
# [25] "Sample_126035" "Sample_126036" "Sample_126037" "Sample_126038" "Sample_126039" "Sample_125994"
# [31] "Sample_125995" "Sample_125996" "Sample_125997" "Sample_125998" "Sample_126022" "Sample_126023"
# [37] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_126014" "Sample_126015" "Sample_126016"
# [43] "Sample_126007" "Sample_126008" "Sample_126009" "Sample_126024" "Sample_126025" "Sample_126026"
# [49] "Sample_126027"

Sample_MetaData_NACC_NoOutliers<-Sample_MetaData[which(Sample_MetaData$SampleID%in%colnames(NACC_RNASeq_Log2_Filtered)), ]
str(Sample_MetaData_NACC_NoOutliers)
# 'data.frame':	46 obs. of  15 variables:
#   $ Sequencing_core_sample_ID: int  126010 126011 126012 126013 126018 126019 126021 125999 126000 126001 ...
# $ Rat_ID                   : chr  "KH L07 A1" "KH L07 A2" "KH L07 B1" "KH L07 B2" ...
# $ Generation               : chr  "F56" "F56" "F56" "F56" ...
# $ Line                     : chr  "bLR" "bLR" "bLR" "bLR" ...
# $ Treatment_group          : chr  "bLR EE" "bLR EE" "bLR EE" "bLR EE" ...
# $ Enrichment               : chr  "EE" "EE" "EE" "EE" ...
# $ SocialDefeat             : chr  "NIL" "NIL" "NIL" "NIL" ...
# $ brain_region             : chr  "Nacc" "Nacc" "Nacc" "Nacc" ...
# $ date_of_sac              : chr  "2-Nov-17" "2-Nov-17" "2-Nov-17" "2-Nov-17" ...
# $ date_of_dissection       : chr  "16/1/2019" "17/1/2019" "16/1/2019" "16/1/2019" ...
# $ hemisphere_dissected     : chr  "right" "left" "right" "left" ...
# $ date_of_RNA_extraction   : chr  "28.2.2019" "26.2.2019" "26.2.2019" "26.2.2019" ...
# $ RNA_conc                 : num  42.7 51.5 56.6 78.1 62.7 27.7 55.4 91.2 40.6 98.9 ...
# $ ratio_260_280            : num  1.98 2.08 2.08 2.06 2.08 1.95 2.09 2.05 1.94 2.06 ...
# $ SampleID                 : chr  "Sample_126010" "Sample_126011" "Sample_126012" "Sample_126013" ...

table(Sample_MetaData_NACC_NoOutliers$brain_region)
# Nacc 
# 46 
#Good.

#But it may still be in the wrong order:
Sample_MetaData_NACC_NoOutliers$Sequencing_core_sample_ID
#Yep.

Sample_MetaData_NACC_NoOutliers_Ordered<-Sample_MetaData_NACC_NoOutliers[order(Sample_MetaData_NACC_NoOutliers$Sequencing_core_sample_ID),]
Sample_MetaData_NACC_NoOutliers_Ordered$Sequencing_core_sample_ID
# [1] 125991 125992 125993 125994 125995 125996 125997 125998 125999 126000 126001 126002 126003 126004 126005 126006 126007 126008 126009 126010 126011 126012
# [23] 126013 126014 126015 126016 126018 126019 126021 126022 126023 126024 126025 126026 126027 126028 126029 126030 126031 126032 126033 126034 126035 126036
# [45] 126037 126039

cbind(colnames(NACC_RNASeq_Log2_Filtered), Sample_MetaData_NACC_NoOutliers_Ordered$SampleID)
#Looks good now.
sum(colnames(NACC_RNASeq_Log2_Filtered)==Sample_MetaData_NACC_NoOutliers_Ordered$SampleID)
#[1] 46
#Yep, they're all in the same order.

cbind(names(NACC_RNASeq_CountsPerSample), Sample_MetaData_NACC_NoOutliers_Ordered$SampleID)
#Same order

Sample_MetaData_NACC_NoOutliers_Ordered$CountsPerSample<-NACC_RNASeq_CountsPerSample

Sample_BehavHormonalData$SampleID<-paste("Sample_", Sample_BehavHormonalData$NACC_Sequencing_core_sample_ID, sep="")
Sample_BehavHormonalData$SampleID
# [1] "Sample_126010" "Sample_126011" "Sample_126012" "Sample_126013" "Sample_126017" "Sample_126018" "Sample_126019" "Sample_126020"
# [9] "Sample_126021" "Sample_125999" "Sample_126000" "Sample_126001" "Sample_126002" "Sample_126003" "Sample_126004" "Sample_126005"
# [17] "Sample_126006" "Sample_126028" "Sample_126029" "Sample_126030" "Sample_126031" "Sample_126032" "Sample_126033" "Sample_126034"
# [25] "Sample_126035" "Sample_126036" "Sample_126037" "Sample_126038" "Sample_126039" "Sample_125994" "Sample_125995" "Sample_125996"
# [33] "Sample_125997" "Sample_125998" "Sample_126022" "Sample_126023" "Sample_125991" "Sample_125992" "Sample_125993" "Sample_126014"
# [41] "Sample_126015" "Sample_126016" "Sample_126007" "Sample_126008" "Sample_126009" "Sample_126024" "Sample_126025" "Sample_126026"
# [49] "Sample_126027"

Sample_BehavHormonalData_NACC_NoOutliers<-Sample_BehavHormonalData[which(Sample_BehavHormonalData$SampleID%in%colnames(NACC_RNASeq_Log2_Filtered)), ]
str(Sample_BehavHormonalData_NACC_NoOutliers)
#'data.frame':	46 obs. of  59 variables:

Sample_BehavHormonalData_NACC_NoOutliers_Ordered<-Sample_BehavHormonalData_NACC_NoOutliers[order(Sample_BehavHormonalData_NACC_NoOutliers$NACC_Sequencing_core_sample_ID),]
Sample_BehavHormonalData_NACC_NoOutliers_Ordered$NACC_Sequencing_core_sample_ID
# [1] 125991 125992 125993 125994 125995 125996 125997 125998 125999 126000 126001 126002 126003 126004 126005 126006 126007 126008
# [19] 126009 126010 126011 126012 126013 126014 126015 126016 126018 126019 126021 126022 126023 126024 126025 126026 126027 126028
# [37] 126029 126030 126031 126032 126033 126034 126035 126036 126037 126039

cbind(colnames(NACC_RNASeq_Log2_Filtered), Sample_BehavHormonalData_NACC_NoOutliers_Ordered$SampleID)
#Looks good now.
sum(colnames(NACC_RNASeq_Log2_Filtered)==Sample_BehavHormonalData_NACC_NoOutliers_Ordered$SampleID)
#[1] 46
#yep

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/AnalysisUpdates_20210129")

write.csv(Sample_BehavHormonalData_NACC_NoOutliers_Ordered, "Sample_BehavHormonalData_NACC_NoOutliers_Ordered.csv")

#Getting the factors set-up properly:
str(Sample_MetaData_NACC_NoOutliers_Ordered)
Sample_MetaData_NACC_NoOutliers_Ordered$Generation<-as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Generation)
levels(Sample_MetaData_NACC_NoOutliers_Ordered$Generation)
#[1] "F49" "F53" "F56"
table(Sample_MetaData_NACC_NoOutliers_Ordered$Generation)
# F49 F53 F56 
# 17   9  20 
#Let's make F56 the reference level
Sample_MetaData_NACC_NoOutliers_Ordered$Generation<-relevel(Sample_MetaData_NACC_NoOutliers_Ordered$Generation, ref="F56")
levels(Sample_MetaData_NACC_NoOutliers_Ordered$Generation)
#[1] "F56" "F49" "F53"

Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment<-as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment)
levels(Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment)
#[1] "EE"  "NIL"
Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment<-relevel(Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment, ref="NIL")
levels(Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment)
#[1] "NIL" "EE" 

Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat<-as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat)
levels(Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat)
#[1] "NIL" "SD" 

#########################################

#Double-checking whether any of the batch-related variables are unevenly distributed by group:

table(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$Generation)
#             F56 F49 F53
# bLR EE         4   6   3
# bLR EE + SD    3   5   3
# bLR NIL        5   2   3
# bLR NIL + SD   8   4   0

#Generation is unevenly distributed. :(
#But maybe not too bad? F49 and F56 are kind-of close...sort-of...

write.csv(table(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$Generation), "NACC_TreatmentGroupByGeneration.csv")

fisher.test(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$Generation)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group and Sample_MetaData_NACC_NoOutliers_Ordered$Generation
# p-value = 0.2514
# alternative hypothesis: two.sided

#But not significantly so...


table(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)

#                 16/1/2019 17/1/2019 21/2/2019
# bLR EE               6         1         6
# bLR EE + SD          4         2         5
# bLR NIL              8         0         2
# bLR NIL + SD         4         4         4

#Hmmm..dissection date may be unevenly distributed across groups.

write.csv(table(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), "NACC_TreatmentGroupByDateOfDissection.csv")

fisher.test(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group and Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection
# p-value = 0.2211
# alternative hypothesis: two.sided

#but not significantly so...

table(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
#               26.2.2019 28.2.2019
# LR EE               9         4
# LR EE + SD          4         7
# LR Nil              6         4
# LR Nil + SD         8         4
#.... and RNA-extraction date too.

write.csv(table(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction), "NACC_TreatmentGroupByDateOfRNAExtraction.csv")


fisher.test(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group and Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction
# p-value = 0.3868
# alternative hypothesis: two.sided

#...But also not significantly so.


#I think they may be partially redundant. Let's check:

table(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
#             26.2.2019 28.2.2019
# 16/1/2019        14         8
# 17/1/2019         5         2
# 21/2/2019         8         9
#Nope. There are 6 distinct batches. Meh.

#I'm still going to make a combined variable for ease of spotting effects later:
Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch<-paste(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
table(Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch)

# 16/1/2019 26.2.2019 16/1/2019 28.2.2019 17/1/2019 26.2.2019 17/1/2019 28.2.2019 21/2/2019 26.2.2019 
# 14                   8                   5                   2                   8 
# 21/2/2019 28.2.2019 
# 9 

table(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch)

#                16/1/2019 26.2.2019 16/1/2019 28.2.2019 17/1/2019 26.2.2019 17/1/2019 28.2.2019
# bLR EE                         5                   1                   1                   0
# bLR EE + SD                    1                   3                   1                   1
# bLR NIL                        4                   4                   0                   0
# bLR NIL + SD                   4                   0                   3                   1
# 
#                     21/2/2019 26.2.2019 21/2/2019 28.2.2019
# bLR EE                         3                   3
# bLR EE + SD                    2                   3
# bLR NIL                        2                   0
# bLR NIL + SD                   1                   3

fisher.test(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group, Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch)

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group and Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch
# p-value = 0.2176
# alternative hypothesis: two.sided



#Out of curiousity, I wonder what things looked like in the original data before Fan removed the outliers:

table(Sample_MetaData$date_of_dissection[Sample_MetaData$brain_region=="Nacc"], Sample_MetaData$date_of_RNA_extraction[Sample_MetaData$brain_region=="Nacc"])
#               26.2.2019 28.2.2019
# 16/1/2019        16         8
# 17/1/2019         5         2
# 21/2/2019         9         9

#Two outliers came from the 16/1 26.2.2019 batch, 1 came from the 21/2 26.2.2019 batch

table(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
#               26.2.2019 28.2.2019
# 16/1/2019        14         8
# 17/1/2019         5         2
# 21/2/2019         8         9

#So the outliers weren't simply a bad batch. Interesting.


#Are these also related to generation?
table(Sample_MetaData_NACC_NoOutliers_Ordered$Generation, Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch)

# 16/1/2019 26.2.2019 16/1/2019 28.2.2019 17/1/2019 26.2.2019 17/1/2019 28.2.2019 21/2/2019 26.2.2019
# F56                   8                   5                   5                   2                   0
# F49                   0                   0                   0                   0                   8
# F53                   6                   3                   0                   0                   0
# 
# 21/2/2019 28.2.2019
# F56                   0
# F49                   9
# F53                   0


table(Sample_MetaData_NACC_NoOutliers_Ordered$Generation, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
#         16/1/2019 17/1/2019 21/2/2019
# F56        13         7         0
# F49         0         0        17
# F53         9         0         0


#Yes, almost exclusively - dissection batch and generation are almost identical, as makes sense.

write.csv(table(Sample_MetaData_NACC_NoOutliers_Ordered$Generation, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), "NACC_GenerationByDateOfDissection.csv")


#Let's see how RNA conc measured at the core relates to RNA conc measured at MBNI:

plot(Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc, Sample_BehavHormonalData_NACC_NoOutliers_Ordered$NACC_RNA_conc)
#exactly 1:1 correlation?
#Looks to me like they didn't re-measure it, they just used Angela's measurements)


anova(lm(RNA_conc~Treatment_group, data=Sample_MetaData_NACC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: RNA_conc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Treatment_group  3  2579.3  859.78  1.2489 0.3042
# Residuals       42 28914.7  688.45  

summary.lm(lm(Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc~Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat))

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -46.161 -17.541  -3.811  14.721  88.065 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                              68.335      5.572  12.263 8.63e-16 ***
#   Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeatSD   -2.474      7.880  -0.314    0.755    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 26.72 on 44 degrees of freedom
# Multiple R-squared:  0.002235,	Adjusted R-squared:  -0.02044 
# F-statistic: 0.09855 on 1 and 44 DF,  p-value: 0.7551

summary.lm(lm(Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc~Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment))
# Call:
#   lm(formula = Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc ~ 
#        Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -40.637 -17.346  -4.188  16.611  81.927 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                            74.473      5.492  13.559   <2e-16 ***
#   Sample_MetaData_NACC_NoOutliers_Ordered$EnrichmentEE  -14.135      7.604  -1.859   0.0697 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 25.76 on 44 degrees of freedom
# Multiple R-squared:  0.07282,	Adjusted R-squared:  0.05175 
# F-statistic: 3.456 on 1 and 44 DF,  p-value: 0.06973

pdf("NAcc_Boxplot_RNAConc_vs_TreatmentGroup.pdf", width=5, height=5)
boxplot(RNA_conc~Treatment_group, data=Sample_MetaData_NACC_NoOutliers_Ordered, cex.axis=0.5, ylab="RNA Concentration")
dev.off()

summary.lm(lm(RNA_conc~Enrichment+SocialDefeat, data=Sample_MetaData_NACC_NoOutliers_Ordered))

# Call:
#   lm(formula = RNA_conc ~ Enrichment + SocialDefeat, data = Sample_MetaData_NACC_NoOutliers_Ordered)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -38.616 -18.560  -2.766  17.081  79.892 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      76.508      6.950  11.008 4.32e-14 ***
#   EnrichmentEE    -14.460      7.700  -1.878   0.0672 .  
# SocialDefeatSD   -3.731      7.693  -0.485   0.6301    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 25.99 on 43 degrees of freedom
# Multiple R-squared:  0.07787,	Adjusted R-squared:  0.03498 
# F-statistic: 1.815 on 2 and 43 DF,  p-value: 0.175

#RNA concentration may be related to enrichment.


anova(lm(ratio_260_280~Treatment_group, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: ratio_260_280
# Df   Sum Sq   Mean Sq F value Pr(>F)
# Treatment_group  3 0.004117 0.0013723  0.5892 0.6255
# Residuals       42 0.097824 0.0023291

#Good, RNA concentration and Quality don't vary with treatment group.

summary.lm(lm(ratio_260_280~Enrichment+SocialDefeat, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.11660 -0.01996  0.00340  0.02993  0.12993 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     2.061080   0.012979 158.797   <2e-16 ***
#   EnrichmentEE   -0.006527   0.014379  -0.454    0.652    
# SocialDefeatSD -0.004481   0.014366  -0.312    0.757    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04853 on 43 degrees of freedom
# Multiple R-squared:  0.006487,	Adjusted R-squared:  -0.03972 
# F-statistic: 0.1404 on 2 and 43 DF,  p-value: 0.8694

anova(lm(CountsPerSample~Treatment_group, data=Sample_MetaData_NACC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# Treatment_group  3 1.5539e+15 5.1797e+14   7.702 0.0003277 ***
#   Residuals       42 2.8246e+15 6.7251e+13                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#But total counts per sample do. Ouch!

pdf("NAcc_Boxplot_CountsPerSample_vs_TreatmentGroup.pdf", width=5, height=5)
boxplot(CountsPerSample~Treatment_group, data=Sample_MetaData_NACC_NoOutliers_Ordered, cex.axis=0.5, ylab="Total Counts Per Sample")
dev.off()

summary.lm(lm(CountsPerSample~Enrichment+SocialDefeat, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -20540527  -4465145   -205265   5089679  18034901 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     35022831    2236054  15.663  < 2e-16 ***
#   EnrichmentEE   -10927086    2477260  -4.411 6.78e-05 ***
#   SocialDefeatSD  -1986716    2474917  -0.803    0.427    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8361000 on 43 degrees of freedom
# Multiple R-squared:  0.3135,	Adjusted R-squared:  0.2815 
# F-statistic: 9.817 on 2 and 43 DF,  p-value: 0.0003078

summary.lm(lm(CountsPerSample~Enrichment*SocialDefeat, data=Sample_MetaData_NACC_NoOutliers_Ordered))

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -18474252  -3729348  -1003178   3575435  20307804 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 32749928    2593287  12.629 6.84e-16 ***
#   EnrichmentEE                -6905797    3449394  -2.002   0.0518 .  
# SocialDefeatSD               2180272    3511327   0.621   0.5380    
# EnrichmentEE:SocialDefeatSD -7981650    4859667  -1.642   0.1080    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8201000 on 42 degrees of freedom
# Multiple R-squared:  0.3549,	Adjusted R-squared:  0.3088 
# F-statistic: 7.702 on 3 and 42 DF,  p-value: 0.0003277

#Is this some sort of artifact due to the order of samples at the sequencing core?

pdf("CountsPerSample_vsSeqCoreID_byEE.pdf", height=5, width=8)
plot(CountsPerSample~Sequencing_core_sample_ID, data=Sample_MetaData_NACC_NoOutliers_Ordered, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment))
dev.off()


#Does RNAconc correlate with purity?

plot(Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc, Sample_BehavHormonalData_NACC_NoOutliers_Ordered$NACC_ratio_260_280)
#Maybe - the samples with the lowest concentration may have either the lowest or highest 260/280, but hard to say, and the there isn't much range in purity (the samples are all very pure - 1.95-2.20)

cor(Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280, Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc) 
#[1] 0.1312447

cor(NACC_RNASeq_CountsPerSample, Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280)
#[1] -0.1359465

cor(NACC_RNASeq_CountsPerSample, Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc)
#[1] 0.1621681

plot(NACC_RNASeq_CountsPerSample~Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc)
#pretty unimpressive - the correlation is mostly driven by the samples with RNA_conc>80

#What about generation?

anova(lm(RNA_conc~Generation, data=Sample_MetaData_NACC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: RNA_conc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Generation  2   633.8  316.91  0.4416 0.6459
# Residuals  43 30860.2  717.68 

anova(lm(ratio_260_280~Generation, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: ratio_260_280
# Df   Sum Sq   Mean Sq F value Pr(>F)
# Generation  2 0.002510 0.0012552  0.5428  0.585
# Residuals  43 0.099431 0.0023123

anova(lm(CountsPerSample~Generation, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value Pr(>F)
# Generation  2 5.0805e+13 2.5402e+13  0.2524 0.7781
# Residuals  43 4.3277e+15 1.0064e+14 

#No relationship with generation.

anova(lm(RNA_conc~DissectionExtractionBatch, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: RNA_conc
# Df  Sum Sq Mean Sq F value Pr(>F)
# DissectionExtractionBatch  5  3904.4  780.87  1.1321 0.3593
# Residuals                 40 27589.7  689.74 

pdf("NAcc_Boxplot_RNA_conc_vs_DissectionExtractionBatch.pdf", width=15, height=5)
boxplot(Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc~Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch, cex.axis=0.5)
dev.off()
#Although it isn't significant, the 2 samples from 17/1/2019 28.2.2019 seem to have noticeably lower RNA concentration.

anova(lm(ratio_260_280~DissectionExtractionBatch, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: ratio_260_280
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# DissectionExtractionBatch  5 0.055721 0.0111443  9.6445 4.255e-06 ***
#   Residuals                 40 0.046220 0.0011555                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#But RNA-quality does very intensely by dissection-extraction batch. Ouch.  Let's take a closer look at that:

pdf("NAcc_Boxplot_ratio_260_280_vs_DissectionExtractionBatch.pdf", width=15, height=5)
boxplot(Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280~Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch, cex.axis=0.5)
dev.off()

#Looks like 17/1/2019 28.2.2019 has much lower RNA quality than the other batches.
#Actually, in general it seems like the 28.2.2019 extractions may have lower RNA Quality. Let's test that

anova(lm(ratio_260_280~date_of_RNA_extraction, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: ratio_260_280
# Df   Sum Sq  Mean Sq F value    Pr(>F)    
# date_of_RNA_extraction  1 0.034832 0.034832  22.838 1.994e-05 ***
#   Residuals              44 0.067109 0.001525                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Yep, definitely.

pdf("NAcc_Boxplot_ratio_260_280_vs_ExtractionBatch.pdf", width=5, height=5)
boxplot(Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction, cex.axis=0.5, ylab="RNA Quality: Ratio 260/280")
dev.off()
#Yep, definitely. But I'm not sure how much this effect size would matter - it is quite small, just consistently different between groups:

summary.lm(lm(lm(ratio_260_280~date_of_RNA_extraction, data=Sample_MetaData_NACC_NoOutliers_Ordered)))

# Call:
#   lm(formula = lm(ratio_260_280 ~ date_of_RNA_extraction, data = Sample_MetaData_NACC_NoOutliers_Ordered))
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.082632 -0.021603  0.001481  0.020453  0.101481 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      2.078519   0.007516 276.548  < 2e-16 ***
#   date_of_RNA_extraction28.2.2019 -0.055887   0.011695  -4.779 1.99e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03905 on 44 degrees of freedom
# Multiple R-squared:  0.3417,	Adjusted R-squared:  0.3267 
# F-statistic: 22.84 on 1 and 44 DF,  p-value: 1.994e-05

anova(lm(CountsPerSample~DissectionExtractionBatch, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value Pr(>F)
# DissectionExtractionBatch  5 6.6187e+14 1.3237e+14  1.4247 0.2364
# Residuals                 40 3.7166e+15 9.2915e+13

pdf("NAcc_Boxplot_CountsPerSample_vs_DissectionExtractionBatch.pdf", width=8, height=5)
boxplot(CountsPerSample~DissectionExtractionBatch, data=Sample_MetaData_NACC_NoOutliers_Ordered, cex.axis=0.5)
dev.off()

#Although it isn't significant, the samples from 16/1/2019 28.2.2019 seem to have noticeably lower CountsPerSample. There are a few with lower Counts per sample in 16/1/2019 26.2.2019 - maybe a dissection batch issue?

pdf("NAcc_Boxplot_CountsPerSample_vs_DissectionBatch.pdf", width=5, height=5)
boxplot(CountsPerSample~date_of_dissection, data=Sample_MetaData_NACC_NoOutliers_Ordered, ylab="Library Size: Counts per Sample", cex.axis=0.5)
dev.off()

anova(lm(CountsPerSample~date_of_dissection, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value Pr(>F)
# date_of_dissection  2 3.1929e+14 1.5965e+14  1.6912 0.1963
# Residuals          43 4.0592e+15 9.4399e+13

#Not sig.

anova(lm(CountsPerSample~date_of_dissection+Treatment_group, data=Sample_MetaData_NACC_NoOutliers_Ordered))

# Analysis of Variance Table
# 
# Response: CountsPerSample
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# date_of_dissection  2 3.1929e+14 1.5965e+14  2.6693   0.08163 .  
# Treatment_group     3 1.6668e+15 5.5561e+14  9.2897 8.702e-05 ***
#   Residuals          40 2.3924e+15 5.9809e+13                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#But almost sig if you also control for the relationship with treatment group...

summary.lm(lm(CountsPerSample~date_of_dissection+Enrichment+SocialDefeat, data=Sample_MetaData_NACC_NoOutliers_Ordered))
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -15432530  -4794919  -1245205   4972271  20031426 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  33026306    2200940  15.006  < 2e-16 ***
#   date_of_dissection17/1/2019   9490029    3574175   2.655   0.0112 *  
#   date_of_dissection21/2/2019   5747556    2584284   2.224   0.0317 *  
#   EnrichmentEE                -11661733    2344651  -4.974 1.22e-05 ***
#   SocialDefeatSD               -4363542    2442516  -1.786   0.0814 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 7757000 on 41 degrees of freedom
# Multiple R-squared:  0.4365,	Adjusted R-squared:  0.3815 
# F-statistic: 7.939 on 4 and 41 DF,  p-value: 7.785e-05


table(Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
#         16/1/2019 17/1/2019 21/2/2019
# NIL        14         1         8
# SD          8         6         9
fisher.test(Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat and Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection
# p-value = 0.07168
# alternative hypothesis: two.sided

#That might explain it.

table(Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
#         16/1/2019 17/1/2019 21/2/2019
# NIL        12         4         6
# EE         10         3        11

# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment and Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection
# p-value = 0.4838
# alternative hypothesis: two.sided

fisher.test(Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)

table(Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
#       26.2.2019 28.2.2019
# NIL        15         8
# SD         12        11

fisher.test(Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat, Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
# Fisher's Exact Test for Count Data
# 
# data:  Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat and Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction
# p-value = 0.5499
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.4513058 6.6514582
# sample estimates:
# odds ratio 
# 1.698445 




#In review: Relationships that might be important:

#CountsPerSample~Enrichment (maybe also dissection)
#ratio_260_280~date_of_RNA_extraction
#RNA conc~Enrichment (trend)
#Dissection~Social Defeat (trend)

#Also: generation and processing batches look unevenly distributed by treatment group, but these relationships aren't significant.


###############################################

#PCA/MDS

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers")

#MDS plots

png("MDSplots_Group_DissectionBatch.png")
par(mfrow=c(2,2))
col.group <- Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group
#levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.batch <- as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
levels(col.batch) <-  brewer.pal(nlevels(col.batch), "Set2")
col.batch <- as.character(col.batch)
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)))
title(main="A. Sample groups")
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)), dim=c(3,4))
title(main="B. Sample groups")
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)), col=col.batch, dim=c(1,2))
title(main="C. Dissection Batch")
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)), col=col.batch, dim=c(3,4))
title(main="D. Dissection Batch")
dev.off()
#Looks like dissection batch 1 (Generation F49) may be related to PC2

png("MDSplots_Group_ExtractionBatch.png")
par(mfrow=c(2,2))
col.group <- Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group
#levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.batch <- as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
levels(col.batch) <-  brewer.pal(nlevels(col.batch), "Set2")
col.batch <- as.character(col.batch)
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)))
title(main="A. Sample groups")
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)), dim=c(3,4))
title(main="B. Sample groups")
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)), col=col.batch, dim=c(1,2))
title(main="C. Extraction Batch")
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)), col=col.batch, dim=c(3,4))
title(main="D. Extraction Batch")
dev.off()
#I'm having a hard time telling if Extraction and/or Dissection Batch are related to many of the top sources of variation in the data - maybe PC2?


png("MDSplots_Group_Generation.png")
par(mfrow=c(2,2))
col.group <- Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group
#levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
#col.group <- as.character(col.group)
col.batch <- Sample_MetaData_NACC_NoOutliers_Ordered$Generation
#levels(col.batch) <-  brewer.pal(nlevels(col.batch), "Set2")
#col.batch <- as.character(col.batch)
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)))
title(main="A. Sample groups")
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)), col=as.numeric(as.factor(col.group)), dim=c(3,4))
title(main="B. Sample groups")
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Generation)), col=as.numeric(as.factor(col.batch)), dim=c(1,2))
title(main="C. Generation")
plotMDS(NACC_RNASeq_Log2_Filtered, labels=as.numeric(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Generation)), col=as.numeric(as.factor(col.batch)), dim=c(3,4))
title(main="D. Generation")
dev.off()
#Generation (which is also dissection batch) looks like it matters for the 2nd & 3rd major gradients of variation in the data. We'll have to take a closer look.


#############################

#Switching over to principal components analysis (PCA)  to see if we can generate some concise co-variates that can be used to control for technical variation.

#Already had some sample-level normalization.
#Before running PCA, I centered and scaled by gene, so that highly variable genes aren't driving the effect.

str(NACC_RNASeq_Log2_Filtered)
#num [1:17765, 1:46] -3.15 -1.24 8.41 6.56 2.41 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17765] "ENSRNOG00000061316" "ENSRNOG00000029897" "ENSRNOG00000014303" "ENSRNOG00000014330" ...
# ..$ : chr [1:46] "Sample_125991" "Sample_125992" "Sample_125993" "Sample_125994" ...
NACC_RNASeq_Log2_Filtered_Scaled<-t(scale(t(NACC_RNASeq_Log2_Filtered), center=T, scale=T))

#Sanity Check:
head(apply(NACC_RNASeq_Log2_Filtered_Scaled, 1, mean))
# ENSRNOG00000061316 ENSRNOG00000029897 ENSRNOG00000014303 ENSRNOG00000014330 ENSRNOG00000049505 
# 9.614207e-17       5.170112e-17       3.054575e-16       2.805396e-15       4.124989e-16 
# ENSRNOG00000014916 
# 5.481726e-16 

head(apply(NACC_RNASeq_Log2_Filtered_Scaled, 1, sd))
# ENSRNOG00000061316 ENSRNOG00000029897 ENSRNOG00000014303 ENSRNOG00000014330 ENSRNOG00000049505 
# 1                  1                  1                  1                  1 
# ENSRNOG00000014916 
# 1 

head(apply(NACC_RNASeq_Log2_Filtered_Scaled, 2, mean))
# Sample_125991 Sample_125992 Sample_125993 Sample_125994 Sample_125995 Sample_125996 
#0.029551749   0.125653136   0.062034804   0.042376246  -0.001903309   0.016930185
#hmmm... I wonder if that amount of variability across samples could be problematic.


pcaNormFiltered<-prcomp(t(NACC_RNASeq_Log2_Filtered_Scaled))
tmp<-pcaNormFiltered$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")

PC1<-pcaNormFiltered$x[,1]
PC2<-pcaNormFiltered$x[,2]

PC3<-pcaNormFiltered$x[,3]
PC4<-pcaNormFiltered$x[,4]

tmp<-data.frame(gene_id=row.names(pcaNormFiltered$rotation),pcaNormFiltered$rotation[,c(1:4)])

NACC_Eigenvectors_Annotated<-join(tmp, NACC_RNASeq_Log2_Annotated[, c(1, 48:50)], by="gene_id", type="left")

write.csv(NACC_Eigenvectors_Annotated, "Eigenvectors_PCA_1_4.csv")

#Output a scree plot for the PCA:
png("PCA Scree Plot1.png")
plot(summary(pcaNormFiltered)$importance[2,]~(c(1:46)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 explains around 15% of the variation in the dataset, all the others are less than 15%. Pretty normal.

png("PCA Scree Plot2.png")
plot(summary(pcaNormFiltered)$importance[3,]~(c(1:46)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

#Quickly coloring plots using group as factor:
levels(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("NACC_PC1vsPC2_vsGroup.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
legend(max(PC2-50), min(PC1)+50, c("LR EE", "LR EE + SD", "LR Nil", "LR Nil + SD"), text.col=c(1, 2, 3, 4), pch=19, col=c(1, 2, 3, 4))
dev.off()
#PC1 isn't driven by outliers or group status
#PC2 looks related to group status.
#Looks like Fan's graph - good. 

png("NACC_PC1_vsGroup.png")
boxplot(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)
dev.off()
#PC1 actually does look related to group. Cool. - in particular, social defeat.

anova(lm(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))

# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value Pr(>F)
# Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group  3  12112  4037.4  1.4762 0.2348
# Residuals                                               42 114865  2734.9 
# 
#But not sig.

png("NACC_PC1_vsSD.png")
boxplot(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat)
dev.off()

anova(lm(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat))
# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat  1   8028  8028.4  2.9698 0.09186 .
# Residuals                                            44 118949  2703.4                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Close!

temp<-data.frame(PC1, Sample_MetaData_NACC_NoOutliers_Ordered)
anova(lm(PC1~Enrichment*SocialDefeat, data=temp))

# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Enrichment               1   1866  1865.6  0.6821 0.41352  
# SocialDefeat             1   8783  8782.8  3.2114 0.08033 .
# Enrichment:SocialDefeat  1   1464  1463.7  0.5352 0.46849  
# Residuals               42 114865  2734.9                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


png("NACC_PC2_vsGroup.png")
boxplot(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)
dev.off()
#That's cool - looks like social defeat is related to PC2- in particular, social defeat.

anova(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))

# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value Pr(>F)
# Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group  3   9583  3194.2  1.4735 0.2355
# Residuals                                               42  91048  2167.8  
# anova(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment))

# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value Pr(>F)
# Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment  1   6067  6066.7  2.5139   0.12
# Residuals                                          44 106185  2413.3 

png("NACC_PC2_vsSD.png")
boxplot(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat)
dev.off()

anova(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat))

# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value Pr(>F)
# Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat  1   2451  2451.4  1.0986 0.3003
# Residuals                                            44  98179  2231.3  

temp<-data.frame(PC2, Sample_MetaData_NACC_NoOutliers_Ordered)
anova(lm(PC2~Enrichment*SocialDefeat, data=temp))

# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value Pr(>F)
# Enrichment               1   4701  4700.7  2.1684 0.1483
# SocialDefeat             1   3101  3101.5  1.4307 0.2384
# Enrichment:SocialDefeat  1   1780  1780.4  0.8213 0.3700
# Residuals               42  91048  2167.8 

png("NACC_PC3_vsGroup.png")
boxplot(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)
dev.off()

anova(lm(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group))
#PC3 is an outlier


png("NACC_PC4_vsGroup.png")
boxplot(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group)
dev.off()
#PC4 also seems to be an outlier

#Is it related to total counts?
png("NACC_PC1vsCountsPerSample.png")
plot(PC1~NACC_RNASeq_CountsPerSample)
dev.off()
#Maybe?

summary.lm(lm(PC1~NACC_RNASeq_CountsPerSample))
# Call:
#   lm(formula = PC1 ~ NACC_RNASeq_CountsPerSample)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -131.75  -17.91    7.19   25.50  132.66 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                  3.650e+01  2.362e+01   1.546    0.129
# NACC_RNASeq_CountsPerSample -1.289e-06  7.883e-07  -1.635    0.109
# 
# Residual standard error: 52.16 on 44 degrees of freedom
# Multiple R-squared:  0.05726,	Adjusted R-squared:  0.03583 
# F-statistic: 2.672 on 1 and 44 DF,  p-value: 0.1092

anova(lm(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment*Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat+NACC_RNASeq_CountsPerSample))
#Counts per sample doesn't seem to mediate the relationship with group 

png("NACC_PC2vsCountsPerSample.png")
plot(PC2~NACC_RNASeq_CountsPerSample)
dev.off()
summary.lm(lm(PC2~NACC_RNASeq_CountsPerSample))

# Call:
#   lm(formula = PC2 ~ NACC_RNASeq_CountsPerSample)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -110.545  -29.801    3.851   37.145   73.930 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                  2.915e+01  2.115e+01   1.378    0.175
# NACC_RNASeq_CountsPerSample -1.029e-06  7.059e-07  -1.458    0.152
# 
# Residual standard error: 46.71 on 44 degrees of freedom
# Multiple R-squared:  0.04608,	Adjusted R-squared:  0.0244 
# F-statistic: 2.125 on 1 and 44 DF,  p-value: 0.152

anova(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$Enrichment*Sample_MetaData_NACC_NoOutliers_Ordered$SocialDefeat+NACC_RNASeq_CountsPerSample))


levels(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Generation))

png("NACC_PC1vsPC2_vsGeneration.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$Generation))
legend(max(PC2-50), min(PC1)+50, c("F49","F53", "F56"), text.col=c(1, 2, 3), pch=19, col=c(1, 2, 3))
dev.off()
#PC2 may be related to generation


png("NACC_PC1_vsGeneration.png")
boxplot(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$Generation)
dev.off()


png("NACC_PC2_vsGeneration.png")
boxplot(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$Generation)
dev.off()
#It is the F49 samples that look different
anova(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$Generation))

# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value  Pr(>F)   
# Sample_MetaData_NACC_NoOutliers_Ordered$Generation  2  25056 12527.8   7.128 0.00212 **
#   Residuals                                          43  75575  1757.6                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


png("NACC_PC3_vsGeneration.png")
boxplot(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$Generation)
dev.off()
#PC3 also looks related to generation - again, it's the F49 samples that look really different.
anova(lm(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$Generation))

# Analysis of Variance Table
# 
# Response: PC3
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Sample_MetaData_NACC_NoOutliers_Ordered$Generation  2  11051  5525.6  3.5633 0.03699 *
#   Residuals                                          43  66679  1550.7                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

png("NACC_PC4_vsGeneration.png")
boxplot(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$Generation)
dev.off()

anova(lm(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$Generation))


levels(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_sac))
#[1] "10-Jan-16" "2-Nov-17"  "24-Oct-17" "28-Oct-17" "6-Feb-17"  "8-Jan-16"  "9-Feb-17"  "9-Jan-16" 

png("NACC_PC1vsPC2_vsDateOfSacrifice.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_sac))
legend(max(PC2-50), min(PC1)+50, c("10-Jan-16", "2-Nov-17",  "24-Oct-17", "28-Oct-17", "6-Feb-17",  "8-Jan-16",  "9-Feb-17",  "9-Jan-16"), text.col=c(1, 2, 3, 4, 5,6,7,8), pch=19, col=c(1, 2, 3, 4, 5,6,7,8))
dev.off()
#PC2 seems related to date of sacrifice too.
#But date of sacrifice is basically each of the treatment subgroups in a generation.

levels(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction))
#[1] "26.2.2019" "28.2.2019"

png("NACC_PC1vsPC2_vsDateOfExtraction.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction))
legend(max(PC2-50), min(PC1)+50, c("26.2.2019", "28.2.2019"), text.col=c(1, 2), pch=19, col=c(1, 2))
dev.off()
#Date of Extraction may be related to PC2.

png("NACC_PC1_vsDateOfExtraction.png")
boxplot(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
dev.off()
#Doesn't look too bad.

anova(lm(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction))

# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value Pr(>F)
# Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction  1     11   10.57  0.0037  0.952
# Residuals                                                      44 126967 2885.61

png("NACC_PC2_vsDateOfExtraction.png")
boxplot(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
dev.off()
#pretty striking
anova(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction))
# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction  1  18106 18106.3  9.6539 0.003305 **
#   Residuals                                                      44  82524  1875.6                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(PC2, Sample_MetaData_NACC_NoOutliers_Ordered, NACC_RNASeq_CountsPerSample)
anova(lm(PC2~Enrichment*SocialDefeat+Generation+date_of_RNA_extraction, data=Temp))
# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Enrichment               1   4701  4700.7  3.0645 0.087886 . 
# SocialDefeat             1   3101  3101.5  2.0219 0.162995   
# Generation               2  21308 10654.1  6.9455 0.002629 **
#   date_of_RNA_extraction   1  10067 10067.0  6.5628 0.014392 * 
#   Enrichment:SocialDefeat  1   1629  1628.8  1.0619 0.309142   
# Residuals               39  59824  1534.0                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#So PC2 is independently related to both Generation and Date of RNA Extraction


png("NACC_PC3_vsDateOfExtraction.png")
boxplot(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
dev.off()

anova(lm(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction))


png("NACC_PC4_vsDateOfExtraction.png")
boxplot(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction)
dev.off()

anova(lm(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_RNA_extraction))


# RNA quality seemed to vary with RNA extraction date.
#So is it related to RNAIntegrity?
png("NACC_PC1vsRNAIntegrity.png")
plot(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280)
dev.off()
summary.lm(lm(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280))
#Nope.

png("NACC_PC2vsRNAIntegrity.png")
plot(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280)
dev.off()
summary.lm(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280))

# Call:
#   lm(formula = PC2 ~ Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -115.417  -25.961    3.291   42.131   67.031 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                              493.0      298.8    1.65    0.106
# Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280   -239.8      145.4   -1.65    0.106
# 
# Residual standard error: 46.41 on 44 degrees of freedom
# Multiple R-squared:  0.05827,	Adjusted R-squared:  0.03687 
# F-statistic: 2.723 on 1 and 44 DF,  p-value: 0.1061

#Maybe? Not sig.
Temp<-data.frame(PC2, Sample_MetaData_NACC_NoOutliers_Ordered, NACC_RNASeq_CountsPerSample)
anova(lm(PC2~Enrichment*SocialDefeat+Generation+date_of_RNA_extraction+ratio_260_280, data=Temp))

# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Enrichment               1   4701  4700.7  2.9883 0.091986 . 
# SocialDefeat             1   3101  3101.5  1.9716 0.168394   
# Generation               2  21308 10654.1  6.7729 0.003049 **
#   date_of_RNA_extraction   1  10067 10067.0  6.3997 0.015681 * 
#   ratio_260_280            1     67    67.3  0.0428 0.837267   
# Enrichment:SocialDefeat  1   1610  1610.1  1.0235 0.318079   
# Residuals               38  59776  1573.0                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Nope - doesn't seem to mediate the effects of RNA extraction date.

png("NACC_PC3vsRNAIntegrity.png")
plot(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280)
dev.off()
summary.lm(lm(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280))
#Nope.

png("NACC_PC4vsRNAIntegrity.png")
plot(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280)
dev.off()
summary.lm(lm(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$ratio_260_280))
#Nope.

#Interesting. So Extraction date correlates with RNA Integrity, but that does not mediate the main effects of Extraction date on gene expression.


png("NACC_PC1vsRNAconc.png")
plot(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc)
dev.off()
#looks like a slight negative correlation
summary.lm(lm(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc))
# Call:
#   lm(formula = PC1 ~ Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -104.931  -26.809   -5.421   22.251  123.026 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                       46.2636    20.4751   2.260   0.0289 *
#   Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc  -0.6895     0.2843  -2.425   0.0195 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 50.45 on 44 degrees of freedom
# Multiple R-squared:  0.1179,	Adjusted R-squared:  0.09787 
# F-statistic: 5.882 on 1 and 44 DF,  p-value: 0.01947

Temp<-data.frame(PC1, Sample_MetaData_NACC_NoOutliers_Ordered, NACC_RNASeq_CountsPerSample)
anova(lm(PC1~Enrichment*SocialDefeat+RNA_conc, data=Temp))

# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Enrichment               1   1866  1865.6  0.7461 0.39273  
# SocialDefeat             1   8783  8782.8  3.5126 0.06804 .
# RNA_conc                 1  11740 11740.4  4.6954 0.03611 *
#   Enrichment:SocialDefeat  1   2073  2072.5  0.8289 0.36792  
# Residuals               41 102516  2500.4                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


anova(lm(PC1~Enrichment+SocialDefeat+RNA_conc, data=Temp))

# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Enrichment    1   1866  1865.6  0.7492 0.39166  
# SocialDefeat  1   8783  8782.8  3.5269 0.06733 .
# RNA_conc      1  11740 11740.4  4.7146 0.03561 *
#   Residuals    42 104589  2490.2                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

png("NACC_PC2vsRNAconc.png")
plot(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc)
dev.off()
#looks like a slight positive correlation
summary.lm(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc))
#Not sig.

png("NACC_PC3vsRNAconc.png")
plot(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc)
dev.off()
summary.lm(lm(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc))

png("NACC_PC4vsRNAconc.png")
plot(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc)
dev.off()
summary.lm(lm(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$RNA_conc))



levels(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection))
#[1] "16/1/2019" "17/1/2019" "21/2/2019"



levels(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch) )
# [1] "16/1/2019 26.2.2019" "16/1/2019 28.2.2019" "17/1/2019 26.2.2019" "17/1/2019 28.2.2019"
# [5] "21/2/2019 26.2.2019" "21/2/2019 28.2.2019"

png("NACC_PC1vsPC2_vsDateOfExtractionDissection.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch))
#I erased the legend because it was going to be too large.
dev.off()

png("NACC_PC1_vsExtractionDissection.png")
boxplot(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch)
dev.off()
#Looks related to dissection batch
anova(lm(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch))
#Not sig.

png("NACC_PC1_vsDissection.png")
boxplot(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
dev.off()
#Yep.
anova(lm(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection))
# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection  2  13122  6561.2   2.478 0.09582 .
# Residuals                                                  43 113855  2647.8                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(PC1, Sample_MetaData_NACC_NoOutliers_Ordered)
summary.lm(lm(PC1~date_of_dissection, data=Temp))

# Call:
#   lm(formula = PC1 ~ date_of_dissection, data = Temp)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -128.369  -19.344    1.245   22.969  111.337 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                  -6.9787    10.9706  -0.636   0.5281  
# date_of_dissection17/1/2019  46.8431    22.3296   2.098   0.0418 *
#   date_of_dissection21/2/2019  -0.4048    16.6164  -0.024   0.9807  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 51.46 on 43 degrees of freedom
# Multiple R-squared:  0.1033,	Adjusted R-squared:  0.06164 
# F-statistic: 2.478 on 2 and 43 DF,  p-value: 0.09582

#Interesting - I would have expected 1/16 and 1/17 to be similar because they are almost the same day. Maybe it is due to the imbalance across treatment groups?

Temp<-data.frame(PC1, Sample_MetaData_NACC_NoOutliers_Ordered, NACC_RNASeq_CountsPerSample)
anova(lm(PC1~Enrichment*SocialDefeat+RNA_conc+date_of_dissection, data=Temp))

# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Enrichment               1   1866  1865.6  0.7887 0.37995  
# SocialDefeat             1   8783  8782.8  3.7129 0.06130 .
# RNA_conc                 1  11740 11740.4  4.9632 0.03173 *
#   date_of_dissection       2   9052  4526.1  1.9134 0.16118  
# Enrichment:SocialDefeat  1   3283  3282.8  1.3878 0.24592  
# Residuals               39  92254  2365.5                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Doesn't mediate the effect of RNAConc

# Call:
#   lm(formula = PC1 ~ Enrichment * SocialDefeat + RNA_conc + date_of_dissection, 
#      data = Temp)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -83.978 -27.473   1.457  29.179 127.669 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                  36.7533    27.6902   1.327   0.1921  
# EnrichmentEE                 -8.7267    21.6136  -0.404   0.6886  
# SocialDefeatSD               -0.7615    22.6053  -0.034   0.9733  
# RNA_conc                     -0.6444     0.2872  -2.244   0.0306 *
#   date_of_dissection17/1/2019  39.1451    22.9455   1.706   0.0960 .
# date_of_dissection21/2/2019  -7.7554    16.4283  -0.472   0.6395  
# EnrichmentEE:SocialDefeatSD  34.8300    29.5660   1.178   0.2459  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 48.64 on 39 degrees of freedom
# Multiple R-squared:  0.2735,	Adjusted R-squared:  0.1617 
# F-statistic: 2.447 on 6 and 39 DF,  p-value: 0.04193

#17-1 is definitely the weird dissection day.

# In a model that also included enrichment, social defeat, and RNA concentration, both the relationship between PC1 and dissection batch and PC1 and social defeat weakens, probably because the two variables are partially collinear.

Temp<-data.frame(PC1, Sample_MetaData_NACC_NoOutliers_Ordered, NACC_RNASeq_CountsPerSample)
summary.lm(lm(PC1~Enrichment*SocialDefeat+RNA_conc+date_of_dissection, data=Temp))

Temp<-data.frame(PC1, Sample_MetaData_NACC_NoOutliers_Ordered, NACC_RNASeq_CountsPerSample)
anova(lm(PC1~Enrichment*SocialDefeat+RNA_conc+CountsPerSample+date_of_dissection, data=Temp))
# Analysis of Variance Table
# 
# Response: PC1
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Enrichment               1   1866  1865.6  0.8207 0.37067  
# SocialDefeat             1   8783  8782.8  3.8639 0.05667 .
# RNA_conc                 1  11740 11740.4  5.1651 0.02879 *
#   CountsPerSample          1   3827  3827.4  1.6838 0.20224  
# date_of_dissection       2  12555  6277.5  2.7617 0.07588 .
# Enrichment:SocialDefeat  1   1831  1830.9  0.8055 0.37511  
# Residuals               38  86375  2273.0                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



png("NACC_PC2_vsExtractionDissection.png")
boxplot(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch)
dev.off()


png("NACC_PC2_vsDissection.png")
boxplot(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
dev.off()

anova(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection))
# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value  Pr(>F)   
# Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection  2  25191 12595.6  7.1794 0.00204 **
#   Residuals                                                  43  75439  1754.4                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Yep dissection matters there too.

Temp<-data.frame(PC2,Sample_MetaData_NACC_NoOutliers_Ordered)
summary.lm(lm(PC2~date_of_dissection, data=Temp))

# Call:
#   lm(formula = PC2 ~ date_of_dissection, data = Temp)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -95.372 -21.933   4.111  26.535  78.201 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   -21.14       8.93  -2.367 0.022490 *  
#   date_of_dissection17/1/2019    15.18      18.18   0.835 0.408413    
# date_of_dissection21/2/2019    50.95      13.53   3.767 0.000497 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 41.89 on 43 degrees of freedom
# Multiple R-squared:  0.2503,	Adjusted R-squared:  0.2155 
# F-statistic: 7.179 on 2 and 43 DF,  p-value: 0.00204

#And in this case, 2/21 appears to be the weird batch.

#But we already knew that generation mattered, and generation and dissecton are closely related. Which matters more?

Temp<-data.frame(PC2, Sample_MetaData_NACC_NoOutliers_Ordered, NACC_RNASeq_CountsPerSample)
anova(lm(PC2~Enrichment*SocialDefeat+date_of_dissection+date_of_RNA_extraction, data=Temp))

# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Enrichment               1   4701  4700.7  3.0411 0.089058 . 
# SocialDefeat             1   3101  3101.5  2.0065 0.164565   
# date_of_dissection       2  20361 10180.5  6.5863 0.003433 **
#   date_of_RNA_extraction   1  10717 10716.9  6.9333 0.012065 * 
#   Enrichment:SocialDefeat  1   1467  1467.3  0.9493 0.335909   
# Residuals               39  60283  1545.7                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(lm(PC2~Enrichment*SocialDefeat+Generation+date_of_RNA_extraction, data=Temp))
# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Enrichment               1   4701  4700.7  3.0645 0.087886 . 
# SocialDefeat             1   3101  3101.5  2.0219 0.162995   
# Generation               2  21308 10654.1  6.9455 0.002629 **
#   date_of_RNA_extraction   1  10067 10067.0  6.5628 0.014392 * 
#   Enrichment:SocialDefeat  1   1629  1628.8  1.0619 0.309142   
# Residuals               39  59824  1534.0                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#The two models seem relatively equivalent. Generation does a slightly better job of explaining the data (lower residual sum of squares)


png("NACC_PC3_vsExtractionDissection.png")
boxplot(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch)
dev.off()
#Maybe a relationship with dissection batch?
anova(lm(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch))
# Analysis of Variance Table
# 
# Response: PC3
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch  5  15689  3137.8   2.023 0.09608 .
# Residuals                                                         40  62042  1551.0                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


png("NACC_PC3_vsDissection.png")
boxplot(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection)
dev.off()

anova(lm(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection))
# Analysis of Variance Table
# 
# Response: PC3
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection  2  12332  6165.9  4.0541 0.02438 *
#   Residuals                                                  43  65399  1520.9                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Yep - dissection.  
Temp<-data.frame(PC3,Sample_MetaData_NACC_NoOutliers_Ordered)
summary.lm(lm(PC3~date_of_dissection, data=Temp))
# Call:
#   lm(formula = PC3 ~ date_of_dissection, data = Temp)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -200.268  -17.237    3.341   20.244   49.918 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                  -15.904      8.315  -1.913  0.06245 . 
# date_of_dissection17/1/2019   17.499     16.923   1.034  0.30691   
# date_of_dissection21/2/2019   35.829     12.594   2.845  0.00677 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 39 on 43 degrees of freedom
# Multiple R-squared:  0.1586,	Adjusted R-squared:  0.1195 
# F-statistic: 4.054 on 2 and 43 DF,  p-value: 0.02438

#2/21 is the weird batch for this one too.


#But we already knew there was a relationship with generation, which is highly related. Which is stronger?
anova(lm(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$Generation))
# Analysis of Variance Table
# 
# Response: PC3
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Sample_MetaData_NACC_NoOutliers_Ordered$Generation  2  11051  5525.6  3.5633 0.03699 *
#   Residuals                                          43  66679  1550.7                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#The relationship with dissection is a better descriptor of the data.

png("NACC_PC4_vsExtractionDissection.png")
boxplot(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$DissectionExtractionBatch)
dev.off()


levels(as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected))
#[1] "left"  "right"

png("NACC_PC1vsPC2_vsHemisphere.png")
plot(PC1~PC2, main="Principal Components Analysis of Filtered Data", col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected))
legend(max(PC2-50), min(PC1)+80, c("left", "right"), text.col=c(1, 2), pch=19, col=c(1, 2))
dev.off()
#Maybe a relationship with PC2?


png("NACC_PC1_vsHemisphere.png")
boxplot(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected)
dev.off()
summary.lm(lm(PC1~Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected))

png("NACC_PC2_vsHemisphere.png")
boxplot(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected)
dev.off()
summary.lm(lm(PC2~Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected))
#Nope, very much not sig.

png("NACC_PC3_vsHemisphere.png")
boxplot(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected)
dev.off()
summary.lm(lm(PC3~Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected))

png("NACC_PC4_vsHemisphere.png")
boxplot(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected)
dev.off()
summary.lm(lm(PC4~Sample_MetaData_NACC_NoOutliers_Ordered$hemisphere_dissected))

#Hemisphere definitely doesn't seem to matter (unsurprisingly)

colnames(Sample_MetaData_NACC_NoOutliers_Ordered)
colnames(Sample_BehavHormonalData_NACC_NoOutliers_Ordered)

write.csv(cor(cbind(PC1, PC2, PC3, PC4, as.matrix(Sample_MetaData_NACC_NoOutliers_Ordered[, c(13:14, 16)]), as.matrix(Sample_BehavHormonalData_NACC_NoOutliers_Ordered[,c(10:21,35:46, 48:51)])), use="pairwise.complete.obs"), "PCA_vsNumericVar_CorMatrix.csv")
#Note:  many of these variables are really really not distributed normally - let's try running a spearman instead.

write.csv(cor(cbind(PC1, PC2, PC3, PC4, as.matrix(Sample_MetaData_NACC_NoOutliers_Ordered[, c(13:14, 16)]), as.matrix(Sample_BehavHormonalData_NACC_NoOutliers_Ordered[,c(10:21,35:46, 48:51)])), use="pairwise.complete.obs", method="spearman"), "PCA_vsNumericVar_CorMatrix_Spearman.csv")

#Interesting. Sadly, many variables that seem to have the strongest relationships are the one's that probably also have puny sample sizes - i.e., probably just noise.

plot(PC1~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$cranky_USVs)
#Huh. That's actually kind-of pretty.
#Too bad that is one of the variables I need to recalculate due to the big generational differences in recording. :(

plot(PC1~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$testosterone)
#The sample size for that actually isn't completely pathetic, but the data is noisy.
#Maybe worth following up on.

plot(PC3~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$AGG_DefeatDay3)
#That is just a super-dinky sample size driving that effect (only 3 aggressive animals)


plot(PC3~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$SUB_DefeatDay1)
#That might actually be something, although it is just the defeated animals.
plot(PC3~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$SUB_DefeatDay1, ylim=c(-100, 100))
#whereas that might potentially be interesting - although a little noisy.


plot(PC3~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$SUB_DefeatDay1, col=as.factor(Sample_BehavHormonalData_NACC_NoOutliers_Ordered$Enrichment))
#...and it isn't driven by enrichment status within the group.
boxplot(PC3~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$Social_Defeat)
#But it would be weird for the relationship to go in the opposite of the general relationship with social defeat, so probably just an artifact.

plot(PC3~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_on_top_of_stimulus_animal_Ethovision)
#A few animals are driving the negative correlation. Again, pretty uninspiring.


plot(PC4~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT)
#Also pretty uninspiring.


plot(PC2~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$IL.6..pg.ml.)
#That may be a negative correlation, but a little noisy.

plot(PC3~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$Oxytocin..pg.ml.)
#Relationship seems driven by an outlier - outlier for both oxytocin and PC3.
plot(PC4~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$Oxytocin..pg.ml.)
#That sample is also an outlier for PC4. Maybe we should just cut them out?


#I should probably take a look at the eigenvectors and see if there are any hints. :(
#Nope, nothing obvious.



#So it looks like the most important independent co-variates are:
# PC1
# social defeat -trend
# RNA-conc*
# date of dissection -trend
# 
# PC2
# enrichment - trend
# generation or date of dissection*
# RNA extraction*
# 
# PC3
# generation or date of dissection*
# 

#So maybe two models:
# y~SD*EE+RNAconc+RNAextract+DissectionDate
# y~SD*EE+RNAconc+RNAextract+Generation

#We should maybe also try versions with fewer of those co-variates, since only CountsPerSample was collinear with treatment group
#Although there was a trend towards a relationship between dissection and social defeat and a trend towards a relationship between RNA concentration and enrichment.

#Speaking of which, I wonder if we should include CountsPerSample for that reason. :(
#No - going back to the results, it doesn't seem strongly related to the PCs - I think the TMM normalization helped with that.



##############################

#Limma-voom analysis:

#

design <- model.matrix(~Enrichment*SocialDefeat+RNA_conc+date_of_RNA_extraction+Generation, data=Sample_MetaData_NACC_NoOutliers_Ordered)

design
# (Intercept) EnrichmentEE SocialDefeatSD RNA_conc date_of_RNA_extraction28.2.2019 GenerationF49
# 68           1            0              0     47.0                               0             0
# 69           1            0              0     67.4                               0             0
# 70           1            0              0    114.9                               1             0
# 61           1            0              0     61.3                               1             0
# 62           1            0              0     50.3                               1             0
# 63           1            0              0     92.2                               0             0
# 64           1            0              0    156.4                               0             0
# 65           1            0              0     58.0                               1             0
# 41           1            0              1     91.2                               0             0
# 42           1            0              1     40.6                               1             0
# 43           1            0              1     98.9                               0             0
# 44           1            0              1     68.8                               0             0
# 45           1            0              1     81.0                               0             0
# 46           1            0              1     97.8                               0             0
# 47           1            0              1     97.8                               0             0
# 48           1            0              1     54.2                               0             0
# 74           1            1              0     63.3                               0             0
# 75           1            1              0     50.7                               0             0
# 76           1            1              0     47.8                               0             0
# 32           1            1              0     42.7                               1             0
# 33           1            1              0     51.5                               0             0
# 34           1            1              0     56.6                               0             0
# 35           1            1              0     78.1                               0             0
# 71           1            1              1     76.6                               1             0
# 72           1            1              1     69.5                               1             0
# 73           1            1              1     55.7                               0             0
# 37           1            1              1     62.7                               1             0
# 38           1            1              1     27.7                               1             0
# 40           1            1              1     55.4                               0             0
# 66           1            0              0     81.8                               0             1
# 67           1            0              0     54.8                               0             1
# 77           1            0              1     90.2                               1             1
# 78           1            0              1     35.9                               1             1
# 79           1            0              1     36.5                               1             1
# 80           1            0              1     61.4                               0             1
# 49           1            1              0     36.1                               0             1
# 50           1            1              0     24.0                               1             1
# 51           1            1              0     89.2                               1             1
# 52           1            1              0    111.8                               1             1
# 53           1            1              0     68.3                               0             1
# 54           1            1              0     67.5                               0             1
# 55           1            1              1     93.4                               1             1
# 56           1            1              1     19.7                               0             1
# 57           1            1              1     79.1                               1             1
# 58           1            1              1     48.6                               1             1
# 60           1            1              1     72.1                               0             1
# GenerationF53 EnrichmentEE:SocialDefeatSD
# 68             1                           0
# 69             1                           0
# 70             1                           0
# 61             0                           0
# 62             0                           0
# 63             0                           0
# 64             0                           0
# 65             0                           0
# 41             0                           0
# 42             0                           0
# 43             0                           0
# 44             0                           0
# 45             0                           0
# 46             0                           0
# 47             0                           0
# 48             0                           0
# 74             1                           0
# 75             1                           0
# 76             1                           0
# 32             0                           0
# 33             0                           0
# 34             0                           0
# 35             0                           0
# 71             1                           1
# 72             1                           1
# 73             1                           1
# 37             0                           1
# 38             0                           1
# 40             0                           1
# 66             0                           0
# 67             0                           0
# 77             0                           0
# 78             0                           0
# 79             0                           0
# 80             0                           0
# 49             0                           0
# 50             0                           0
# 51             0                           0
# 52             0                           0
# 53             0                           0
# 54             0                           0
# 55             0                           1
# 56             0                           1
# 57             0                           1
# 58             0                           1
# 60             0                           1
# 
# attr(,"assign")
# [1] 0 1 2 3 4 5 5 6
# attr(,"contrasts")
# attr(,"contrasts")$Enrichment
# [1] "contr.treatment"
# 
# attr(,"contrasts")$SocialDefeat
# [1] "contr.treatment"
# 
# attr(,"contrasts")$date_of_RNA_extraction
# [1] "contr.treatment"
# 
# attr(,"contrasts")$Generation
# [1] "contr.treatment"


#Man that seems like a lot of df for a dinky sample size (46 samples, 8 df in model)...

v <- voom(NACC_RNASeq_dge_TMM, design, plot=TRUE)
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

#     (Intercept) EnrichmentEE SocialDefeatSD RNA_conc date_of_RNA_extraction28.2.2019 GenerationF49
# -1        3673            5              1        1                               0           365
# 0         1226        17758          17761    17755                           17765         16827
# 1        12866            2              3        9                               0           573
#       GenerationF53 EnrichmentEE:SocialDefeatSD
# -1             3                           0
# 0          17749                       17765
# 1             13                           0

#1 means upregulated, -1 means downregulated, threshold = FDR<0.05. Therefore, using this analysis there are 5 genes down-regulated and 2 upregulated with enrichment, 1 gene downregulated and 3 upregulated with social defeat, many genes that differ by generation (vs. F49 in particular), and only one gene that differs with RNAconc.

toptable(efit, coef=2)
#Annotation is in ENSEMBL - will need to look at this offline in excel with additional annotation to be able to interpret it.

write.fit(efit, adjust="BH", file="Limma_results_Model_EESD_RNAconc_RNAextract_Gen.txt")

write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")

##################

#Using dissection instead of generation:

design <- model.matrix(~Enrichment*SocialDefeat+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=Sample_MetaData_NACC_NoOutliers_Ordered)

design

# (Intercept) EnrichmentEE SocialDefeatSD RNA_conc date_of_RNA_extraction28.2.2019
# 68           1            0              0     47.0                               0
# 69           1            0              0     67.4                               0
# 70           1            0              0    114.9                               1
# 61           1            0              0     61.3                               1
# 62           1            0              0     50.3                               1
# 63           1            0              0     92.2                               0
# 64           1            0              0    156.4                               0
# 65           1            0              0     58.0                               1
# 41           1            0              1     91.2                               0
# 42           1            0              1     40.6                               1
# 43           1            0              1     98.9                               0
# 44           1            0              1     68.8                               0
# 45           1            0              1     81.0                               0
# 46           1            0              1     97.8                               0
# 47           1            0              1     97.8                               0
# 48           1            0              1     54.2                               0
# 74           1            1              0     63.3                               0
# 75           1            1              0     50.7                               0
# 76           1            1              0     47.8                               0
# 32           1            1              0     42.7                               1
# 33           1            1              0     51.5                               0
# 34           1            1              0     56.6                               0
# 35           1            1              0     78.1                               0
# 71           1            1              1     76.6                               1
# 72           1            1              1     69.5                               1
# 73           1            1              1     55.7                               0
# 37           1            1              1     62.7                               1
# 38           1            1              1     27.7                               1
# 40           1            1              1     55.4                               0
# 66           1            0              0     81.8                               0
# 67           1            0              0     54.8                               0
# 77           1            0              1     90.2                               1
# 78           1            0              1     35.9                               1
# 79           1            0              1     36.5                               1
# 80           1            0              1     61.4                               0
# 49           1            1              0     36.1                               0
# 50           1            1              0     24.0                               1
# 51           1            1              0     89.2                               1
# 52           1            1              0    111.8                               1
# 53           1            1              0     68.3                               0
# 54           1            1              0     67.5                               0
# 55           1            1              1     93.4                               1
# 56           1            1              1     19.7                               0
# 57           1            1              1     79.1                               1
# 58           1            1              1     48.6                               1
# 60           1            1              1     72.1                               0
# date_of_dissection17/1/2019 date_of_dissection21/2/2019 EnrichmentEE:SocialDefeatSD
# 68                           0                           0                           0
# 69                           0                           0                           0
# 70                           0                           0                           0
# 61                           0                           0                           0
# 62                           0                           0                           0
# 63                           0                           0                           0
# 64                           0                           0                           0
# 65                           0                           0                           0
# 41                           0                           0                           0
# 42                           1                           0                           0
# 43                           0                           0                           0
# 44                           0                           0                           0
# 45                           1                           0                           0
# 46                           1                           0                           0
# 47                           1                           0                           0
# 48                           0                           0                           0
# 74                           0                           0                           0
# 75                           0                           0                           0
# 76                           0                           0                           0
# 32                           0                           0                           0
# 33                           1                           0                           0
# 34                           0                           0                           0
# 35                           0                           0                           0
# 71                           0                           0                           1
# 72                           0                           0                           1
# 73                           0                           0                           1
# 37                           0                           0                           1
# 38                           1                           0                           1
# 40                           1                           0                           1
# 66                           0                           1                           0
# 67                           0                           1                           0
# 77                           0                           1                           0
# 78                           0                           1                           0
# 79                           0                           1                           0
# 80                           0                           1                           0
# 49                           0                           1                           0
# 50                           0                           1                           0
# 51                           0                           1                           0
# 52                           0                           1                           0
# 53                           0                           1                           0
# 54                           0                           1                           0
# 55                           0                           1                           1
# 56                           0                           1                           1
# 57                           0                           1                           1
# 58                           0                           1                           1
# 60                           0                           1                           1
# attr(,"assign")
# [1] 0 1 2 3 4 5 5 6
# attr(,"contrasts")
# attr(,"contrasts")$Enrichment
# [1] "contr.treatment"
# 
# attr(,"contrasts")$SocialDefeat
# [1] "contr.treatment"
# 
# attr(,"contrasts")$date_of_RNA_extraction
# [1] "contr.treatment"
# 
# attr(,"contrasts")$date_of_dissection
# [1] "contr.treatment"

v <- voom(NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_Cov_RNAConcExtractDissection_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

#       (Intercept) EnrichmentEE SocialDefeatSD RNA_conc date_of_RNA_extraction28.2.2019
# -1        3727            4              2        1                               0
# 0         1157        17759          17761    17759                           17765
# 1        12881            2              2        5                               0
# date_of_dissection17/1/2019 date_of_dissection21/2/2019 EnrichmentEE:SocialDefeatSD
# -1                           0                        1010                           0
# 0                        17765                       15691                       17765
# 1                            0                        1064                           0

#Wow - one of those dissection dates really really matters.
#Date of RNA extraction doesn't seem to be strongly related to anything.

topTable(efit, coef=2)

#                         logFC  AveExpr         t      P.Value   adj.P.Val         B
# ENSRNOG00000020084 -0.4747789 3.318758 -6.414356 9.526874e-08 0.001692449 7.7122940
# ENSRNOG00000034126 -0.7557070 4.263064 -6.203590 1.923939e-07 0.001708939 7.1008828
# ENSRNOG00000020073  0.7280289 3.326349  5.775905 7.999075e-07 0.003552589 5.7246974
# ENSRNOG00000054953 -1.4335117 1.290944 -5.962860 4.292734e-07 0.002542014 5.7106858
# ENSRNOG00000019276 -0.3418350 5.723055 -5.623126 1.328880e-06 0.004721511 5.2275926
# ENSRNOG00000046264  0.4036324 3.020151  5.040181 9.078238e-06 0.026879149 3.4719236

topTable(efit, coef=3)

#                     logFC    AveExpr         t      P.Value   adj.P.Val        B
# ENSRNOG00000015248  1.6120580  1.7976337  5.797414 7.446765e-07 0.008481012 5.589978
# ENSRNOG00000029386  0.9975238  3.4069905  5.710969 9.926471e-07 0.008481012 5.545134
# ENSRNOG00000012546 -0.4169003  3.4445073 -5.002844 1.025527e-05 0.045546235 3.343730

write.fit(efit, adjust="BH", file="Limma_results_Model_EESD_RNAconc_RNAextract_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")

##############

#Trying out reduced models to determine the sensitivity of our results to model specification?


#What if we cut out RNA-extraction since it seems like a weaker co-variate? (only weakly related to PC2, no genes with FDR<0.05):

design <- model.matrix(~Enrichment*SocialDefeat+RNA_conc+date_of_dissection, data=Sample_MetaData_NACC_NoOutliers_Ordered)

design

v <- voom(NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_Cov_RNAConcDissection_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD RNA_conc date_of_dissection17/1/2019 date_of_dissection21/2/2019
# -1        3787            4              2        0                           0                        1316
# 0         1116        17759          17761    17764                       17765                       15166
# 1        12862            2              2        1                           0                        1283
# EnrichmentEE:SocialDefeatSD
# -1                           0
# 0                        17765
# 1                            0

#Strangely, that mostly affected the relationships with RNA-Conc and date of dissection.

topTable(efit, coef=2)
# logFC  AveExpr         t      P.Value   adj.P.Val        B
# ENSRNOG00000020084 -0.4747104 3.318758 -6.521548 6.065862e-08 0.001077600 8.097658
# ENSRNOG00000034126 -0.7755385 4.263064 -6.295763 1.298039e-07 0.001152983 7.462840
# ENSRNOG00000020073  0.7289339 3.326349  5.881340 5.239936e-07 0.002327187 6.092169
# ENSRNOG00000054953 -1.4603387 1.290944 -6.003225 3.477351e-07 0.002059171 5.788799
# ENSRNOG00000019276 -0.3463881 5.723055 -5.766001 7.720043e-07 0.002742931 5.756421
# ENSRNOG00000046264  0.4061049 3.020151  5.158159 5.865706e-06 0.017367377 3.860436

topTable(efit, coef=3)

#                     logFC    AveExpr         t      P.Value  adj.P.Val         B
# ENSRNOG00000015248  1.6730685  1.7976337  5.975543 3.816889e-07 0.00597786 6.1714325
# ENSRNOG00000029386  1.0018227  3.4069905  5.806873 6.729930e-07 0.00597786 5.9043995
# ENSRNOG00000012546 -0.4201364  3.4445073 -5.085629 7.454186e-06 0.03310590 3.6396412

write.fit(efit, adjust="BH", file="Limma_results_Model_EESD_RNAconc_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")


#What if we cut out RNA-extraction since it seems like a weaker co-variate? (only weakly related to PC2, no genes with FDR<0.05):

design <- model.matrix(~Enrichment*SocialDefeat+RNA_conc+date_of_dissection, data=Sample_MetaData_NACC_NoOutliers_Ordered)

design

v <- voom(NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_Cov_RNAConcDissection_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD RNA_conc date_of_dissection17/1/2019 date_of_dissection21/2/2019
# -1        3787            4              2        0                           0                        1316
# 0         1116        17759          17761    17764                       17765                       15166
# 1        12862            2              2        1                           0                        1283
# EnrichmentEE:SocialDefeatSD
# -1                           0
# 0                        17765
# 1                            0

#Strangely, that mostly affected the relationships with RNA-Conc and date of dissection.

topTable(efit, coef=2)
# logFC  AveExpr         t      P.Value   adj.P.Val        B
# ENSRNOG00000020084 -0.4747104 3.318758 -6.521548 6.065862e-08 0.001077600 8.097658
# ENSRNOG00000034126 -0.7755385 4.263064 -6.295763 1.298039e-07 0.001152983 7.462840
# ENSRNOG00000020073  0.7289339 3.326349  5.881340 5.239936e-07 0.002327187 6.092169
# ENSRNOG00000054953 -1.4603387 1.290944 -6.003225 3.477351e-07 0.002059171 5.788799
# ENSRNOG00000019276 -0.3463881 5.723055 -5.766001 7.720043e-07 0.002742931 5.756421
# ENSRNOG00000046264  0.4061049 3.020151  5.158159 5.865706e-06 0.017367377 3.860436

topTable(efit, coef=3)

#                     logFC    AveExpr         t      P.Value  adj.P.Val         B
# ENSRNOG00000015248  1.6730685  1.7976337  5.975543 3.816889e-07 0.00597786 6.1714325
# ENSRNOG00000029386  1.0018227  3.4069905  5.806873 6.729930e-07 0.00597786 5.9043995
# ENSRNOG00000012546 -0.4201364  3.4445073 -5.085629 7.454186e-06 0.03310590 3.6396412

write.fit(efit, adjust="BH", file="Limma_results_Model_EESD_RNAconc_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")

#########################

#One more try, even simpler design:

design <- model.matrix(~Enrichment*SocialDefeat+date_of_dissection, data=Sample_MetaData_NACC_NoOutliers_Ordered)

design

v <- voom(NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_Cov_Dissection_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD date_of_dissection17/1/2019 date_of_dissection21/2/2019
# -1        4203            4              1                           0                        1201
# 0          576        17760          17762                       17765                       15402
# 1        12986            1              2                           0                        1162
# EnrichmentEE:SocialDefeatSD
# -1                           0
# 0                        17765
# 1                            0


topTable(efit, coef=2)

#                         logFC     AveExpr         t      P.Value    adj.P.Val         B
# ENSRNOG00000034126 -0.7750333  4.26306412 -6.481965 6.298403e-08 0.0005594556 8.1526113
# ENSRNOG00000020084 -0.4630315  3.31875818 -6.491066 6.106155e-08 0.0005594556 8.1379051
# ENSRNOG00000054953 -1.4207015  1.29094408 -5.902017 4.534607e-07 0.0026852433 5.7551052
# ENSRNOG00000019276 -0.3349146  5.72305500 -5.697189 9.083406e-07 0.0040341677 5.5671680
# ENSRNOG00000020073  0.6906459  3.32634940  5.629034 1.143950e-06 0.0040644556 5.3949862

topTable(efit, coef=3)

# logFC     AveExpr         t      P.Value   adj.P.Val        B
# ENSRNOG00000015248  1.6818241  1.79763367  6.062285 2.629754e-07 0.004671757 6.441264
# ENSRNOG00000029386  0.9974775  3.40699051  5.833531 5.721647e-07 0.005082253 6.044016
# ENSRNOG00000012546 -0.4178142  3.44450728 -5.088037 7.035137e-06 0.041659736 3.696667

write.fit(efit, adjust="BH", file="Limma_results_Model_EESD_Dissect_SampleCounts.txt")



#########################

#Simplest model:

design <- model.matrix(~Enrichment*SocialDefeat, data=Sample_MetaData_NACC_NoOutliers_Ordered)

design

v <- voom(NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("FullModel_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

#       (Intercept) EnrichmentEE SocialDefeatSD EnrichmentEE:SocialDefeatSD
# -1        4214            5              2                           0
# 0          574        17759          17759                       17765
# 1        12977            1              4                           0


topTable(efit, coef=2)

#                       logFC     AveExpr         t      P.Value    adj.P.Val        B
# ENSRNOG00000020084 -0.4582719  3.31875818 -6.778126 1.912312e-08 0.0003397223 9.240607
# ENSRNOG00000034126 -0.7573953  4.26306412 -6.371376 7.845572e-08 0.0006968829 7.932173
# ENSRNOG00000019276 -0.3358228  5.72305500 -5.952466 3.349942e-07 0.0014877932 6.503507
# ENSRNOG00000054953 -1.4065343  1.29094408 -5.994628 2.895397e-07 0.0014877932 6.235834
# ENSRNOG00000020073  0.6818755  3.32634940  5.839055 4.956636e-07 0.0017610926 6.176401
# ENSRNOG00000039455 -0.5045962  4.02550951 -5.164973 4.981060e-06 0.0147480871 3.979227

topTable(efit, coef=3)

#                         logFC   AveExpr         t      P.Value   adj.P.Val        B
# ENSRNOG00000029386  0.9798264  3.406991  6.253133 1.182440e-07 0.001050302 7.453252
# ENSRNOG00000015248  1.7061597  1.797634  6.348271 8.500495e-08 0.001050302 7.217284
# ENSRNOG00000039744  0.5921371  1.660975  5.269372 3.495170e-06 0.015522923 4.090003
# ENSRNOG00000008970 -0.4585400  7.135523 -5.011683 8.356934e-06 0.029692188 3.464244
# ENSRNOG00000045924  0.8202470  1.873462  4.890367 1.255409e-05 0.037170563 3.015210


write.fit(efit, adjust="BH", file="Limma_results_Model_EESD.txt")


#########################

#The top genes appear to be very similar across models. Excellent - that means that the results are robust.

#Here are the results from the simplest model:

#Most of the top (FDR<0.10) enrichment-related genes are protocadherins:
# Pcdhb5
# Pcdhb6
# Pcdhb8
# Pcdhga2
# Pcdhgb6
# Dele1
#With the exception of Pcdhb8, they are all down-regulated.

#The top (FDR<0.10) social defeat related genes I don't recognize: 
# RT1-N2
# AC120486.2
# Abca12
# RT1-CE4
# Pcdh17
# RT1-T24-3
#apparently the RT-1s are peptide-antigen binding and orthologous to HLA-E. interesting.

#And here are the results from the most complex model:
#Most of the top (FDR<0.10) enrichment-related genes are protocadherins:
# Pcdhb5
# Pcdhb6
# Pcdhb8
# Pcdhga2
# Dele1
# Pcdhb7
# AC120486.2

#The top (FDR<0.10) social defeat related genes I don't recognize: 
# AC120486.2
# RT1-N2
# Abca12
# Frmpd1
# LOC360919
# Polr3g

#The relationships with dissection are really strong (many p-values are just listed as 0)



###################

#The results look almost identical to the version without the technical co-variates. So they're robust. Good.

#I reran the most complete model above.
#y~EE*SD+RNAconc+RNAExtract+DateOfDissect

#And then extracted the parts that I need to do other things:

str(efit)#The parts of the object don't appear named. 

head(efit$coefficients)
head(efit$t)
head(efit$cov.coefficients)
head(efit$p.value)

NAcc_Coefficients_EEbySD_RNAConcExtractDissect<-efit$coefficients
#Annotation

NAcc_LimmaResults_EEbySD_RNAConcExtractDissect<-data.frame(NACC_RNASeq_Log2_Annotated$gene_id, NACC_RNASeq_Log2_Annotated$gene_symbol, efit$coefficients, efit$t, efit$p.value, stringsAsFactors = F)
str(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect)

write.csv(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, "NAcc_LimmaResults_EEbySD_RNAConcExtractDissect.csv")

########################


#Looking at the relationship between variables:

str(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect)

pdf("Plot_SDvsEE_coefficients.pdf", width=5, height=5)
plot(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_EEbySD_RNAConcExtractDissect,  xlab="EE coefficients (log2 fold change)", ylab="SD coefficients (log2 fold change)")
BestFitLine<-lm(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_EEbySD_RNAConcExtractDissect)
abline(BestFitLine, col=1, lwd=2)
points(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_EEbySD_RNAConcExtractDissect[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE.2<0.001,], col=3, pch=16)
points(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_EEbySD_RNAConcExtractDissect[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD.2<0.001,], col=2, pch=16)
points(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_EEbySD_RNAConcExtractDissect[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD.2<0.001 & NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE.2<0.001,], col=7, pch=16)
mtext(paste("R-squared=", round(summary.lm(BestFitLine)$r.squared,digits=2), sep="")) 
dev.off()

#Volcano plots:

tiff("VolcanoPlot_SD.tiff", width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, plot(SocialDefeatSD, -log10(SocialDefeatSD.2), pch=19, main="Effect of Social Defeat", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="SD coefficient (log2 fold change)", ylab="SD -log10(p-value)"))
with(subset(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, abs(SocialDefeatSD)>1), points(SocialDefeatSD, -log10(SocialDefeatSD.2), pch=19, col="red", cex=0.6))
with(subset(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, SocialDefeatSD.2<.001), points(SocialDefeatSD, -log10(SocialDefeatSD.2), pch=19, col="green", cex=0.6))
with(subset(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, abs(SocialDefeatSD)>1 & SocialDefeatSD.2<.001), points(SocialDefeatSD, -log10(SocialDefeatSD.2), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "p < 0.001", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()


tiff("VolcanoPlot_EE.tiff", width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, plot(EnrichmentEE, -log10(EnrichmentEE.2), pch=19, main="Effect of Enrichment", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="EE coefficient (log2 fold change)", ylab="EE -log10(p-value)"))
with(subset(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, abs(EnrichmentEE)>1), points(EnrichmentEE, -log10(EnrichmentEE.2), pch=19, col="red", cex=0.6))
with(subset(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, EnrichmentEE.2<.001), points(EnrichmentEE, -log10(EnrichmentEE.2), pch=19, col="green", cex=0.6))
with(subset(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, abs(EnrichmentEE)>1 & EnrichmentEE.2<.001), points(EnrichmentEE, -log10(EnrichmentEE.2), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "p < 0.001", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()

################################################

#I'm going to run this one more time, but without the interaction term because I'm a little nervous that the correlation between the EE and SD effects is created artificially by some oddity about the LR Nil group. 
#This is also hypothetically justified by the lack of any genes showing a significant interaction term 
#Our behavioral data would suggest otherwise, as would the later cell type results, but the fact that the interaction terms seem to be such that they are typically a reversal of the similar effects of EE and SD, it suggests to me that the interaction terms are more likely to reflect some sort of regression to the mean (or lack of replication) due to the LR Nil group being odd in some way.

design <- model.matrix(~Enrichment+SocialDefeat+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=Sample_MetaData_NACC_NoOutliers_Ordered)

design

# (Intercept) EnrichmentEE SocialDefeatSD RNA_conc date_of_RNA_extraction28.2.2019
# 68           1            0              0     47.0                               0
# 69           1            0              0     67.4                               0
# 70           1            0              0    114.9                               1
# 61           1            0              0     61.3                               1
# 62           1            0              0     50.3                               1
# 63           1            0              0     92.2                               0
# 64           1            0              0    156.4                               0
# 65           1            0              0     58.0                               1
# 41           1            0              1     91.2                               0
# 42           1            0              1     40.6                               1
# 43           1            0              1     98.9                               0
# 44           1            0              1     68.8                               0
# 45           1            0              1     81.0                               0
# 46           1            0              1     97.8                               0
# 47           1            0              1     97.8                               0
# 48           1            0              1     54.2                               0
# 74           1            1              0     63.3                               0
# 75           1            1              0     50.7                               0
# 76           1            1              0     47.8                               0
# 32           1            1              0     42.7                               1
# 33           1            1              0     51.5                               0
# 34           1            1              0     56.6                               0
# 35           1            1              0     78.1                               0
# 71           1            1              1     76.6                               1
# 72           1            1              1     69.5                               1
# 73           1            1              1     55.7                               0
# 37           1            1              1     62.7                               1
# 38           1            1              1     27.7                               1
# 40           1            1              1     55.4                               0
# 66           1            0              0     81.8                               0
# 67           1            0              0     54.8                               0
# 77           1            0              1     90.2                               1
# 78           1            0              1     35.9                               1
# 79           1            0              1     36.5                               1
# 80           1            0              1     61.4                               0
# 49           1            1              0     36.1                               0
# 50           1            1              0     24.0                               1
# 51           1            1              0     89.2                               1
# 52           1            1              0    111.8                               1
# 53           1            1              0     68.3                               0
# 54           1            1              0     67.5                               0
# 55           1            1              1     93.4                               1
# 56           1            1              1     19.7                               0
# 57           1            1              1     79.1                               1
# 58           1            1              1     48.6                               1
# 60           1            1              1     72.1                               0
# date_of_dissection17/1/2019 date_of_dissection21/2/2019
# 68                           0                           0
# 69                           0                           0
# 70                           0                           0
# 61                           0                           0
# 62                           0                           0
# 63                           0                           0
# 64                           0                           0
# 65                           0                           0
# 41                           0                           0
# 42                           1                           0
# 43                           0                           0
# 44                           0                           0
# 45                           1                           0
# 46                           1                           0
# 47                           1                           0
# 48                           0                           0
# 74                           0                           0
# 75                           0                           0
# 76                           0                           0
# 32                           0                           0
# 33                           1                           0
# 34                           0                           0
# 35                           0                           0
# 71                           0                           0
# 72                           0                           0
# 73                           0                           0
# 37                           0                           0
# 38                           1                           0
# 40                           1                           0
# 66                           0                           1
# 67                           0                           1
# 77                           0                           1
# 78                           0                           1
# 79                           0                           1
# 80                           0                           1
# 49                           0                           1
# 50                           0                           1
# 51                           0                           1
# 52                           0                           1
# 53                           0                           1
# 54                           0                           1
# 55                           0                           1
# 56                           0                           1
# 57                           0                           1
# 58                           0                           1
# 60                           0                           1
# attr(,"assign")
# [1] 0 1 2 3 4 5 5
# attr(,"contrasts")
# attr(,"contrasts")$Enrichment
# [1] "contr.treatment"
# 
# attr(,"contrasts")$SocialDefeat
# [1] "contr.treatment"
# 
# attr(,"contrasts")$date_of_RNA_extraction
# [1] "contr.treatment"
# 
# attr(,"contrasts")$date_of_dissection
# [1] "contr.treatment"

v <- voom(NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("JustEESD_Cov_RNAConcExtractDissection_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD RNA_conc date_of_RNA_extraction28.2.2019
# -1        3777            5              0        1                               0
# 0         1100        17760          17761    17757                           17765
# 1        12888            0              4        7                               0
# date_of_dissection17/1/2019 date_of_dissection21/2/2019
# -1                           0                         907
# 0                        17765                       15953
# 1                            0                         905


topTable(efit, coef=2)

#                         logFC    AveExpr         t      P.Value  adj.P.Val        B
# ENSRNOG00000054953 -1.0491970  1.2909441 -5.628825 1.220175e-06 0.01416944 5.246229
# ENSRNOG00000034126 -0.5316513  4.2630641 -5.548778 1.595209e-06 0.01416944 4.976021
# ENSRNOG00000032884 -2.2359303  0.1491942 -5.290316 3.777315e-06 0.02236800 3.954280
# ENSRNOG00000020084 -0.3075369  3.3187582 -5.041666 8.602478e-06 0.03056460 3.444040
# ENSRNOG00000058611 -0.4741330  1.6293954 -4.736744 2.334056e-05 0.06910750 2.598366
# ENSRNOG00000018177 -1.4875051 -1.6998426 -5.115406 6.744436e-06 0.02995373 2.511993
# ENSRNOG00000052687 -0.2700643  8.2508644 -4.598559 3.650876e-05 0.09265401 1.775663
# ENSRNOG00000047581  0.7973627  0.7771356  4.463459 5.634112e-05 0.09786975 1.757809
# ENSRNOG00000020073  0.4481302  3.3263494  4.478444 5.370369e-05 0.09786975 1.700259
# ENSRNOG00000008759  0.2878583  1.9773957  4.413376 6.610960e-05 0.09786975 1.634802

#The top genes are still a bunch of protocadherins:
# Scn11a
# Pcdhb6
# Pcdhga2
# Dhrs2
# Pcdhb5
# AABR07027137.1
# Megf8
# Pcdhb8
# AABR07066700.1
# Ccl3
# Csf3r
# Mbd1

topTable(efit, coef=3)

#                         logFC    AveExpr         t      P.Value  adj.P.Val         B
# ENSRNOG00000029386  0.7288571  3.4069905  5.643825 1.160352e-06 0.01823770 5.3728736
# ENSRNOG00000039744  0.4885751  1.6609750  5.473269 2.053217e-06 0.01823770 4.8012913
# ENSRNOG00000000512  0.5953955  4.6896305  4.965120 1.106619e-05 0.04914772 3.0930014
# ENSRNOG00000000723  0.9560568 -0.9603081  5.050344 8.359919e-06 0.04914772 2.6748390
# ENSRNOG00000001854 -0.2303343  5.7485111 -4.804911 1.869561e-05 0.06642549 2.5191091
# ENSRNOG00000003203  0.1059319  6.0181938  4.380285 7.345378e-05 0.21748441 1.1818284
# ENSRNOG00000047706 -1.2884466  0.3144700 -4.205845 1.274626e-04 0.25159698 0.9435905
# ENSRNOG00000061508 -0.1927740  5.2810710 -4.240094 1.144551e-04 0.25159698 0.8032221
# ENSRNOG00000012301  0.5799630  2.0376598  3.987156 2.516298e-04 0.26101225 0.4164571
# ENSRNOG00000017319 -0.2494767  4.9482079 -4.071320 1.939743e-04 0.26101225 0.3273488

#The top genes are still a lot of RT-1s:
# RT1-N2
# RT1-CE4
# Slc26a8
# RT1-CE5
# Tmtc1

write.fit(efit, adjust="BH", file="Limma_results_Model_onlyEESD_RNAconc_RNAextract_Dissect.txt")

#write.csv(NACC_RNASeq_Annotation_noLowHits, "NACC_RNASeq_Annotation_noLowHits.csv")

NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect<-data.frame(NACC_RNASeq_Log2_Annotated$gene_id, NACC_RNASeq_Log2_Annotated$gene_symbol, efit$coefficients, efit$t, efit$p.value, stringsAsFactors = F)
str(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect)

write.csv(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, "NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect.csv")

########################


#Looking at the relationship between variables:

str(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect)

pdf("Plot_SDvsEE_coefficients_noInteraction.pdf", width=5, height=5)
plot(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect,  xlab="EE coefficients (log2 fold change)", ylab="SD coefficients (log2 fold change)")
BestFitLine<-lm(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect)
abline(BestFitLine, col=1, lwd=2)
points(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001,], col=3, pch=16)
points(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001,], col=2, pch=16)
points(SocialDefeatSD~EnrichmentEE, data=NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001 & NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001,], col=7, pch=16)
mtext(paste("R-squared=", round(summary.lm(BestFitLine)$r.squared,digits=2), sep="")) 
dev.off()
#Yep, there is *no* relationship between the two variables once they do not share a control group. 
#That is disconcerting, but perhaps not unexpected

#Volcano plots:

tiff("VolcanoPlot_SD_noInteraction.tiff", width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, plot(SocialDefeatSD, -log10(SocialDefeatSD.2), pch=19, main="Effect of Social Defeat", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="SD coefficient (log2 fold change)", ylab="SD -log10(p-value)"))
with(subset(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, abs(SocialDefeatSD)>1), points(SocialDefeatSD, -log10(SocialDefeatSD.2), pch=19, col="red", cex=0.6))
with(subset(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, SocialDefeatSD.2<.001), points(SocialDefeatSD, -log10(SocialDefeatSD.2), pch=19, col="green", cex=0.6))
with(subset(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, abs(SocialDefeatSD)>1 & SocialDefeatSD.2<.001), points(SocialDefeatSD, -log10(SocialDefeatSD.2), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "p < 0.001", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()


tiff("VolcanoPlot_EE_noInteraction.tiff", width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, plot(EnrichmentEE, -log10(EnrichmentEE.2), pch=19, main="Effect of Enrichment", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="EE coefficient (log2 fold change)", ylab="EE -log10(p-value)"))
with(subset(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, abs(EnrichmentEE)>1), points(EnrichmentEE, -log10(EnrichmentEE.2), pch=19, col="red", cex=0.6))
with(subset(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, EnrichmentEE.2<.001), points(EnrichmentEE, -log10(EnrichmentEE.2), pch=19, col="green", cex=0.6))
with(subset(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, abs(EnrichmentEE)>1 & EnrichmentEE.2<.001), points(EnrichmentEE, -log10(EnrichmentEE.2), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "p < 0.001", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()

#The fold changes for both variables are definitely a lot more moderate now.

################################################

#And one more time: Running the model without the interaction term with fewer co-variates to check for robustness:

design <- model.matrix(~Enrichment+SocialDefeat+date_of_dissection, data=Sample_MetaData_NACC_NoOutliers_Ordered)

design

# (Intercept) EnrichmentEE SocialDefeatSD date_of_dissection17/1/2019 date_of_dissection21/2/2019
# 68           1            0              0                           0                           0
# 69           1            0              0                           0                           0
# 70           1            0              0                           0                           0
# 61           1            0              0                           0                           0
# 62           1            0              0                           0                           0
# 63           1            0              0                           0                           0
# 64           1            0              0                           0                           0
# 65           1            0              0                           0                           0
# 41           1            0              1                           0                           0
# 42           1            0              1                           1                           0
# 43           1            0              1                           0                           0
# 44           1            0              1                           0                           0
# 45           1            0              1                           1                           0
# 46           1            0              1                           1                           0
# 47           1            0              1                           1                           0
# 48           1            0              1                           0                           0
# 74           1            1              0                           0                           0
# 75           1            1              0                           0                           0
# 76           1            1              0                           0                           0
# 32           1            1              0                           0                           0
# 33           1            1              0                           1                           0
# 34           1            1              0                           0                           0
# 35           1            1              0                           0                           0
# 71           1            1              1                           0                           0
# 72           1            1              1                           0                           0
# 73           1            1              1                           0                           0
# 37           1            1              1                           0                           0
# 38           1            1              1                           1                           0
# 40           1            1              1                           1                           0
# 66           1            0              0                           0                           1
# 67           1            0              0                           0                           1
# 77           1            0              1                           0                           1
# 78           1            0              1                           0                           1
# 79           1            0              1                           0                           1
# 80           1            0              1                           0                           1
# 49           1            1              0                           0                           1
# 50           1            1              0                           0                           1
# 51           1            1              0                           0                           1
# 52           1            1              0                           0                           1
# 53           1            1              0                           0                           1
# 54           1            1              0                           0                           1
# 55           1            1              1                           0                           1
# 56           1            1              1                           0                           1
# 57           1            1              1                           0                           1
# 58           1            1              1                           0                           1
# 60           1            1              1                           0                           1
# attr(,"assign")
# [1] 0 1 2 3 3
# attr(,"contrasts")
# attr(,"contrasts")$Enrichment
# [1] "contr.treatment"
# 
# attr(,"contrasts")$SocialDefeat
# [1] "contr.treatment"
# 
# attr(,"contrasts")$date_of_dissection
# [1] "contr.treatment"

v <- voom(NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("JustEESD_Cov_Dissection_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD date_of_dissection17/1/2019 date_of_dissection21/2/2019
# -1        4237            5              1                           0                        1062
# 0          527        17760          17760                       17765                       15704
# 1        13001            0              4                           0                         999


topTable(efit, coef=2)

# logFC    AveExpr         t      P.Value  adj.P.Val        B
# ENSRNOG00000034126 -0.5240861  4.2630641 -5.487675 1.740237e-06 0.02256283 4.879211
# ENSRNOG00000054953 -0.9932969  1.2909441 -5.238720 4.044308e-06 0.02256283 4.158291
# ENSRNOG00000032884 -2.1387615  0.1491942 -5.171035 5.080288e-06 0.02256283 3.707311
# ENSRNOG00000020084 -0.2959155  3.3187582 -5.066703 7.211954e-06 0.02562407 3.599543
# ENSRNOG00000018177 -1.4872643 -1.6998426 -5.354715 2.732495e-06 0.02256283 3.181206
# ENSRNOG00000058611 -0.4638156  1.6293954 -4.756883 2.021710e-05 0.05985947 2.722668
# ENSRNOG00000047581  0.7884975  0.7771356  4.557685 3.887071e-05 0.08631728 2.083612
# ENSRNOG00000008759  0.2790663  1.9773957  4.480636 4.994403e-05 0.09858396 1.882509
# ENSRNOG00000052687 -0.2622196  8.2508644 -4.588169 3.518844e-05 0.08631728 1.808815
# ENSRNOG00000000512  0.4946815  4.6896305  4.440667 5.684997e-05 0.10099397 1.467092

#An almost identical list to what we found with the model with all co-variates

topTable(efit, coef=3)

#                       logFC    AveExpr         t      P.Value  adj.P.Val         B
# ENSRNOG00000029386  0.7155352  3.4069905  5.647691 1.008940e-06 0.01792382 5.5009191
# ENSRNOG00000039744  0.4794443  1.6609750  5.267565 3.669107e-06 0.02155536 4.2564063
# ENSRNOG00000000512  0.6075125  4.6896305  5.184606 4.853445e-06 0.02155536 3.8933357
# ENSRNOG00000000723  0.9640154 -0.9603081  5.211521 4.432760e-06 0.02155536 3.1030266
# ENSRNOG00000001854 -0.2301066  5.7485111 -4.920840 1.173958e-05 0.04171074 2.9798032
# ENSRNOG00000003203  0.1042834  6.0181938  4.387706 6.745495e-05 0.17715460 1.2785482
# ENSRNOG00000047706 -1.2701816  0.3144700 -4.234893 1.100764e-04 0.21727849 1.0461109
# ENSRNOG00000061508 -0.1907008  5.2810710 -4.268262 9.896239e-05 0.21727849 0.9561600
# ENSRNOG00000017319 -0.2508053  4.9482079 -4.185531 1.287778e-04 0.22045841 0.7317116
# ENSRNOG00000012301  0.5787940  2.0376598  4.074441 1.828728e-04 0.26472411 0.707471

#Yep, this list is almost identical to what we found with all covariates as well.

write.fit(efit, adjust="BH", file="Limma_results_Model_onlyEESD_Dissect.txt")

#And with no covariates:

design <- model.matrix(~Enrichment+SocialDefeat, data=Sample_MetaData_NACC_NoOutliers_Ordered)

design

# (Intercept) EnrichmentEE SocialDefeatSD
# 68           1            0              0
# 69           1            0              0
# 70           1            0              0
# 61           1            0              0
# 62           1            0              0
# 63           1            0              0
# 64           1            0              0
# 65           1            0              0
# 41           1            0              1
# 42           1            0              1
# 43           1            0              1
# 44           1            0              1
# 45           1            0              1
# 46           1            0              1
# 47           1            0              1
# 48           1            0              1
# 74           1            1              0
# 75           1            1              0
# 76           1            1              0
# 32           1            1              0
# 33           1            1              0
# 34           1            1              0
# 35           1            1              0
# 71           1            1              1
# 72           1            1              1
# 73           1            1              1
# 37           1            1              1
# 38           1            1              1
# 40           1            1              1
# 66           1            0              0
# 67           1            0              0
# 77           1            0              1
# 78           1            0              1
# 79           1            0              1
# 80           1            0              1
# 49           1            1              0
# 50           1            1              0
# 51           1            1              0
# 52           1            1              0
# 53           1            1              0
# 54           1            1              0
# 55           1            1              1
# 56           1            1              1
# 57           1            1              1
# 58           1            1              1
# 60           1            1              1
# attr(,"assign")
# [1] 0 1 2
# attr(,"contrasts")
# attr(,"contrasts")$Enrichment
# [1] "contr.treatment"
# 
# attr(,"contrasts")$SocialDefeat
# [1] "contr.treatment"

v <- voom(NACC_RNASeq_dge_TMM, design, plot=TRUE)
v

vfit <- lmFit(v, design)
str(vfit)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
str(efit)

png("JustEESD_MeanVarianceTrend.png")
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()

dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) EnrichmentEE SocialDefeatSD
# -1        4277            5              7
# 0          481        17758          17750
# 1        13007            2              8

topTable(efit, coef=2)
# 
#                         logFC    AveExpr         t      P.Value  adj.P.Val        B
# ENSRNOG00000034126 -0.5092704  4.2630641 -5.416680 2.005076e-06 0.01307206 4.746742
# ENSRNOG00000054953 -0.9811978  1.2909441 -5.388266 2.211137e-06 0.01307206 4.684832
# ENSRNOG00000020084 -0.2949332  3.3187582 -5.337181 2.635698e-06 0.01307206 4.561142
# ENSRNOG00000047581  0.8859291  0.7771356  5.305027 2.943329e-06 0.01307206 4.338905
# ENSRNOG00000000512  0.5621524  4.6896305  4.868029 1.300390e-05 0.03850238 2.902246
# ENSRNOG00000018177 -1.3419987 -1.6998426 -4.908212 1.135793e-05 0.03850238 2.085687
# ENSRNOG00000029427 -0.6530607  0.5253441 -4.534649 3.948265e-05 0.08767617 2.021572
# ENSRNOG00000020073  0.4077232  3.3263494  4.451081 5.195961e-05 0.10256250 1.720813
# ENSRNOG00000008759  0.2583222  1.9773957  4.347195 7.292775e-05 0.12955614 1.530841
# ENSRNOG00000061359 -2.7009679 -2.4802048 -4.779814 1.748432e-05 0.04437272 1.189930

#Similar to the version with all co-variates, there are a few switches.

topTable(efit, coef=3)

# logFC  AveExpr         t      P.Value    adj.P.Val        B
# ENSRNOG00000029386  0.7283011 3.406991  6.419635 6.096222e-08 0.0005414969 8.161275
# ENSRNOG00000039744  0.5295106 1.660975  6.434928 5.778041e-08 0.0005414969 8.098309
# ENSRNOG00000059746 -0.5237156 2.860831 -5.256110 3.480700e-06 0.0206115447 4.320365
# ENSRNOG00000047706 -1.3365745 0.314470 -5.116481 5.607539e-06 0.0249044841 3.680998
# ENSRNOG00000000512  0.5646852 4.689631  4.887592 1.217517e-05 0.0432583694 2.950985
# ENSRNOG00000011617  0.2218770 3.891533  4.760559 1.864782e-05 0.0446823476 2.618598
# ENSRNOG00000021536  0.3766528 4.021490  4.670322 2.519528e-05 0.0446823476 2.316260
# ENSRNOG00000017319 -0.2544681 4.948208 -4.692434 2.340784e-05 0.0446823476 2.298914
# ENSRNOG00000011774  0.8374078 1.495489  4.581114 3.386825e-05 0.0446823476 2.243310
# ENSRNOG00000017869  0.2954544 3.775892  4.635786 2.825829e-05 0.0446823476 2.233987

#Similar to the version wtih all co-variates, but definitely has some switches.

write.fit(efit, adjust="BH", file="Limma_results_Model_onlyEESD.txt")

NAcc_LimmaResults_onlyEESD<-data.frame(NACC_RNASeq_Log2_Annotated$gene_id, NACC_RNASeq_Log2_Annotated$gene_symbol, efit$coefficients, efit$t, efit$p.value, stringsAsFactors = F)
str(NAcc_LimmaResults_onlyEESD)

write.csv(NAcc_LimmaResults_onlyEESD, "NAcc_LimmaResults_onlyEESD_noCov.csv")

################################################


#I would like a better gut sense of how the presence or absence of covariates is affecting the results. Is it just reducing power or actually changing the nature of what we are finding?


NAcc_LimmaResults_onlyEESD
NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect

pdf("Plot_EEvsEEWCov_coefficients_noInteraction.pdf", width=5, height=5)
plot(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE~NAcc_LimmaResults_onlyEESD$EnrichmentEE, xlab="EE coefficients: no covariates", ylab="EE coefficients: w/ covariates")
BestFitLine<-lm(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE~NAcc_LimmaResults_onlyEESD$EnrichmentEE)
abline(BestFitLine, col=1, lwd=2)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001]~NAcc_LimmaResults_onlyEESD$EnrichmentEE[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001], col=3, pch=16)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE[NAcc_LimmaResults_onlyEESD$EnrichmentEE.2<0.001]~NAcc_LimmaResults_onlyEESD$EnrichmentEE[NAcc_LimmaResults_onlyEESD$EnrichmentEE.2<0.001], col=2, pch=16)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001 & NAcc_LimmaResults_onlyEESD$EnrichmentEE.2<0.001]~NAcc_LimmaResults_onlyEESD$EnrichmentEE[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001 & NAcc_LimmaResults_onlyEESD$EnrichmentEE.2<0.001], col=7, pch=16)
mtext(paste("R-squared=", round(summary.lm(BestFitLine)$r.squared,digits=2), sep="")) 
dev.off()

pdf("Plot_SDvsSDWCov_coefficients_noInteraction.pdf", width=5, height=5)
plot(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD~NAcc_LimmaResults_onlyEESD$SocialDefeatSD, xlab="SD coefficients: no covariates", ylab="SD coefficients: w/ covariates")
BestFitLine<-lm(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD~NAcc_LimmaResults_onlyEESD$SocialDefeatSD)
abline(BestFitLine, col=1, lwd=2)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001]~NAcc_LimmaResults_onlyEESD$SocialDefeatSD[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001], col=3, pch=16)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD[NAcc_LimmaResults_onlyEESD$SocialDefeatSD.2<0.001]~NAcc_LimmaResults_onlyEESD$SocialDefeatSD[NAcc_LimmaResults_onlyEESD$SocialDefeatSD.2<0.001], col=2, pch=16)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001 & NAcc_LimmaResults_onlyEESD$SocialDefeatSD.2<0.001]~NAcc_LimmaResults_onlyEESD$SocialDefeatSD[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001 & NAcc_LimmaResults_onlyEESD$SocialDefeatSD.2<0.001], col=7, pch=16)
mtext(paste("R-squared=", round(summary.lm(BestFitLine)$r.squared,digits=2), sep="")) 
dev.off()

#Interesting - the presence or absence of covariates in the model seems to matter more for the social defeat results vs. the enrichment results.
#Probably because social defeat (when considered alone) is almost collinear with dissection


#Let's plot a comparison between the results with or without an interaction term:


NAcc_LimmaResults_EEbySD_RNAConcExtractDissect

pdf("Plot_EENoInteractvsEEwInteract_coefficients_wCovariates.pdf", width=5, height=5)
plot(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE, xlab="EE coefficients: with interaction term", ylab="EE coefficients: no interaction term")
BestFitLine<-lm(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE)
abline(BestFitLine, col=1, lwd=2)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001]~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001], col=3, pch=16)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE.2<0.001]~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE.2<0.001], col=2, pch=16)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001 & NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE.2<0.001]~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$EnrichmentEE.2<0.001 & NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$EnrichmentEE.2<0.001], col=7, pch=16)
mtext(paste("R-squared=", round(summary.lm(BestFitLine)$r.squared,digits=2), sep="")) 
dev.off()

pdf("Plot_SDNoInteractvsSDwInteract_coefficients_wCovariates.pdf", width=5, height=5)
plot(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD, xlab="SD coefficients: with interaction term", ylab="SD coefficients: no interaction term")
BestFitLine<-lm(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD)
abline(BestFitLine, col=1, lwd=2)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001]~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001], col=3, pch=16)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD.2<0.001]~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD.2<0.001], col=2, pch=16)
points(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001 & NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD.2<0.001]~NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect$SocialDefeatSD.2<0.001 & NAcc_LimmaResults_EEbySD_RNAConcExtractDissect$SocialDefeatSD.2<0.001], col=7, pch=16)
mtext(paste("R-squared=", round(summary.lm(BestFitLine)$r.squared,digits=2), sep="")) 
dev.off()

#EE is much more sensitive to the presence or absence of an interaction term.


#################################################
#Plot some of the top results to confirm direction of effect, etc.
#I should also do a mediocre attempt to make sure that the effects aren't due to batch confounds.

head(NACC_RNASeq_Log2_Annotated)
colnames(NACC_RNASeq_Log2_Annotated)

head(NACC_RNASeq_Log2_Annotated$gene_symbol)

GroupingVariable<-Sample_MetaData_NACC_NoOutliers_Ordered$Treatment_group
table(GroupingVariable)
levels(GroupingVariable)
GroupingVariable_Factor<-factor(GroupingVariable, levels=c("bLR NIL", "bLR EE", "bLR NIL + SD", "bLR EE + SD"))
levels(GroupingVariable_Factor)

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhb6"),]

pdf("Pcdhb6_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhb6 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="AC120486.2"),]

pdf("AC120486_2_Prettier_vs_Group_v3.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="AC120486.2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhga2"),]

pdf("Pcdhga2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhga2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhb5"),]

pdf("Pcdhb5_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhb5 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhb8"),]

pdf("Pcdhb8_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhb8 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Dele1"),]

pdf("Dele1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Dele1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhb7"),]

pdf("Pcdhb7_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhb7 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="RT1-N2"),]

pdf("RT1-N2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="RT1-N2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="RT1-CE4"),]

pdf("RT1-CE4_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="RT1-CE4 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="RT1-CE5"),]

pdf("RT1-CE5_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="RT1-CE5 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Slc26a8"),]

pdf("Slc26a8_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Slc26a8 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Tmtc1"),]

pdf("Tmtc1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Tmtc1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Scn11a"),]

pdf("Scn11a_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Scn11a Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Dhrs2"),]

pdf("Dhrs2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Dhrs2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="AABR07027137.1"),]

pdf("AABR07027137.1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="AABR07027137.1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Megf8"),]

pdf("Megf8_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Megf8 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Scn11a"),]

pdf("Scn11a_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Scn11a Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Pcdh17"),]

pdf("Pcdh17_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdh17 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Abca12"),]

pdf("Abca12_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Abca12 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


#Top HRLR DE genes that also show effects of SD and EE 

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Tmem144"),]

pdf("Tmem144_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Tmem144 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Pcdhga1"),]

pdf("Pcdhga1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Pcdhga1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Etv4"),]

pdf("Etv4_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Etv4 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()
#I wonder why that doesn't seem to show any effects?


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Uhrf1"),]

pdf("Uhrf1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Uhrf1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Tdg"),]

pdf("Tdg_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Tdg Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Tubg1"),]

pdf("Tubg1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Tubg1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Slc19a3"),]

pdf("Slc19a3_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Slc19a3 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Prss55"),]

pdf("Prss55_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Prss55 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()



#For NIDA-related stuff:

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Oprm1"),]

pdf("Oprm1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Oprm1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()



GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Hba-a2"),][1,]

pdf("Hba-a2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Hba-a2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Hba-a2"),][2,]

pdf("Hba-a2v2_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Hba-a2 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Hbb"),]

pdf("Hbb_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Hbb Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Htra1"),]

pdf("Htra1_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Htra1 Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-NACC_RNASeq_Log2_Filtered[which(NACC_RNASeq_Log2_Annotated$gene_symbol=="Fos"),]

pdf("Fos_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Fos Log2 Cpm", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()

###########################

#Cell type analyses:


library(BrainInABlender)

#Hopefully the extremely minimal filtering used with the data doesn't mess this up too much - we may want to try re-doing the estimates with a more filtered dataset.

Temp<-data.frame(NACC_RNASeq_Log2_Annotated$gene_symbol, NACC_RNASeq_Log2_Filtered)
colnames(Temp)
CellTypeOutput<-Sir_UnMixALot(Temp, dataColumns = c(2:47), geneColumn = 1, species="mouse")
rm(Temp)


summary.lm(lm(y~Enrichment*SocialDefeat+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=Temp))


#Using the model with an EE*SD interaction term and covariates:

for(i in c(1:10)){
  Temp<-data.frame(y=CellTypeOutput$AveragePrimary_CellTypeIndex[i,], Sample_MetaData_NACC_NoOutliers_Ordered)
  Model<-lm(y~Enrichment*SocialDefeat+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=Temp)
  OutputtingStats<-file("NAcc_LM_CellType_EEbySD_GenRNAConcExtract.txt")
  stats_output <- c(
    print(row.names(CellTypeOutput$AveragePrimary_CellTypeIndex)[i]),
    capture.output(summary.lm(Model)),
    print("******************************************************************************")
  )
  cat(stats_output, file="NAcc_LM_CellType_EEbySD_GenRNAConcExtract.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#date_of_dissection21/2/2019 appears to have less astrocytes, less RBC, and more neurons

#The only cell type robustly affected by EE or SD is the general neuron signature... which is unfortunately somewhat inadequate for this brain region. Similar to the expression results, EE and SD cause changes in the same direction (decrease) and the decrease caused by SD is diminished by EE.

#Enrichment also seems to potentially increase RBC content

# Neuron_All
# 
# Call:
#   lm(formula = y ~ Enrichment * SocialDefeat + RNA_conc + date_of_RNA_extraction + 
#        date_of_dissection, data = Temp)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.95862 -0.12919  0.03594  0.13704  0.61432 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                     -0.362525   0.169473  -2.139 0.038912 *  
#   EnrichmentEE                    -0.162390   0.128334  -1.265 0.213441    
# SocialDefeatSD                  -0.358016   0.133581  -2.680 0.010821 *  
#   RNA_conc                         0.005476   0.001701   3.219 0.002636 ** 
#   date_of_RNA_extraction28.2.2019  0.017703   0.091062   0.194 0.846894    
# date_of_dissection17/1/2019      0.110212   0.135787   0.812 0.422046    
# date_of_dissection21/2/2019      0.378531   0.097914   3.866 0.000419 ***
#   EnrichmentEE:SocialDefeatSD      0.396724   0.178365   2.224 0.032150 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2871 on 38 degrees of freedom
# Multiple R-squared:  0.4521,	Adjusted R-squared:  0.3512 
# F-statistic:  4.48 on 7 and 38 DF,  p-value: 0.001023


# RBC
# 
# Call:
#   lm(formula = y ~ Enrichment * SocialDefeat + RNA_conc + date_of_RNA_extraction + 
#        date_of_dissection, data = Temp)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.8013 -0.3301 -0.0228  0.2130  1.0720 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      0.081096   0.289678   0.280 0.781032    
# EnrichmentEE                     0.470718   0.219359   2.146 0.038333 *  
#   SocialDefeatSD                  -0.140106   0.228328  -0.614 0.543125    
# RNA_conc                         0.002955   0.002908   1.016 0.316044    
# date_of_RNA_extraction28.2.2019 -0.173636   0.155650  -1.116 0.271619    
# date_of_dissection17/1/2019     -0.412467   0.232098  -1.777 0.083554 .  
# date_of_dissection21/2/2019     -0.714224   0.167363  -4.268 0.000127 ***
#   EnrichmentEE:SocialDefeatSD     -0.236021   0.304876  -0.774 0.443630    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4908 on 38 degrees of freedom
# Multiple R-squared:  0.4743,	Adjusted R-squared:  0.3774 
# F-statistic: 4.897 on 7 and 38 DF,  p-value: 0.0005194



#Using the model without an interaction term, and covariates:
#Note: this doesn't do anything to account for non-normal distribution (unlike limma). fix?

for(i in c(1:10)){
  Temp<-data.frame(y=CellTypeOutput$AveragePrimary_CellTypeIndex[i,], Sample_MetaData_NACC_NoOutliers_Ordered)
  Model<-lm(y~Enrichment+SocialDefeat+RNA_conc+date_of_RNA_extraction+date_of_dissection, data=Temp)
  OutputtingStats<-file("NAcc_LM_CellType_onlyEESD_DissectRNAConcExtract.txt")
  stats_output <- c(
    print(row.names(CellTypeOutput$AveragePrimary_CellTypeIndex)[i]),
    capture.output(summary.lm(Model)),
    print("******************************************************************************")
  )
  cat(stats_output, file="NAcc_LM_CellType_onlyEESD_GenRNAConcExtract.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#note: the 21/2 date of dissection (F49 generation) has a really different cell type balance than the other two


# ******************************************************************************
#   Oligodendrocyte_Immature
# 
# Call:
#   lm(formula = y ~ Enrichment + SocialDefeat + RNA_conc + date_of_RNA_extraction + 
#        date_of_dissection, data = Temp)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.44349 -0.15301  0.01231  0.15786  0.44076 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                      0.041434   0.130069   0.319   0.7518  
# EnrichmentEE                    -0.056549   0.073659  -0.768   0.4473  
# SocialDefeatSD                  -0.181950   0.075015  -2.426   0.0200 *
#   RNA_conc                         0.002292   0.001390   1.649   0.1072  
# date_of_RNA_extraction28.2.2019 -0.160207   0.072944  -2.196   0.0341 *
#   date_of_dissection17/1/2019     -0.148306   0.109264  -1.357   0.1825  
# date_of_dissection21/2/2019      0.037964   0.079079   0.480   0.6339  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.235 on 39 degrees of freedom
# Multiple R-squared:  0.3733,	Adjusted R-squared:  0.2769 
# F-statistic: 3.872 on 6 and 39 DF,  p-value: 0.003985
# 
# ******************************************************************************
#   RBC
# 
# Call:
#   lm(formula = y ~ Enrichment + SocialDefeat + RNA_conc + date_of_RNA_extraction + 
#        date_of_dissection, data = Temp)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.84952 -0.33863 -0.04658  0.21681  1.10375 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      0.159054   0.270212   0.589 0.559507    
# EnrichmentEE                     0.349644   0.153022   2.285 0.027839 *  
#   SocialDefeatSD                  -0.268708   0.155839  -1.724 0.092578 .  
# RNA_conc                         0.002813   0.002887   0.974 0.335911    
# date_of_RNA_extraction28.2.2019 -0.198418   0.151538  -1.309 0.198074    
# date_of_dissection17/1/2019     -0.379527   0.226990  -1.672 0.102530    
# date_of_dissection21/2/2019     -0.693145   0.164283  -4.219 0.000141 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4883 on 39 degrees of freedom
# Multiple R-squared:  0.466,	Adjusted R-squared:  0.3838 
# F-statistic: 5.672 on 6 and 39 DF,  p-value: 0.0002604
# 
# ******************************************************************************


#I'm hesitant to remove any of the co-variates, because they all seem to matter in this output.
#Although out of all of them, RNAconc seems to matter the least.


#Recycling some code to make graphs:

row.names(CellTypeOutput$AveragePrimary_CellTypeIndex)

GeneY<-CellTypeOutput$AveragePrimary_CellTypeIndex[10,]

pdf("RBC_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="RBC index", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-CellTypeOutput$AveragePrimary_CellTypeIndex[5,]

pdf("NeuronAll_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Neuron_All index", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


GeneY<-CellTypeOutput$AveragePrimary_CellTypeIndex[9,]

pdf("OligodendrocyteImmature_Prettier_vs_Group_v4.pdf", width=7, height=7)
boxplot(GeneY~GroupingVariable_Factor, ylab="Immature Oligodendrocyte index", col=c("darkorange2", "burlywood1", "red", "pink"))
stripchart(GeneY~GroupingVariable_Factor, vertical = TRUE,  method = "jitter", add = TRUE, pch = c(20, 1, 17, 2), cex=2, cex.axis=10, cex.lab=10, col = 'black')
dev.off()


dim(CellTypeOutput$AveragePrimary_CellTypeIndex)
#[1] 10 46
dim(Sample_BehavHormonalData_NACC_NoOutliers_Ordered)
#[1] 46 59

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/Detailed_CellTypeAnalysisOutput")

write.csv(cor(cbind(PC1, PC2, PC3, PC4, t(CellTypeOutput$AveragePrimary_CellTypeIndex), as.matrix(Sample_MetaData_NACC_NoOutliers_Ordered[, c(13:14, 16)]), as.matrix(Sample_BehavHormonalData_NACC_NoOutliers_Ordered[,c(10:21,35:46, 48:51)])), use="pairwise.complete.obs"), "PCA_vsNumericVar_CorMatrix.csv")
#Note:  many of these variables are really really not distributed normally - let's try running a spearman instead.

write.csv(cor(cbind(PC1, PC2, PC3, PC4, t(CellTypeOutput$AveragePrimary_CellTypeIndex), as.matrix(Sample_MetaData_NACC_NoOutliers_Ordered[, c(13:14, 16)]), as.matrix(Sample_BehavHormonalData_NACC_NoOutliers_Ordered[,c(10:21,35:46, 48:51)])), use="pairwise.complete.obs", method="spearman"), "PCA_vsNumericVar_CorMatrix_Spearman.csv")

#As expected, the top principal components of variation are strongly correlated with cell type.

row.names(CellTypeOutput$AveragePrimary_CellTypeIndex)

pdf("MuralvsPC1.pdf", width=5, height=6)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[4,]~PC1, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Mural Index", xlab="PC1")
dev.off()

pdf("InterneuronvsPC1.pdf", width=5, height=6)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[6,]~PC1, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Interneuron Index", xlab="PC1")
dev.off()

pdf("NeuronvsPC1.pdf", width=5, height=6)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[5,]~PC1, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Neuron Index", xlab="PC1")
dev.off()

pdf("AstrocytevsPC2.pdf", width=5, height=6)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[1,]~PC2, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Astrocyte Index", xlab="PC2")
dev.off()

pdf("MicrogliavsPC2.pdf", width=5, height=6)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[3,]~PC2, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Microglia Index", xlab="PC2")
dev.off()

pdf("OligodendrocytevsPC2.pdf", width=5, height=6)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[8,]~PC2, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Oligodendrocyte Index", xlab="PC2")
dev.off()

pdf("NeuronvsPC2.pdf", width=5, height=6)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[5,]~PC2, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Neuron Index", xlab="PC2")
dev.off()

pdf("ProjectionNeuronvsPC2.pdf", width=5, height=6)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[7,]~PC2, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Projection Neuron Index", xlab="PC2")
dev.off()

#Anytime NIL looks totally different from the other groups, it could be these cell types that differ by dissection date. Controlling for that particular dissection date is really important.
#but maybe not 1/16 vs. 1/17?
#Meh.


#The cell type indices seem strongly related to behavior in the basic correlation plots - artifact due to technical variables? (e.g., dissection group)


pdf("MuralvsCORT.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[4,]~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Mural Cell Index", xlab="Corticosterone")
dev.off()
#pretty.

pdf("EndothelialvsCORT.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[2,]~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Endothelial Index", xlab="Corticosterone")
dev.off()
#reasonably convincing, but has an outlier.

pdf("NeuronvsCORT.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[5,]~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$CORT, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Neuron Index", xlab="Corticosterone")
dev.off()
#completely unconvincing

pdf("MicrogliavsIL6.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[3,]~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$IL.6..pg.ml., col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Microglia Index", xlab="IL6 (pg/mL)")
dev.off()
#not particularly convincing

pdf("OligodendrocytesvsOxytocin.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[8,]~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$Oxytocin..pg.ml., col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Oligodendrocyte Index", xlab="Oxytocin (pg/mL)")
dev.off()
#There actually could be something there.


#Note: Both Mural and Endothelial cell indices correlate with counts per sample 

pdf("MuralvsCountsPerSample.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[4,]~Sample_MetaData_NACC_NoOutliers_Ordered$CountsPerSample, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Mural Cell Index", xlab="Counts per Sample")
dev.off()

pdf("EndothelialvsCountsPerSample.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[2,]~Sample_MetaData_NACC_NoOutliers_Ordered$CountsPerSample, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Endothelial Index", xlab="Counts per Sample")
dev.off()


#Correlations with behavior

pdf("AstrocytevsDistance.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[1,]~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$distance_open_field, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Astrocyte Index", xlab="Distance in the Open Field")
dev.off()
#Big dissection group effect driving the relationship due to last dissection group. Probably true for other correlations withing the full sample as well.

pdf("EndothelialvsTimeOnTopVideo.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[2,]~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_on_top_of_stimulus_animal_Video, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Endothelial Index", xlab="SI Test: Time on Top (Video)")
dev.off()
#That is completely unconvincing - profound floor effect.

pdf("EndothelialvsTimeApproaching.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[2,]~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_approaching_stimulus_animal_Video, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="Endothelial Index", xlab="SI Test: Time Approaching (Video)")
dev.off()
#not very convincing.

pdf("RBCvsTimeOnTopVideo.pdf", width=5, height=5)
plot(CellTypeOutput$AveragePrimary_CellTypeIndex[10,]~Sample_BehavHormonalData_NACC_NoOutliers_Ordered$time_on_top_of_stimulus_animal_Video, col=as.factor(Sample_MetaData_NACC_NoOutliers_Ordered$date_of_dissection), ylab="RBC Index", xlab="SI Test: Time on Top (Video)")
dev.off()
#That is completely unconvincing - profound floor effect.


############################################################

#Functional Ontology:

#It looks like GSEA has updated their software and version of code available for R - I can't remember why we originally switched over to using fGSEA vs. GSEA - I think GSEA was being slow and buggy. I'm going to try out the new software and see if it works.

library("devtools")
install_github("GSEA-MSigDB/GSEA_R")
# ERROR: this R is version 3.4.1, package 'GSEA' requires R >= 3.6.0
# Installation failed: Command failed (1)
# Warning message:
#   In strptime(x, fmt, tz = "GMT") :
#   unknown timezone 'zone/tz/2021a.1.0/zoneinfo/America/Detroit'

#Sigh. That is a whole new project to tackle - updating everything. 

#I downloaded the app and then remembered our issue with it before - the tool tries to perform its own DE analysis without any concern for co-variates (i.e., you can only load in expression data and it runs the most basic model). Total crap.
#There is a version called "pre-ranked" that uses a pre-ranked list (e.g., our own DE genes list).  I think this was true before, and I can't remember why we didn't end up using it
#Well, to begin with, it only goes up to 1000 permutations (ouch!)
#And it is *crazy* slow, especially for such a pathetic number of permutations. I've had a single analysis running for 15 min now - I'm not sure if it has just crapped out on me.
#Also sounds like the algorithm only uses rank in the calculations (that may be true for fgsea too, I don't remember) - which seems like it could create some weird artifacts if there is a strong positive or negative skew in the log2FC results.

#Instead, I just downloaded the new MsigDb database (v7.3, downloaded 2021-03-25), curated it, added in our brain cell type signatures, stress signatures, NACC and HC cell type and co-expression clusters. 
#Then converted our results to human orthologs.
#See other code document for details: "Code_ImprovingGMTforSTRHC.R"


#And then ran fGSEA.

library(fgsea)

#I quickly upated the genes with date-related names offline to match the new (better!) nomenclature (e.g., March1->Marchf1, Sept1->Septin1)

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")
NACC_EESD_Combined<-read.csv("NACC_EESD_Combined_DateGenesFixed.csv", header=TRUE, stringsAsFactors = FALSE)
dim(NACC_EESD_Combined)
#[1] 17765    28

colnames(NACC_EESD_Combined)

sum(is.na(NACC_EESD_Combined$gene_symbol))
#[1] 766

NACC_EESD_Combined_noNA<-NACC_EESD_Combined[is.na(NACC_EESD_Combined$gene_symbol)==FALSE,]
dim(NACC_EESD_Combined_noNA)
#[1] 16999    28

sum(duplicated(NACC_EESD_Combined_noNA$gene_symbol))
#[1] 143

NAcc_EE_Betas_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$Coef.EnrichmentEE_MainEffectsModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_EE_Betas_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

length(NAcc_EE_Betas_forGSEA)
#[1] 16856

NAcc_EE_Betas_forGSEARanked<-NAcc_EE_Betas_forGSEA[order(NAcc_EE_Betas_forGSEA)]
head(NAcc_EE_Betas_forGSEARanked)
# Scn11a AABR07028009.1            Ak7     AC128059.4   LOC103693430        Tmem212 
# -2.236         -2.169         -2.026         -1.862         -1.845         -1.836 


NAcc_SD_Betas_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$Coef.SocialDefeatSD_MainEffectsModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_SD_Betas_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

NAcc_SD_Betas_forGSEARanked<-NAcc_SD_Betas_forGSEA[order(NAcc_SD_Betas_forGSEA)]
head(NAcc_SD_Betas_forGSEARanked)
# AABR07044362.6 AABR07044362.7 AABR07044362.1     AC120486.2   LOC103690108           Olah 
# -3.376         -2.786         -2.044         -1.586         -1.288         -1.214


NAcc_EEinteract_Betas_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$Coef.EnrichmentEE_InteractionModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_EEinteract_Betas_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

NAcc_EEinteract_Betas_forGSEARanked<-NAcc_EEinteract_Betas_forGSEA[order(NAcc_EEinteract_Betas_forGSEA)]
head(NAcc_EEinteract_Betas_forGSEARanked)
# AABR07028009.1   LOC103693430         Atp12a            Ak7        Tmem212     AC128059.4 
# -2.472         -2.425         -2.257         -2.155         -2.125         -2.093 

NAcc_SDinteract_Betas_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$Coef.SocialDefeatSD_InteractionModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_SDinteract_Betas_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

NAcc_SDinteract_Betas_forGSEARanked<-NAcc_SDinteract_Betas_forGSEA[order(NAcc_SDinteract_Betas_forGSEA)]
head(NAcc_SDinteract_Betas_forGSEARanked)
# AABR07044362.6 AABR07044362.7     AC120486.2 AABR07044362.1           Rln1 AABR07037356.1 
# -4.113         -3.568         -2.666         -2.426         -2.105         -1.980 


NAcc_EEbySD_Betas_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$Coef.EnrichmentEE.SocialDefeatSD_InteractionModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_EEbySD_Betas_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

NAcc_EEbySD_Betas_forGSEARanked<-NAcc_EEbySD_Betas_forGSEA[order(NAcc_EEbySD_Betas_forGSEA)]
head(NAcc_EEbySD_Betas_forGSEARanked)
# Slc25a54 AABR07024869.1      LOC500567 AABR07058934.1           Tbx4     AC139608.1 
# -2.795         -2.604         -1.912         -1.897         -1.871         -1.769 


#Let's start out with a basic Gene Ontology analysis (C5)

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("c5.all.v7.3.symbols_RatOrtholog.gmt.txt")
str(GMT_ForRats)
# List of 14996
# $ GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE                                                                                                                                                       : chr [1:2081] "Akt3" "Ppargc1a" "Polg2" "Parp1" ...
# $ GOBP_REPRODUCTION                                                                                                                                                                           : chr [1:2081] "Ada" "Gnpda1" "Zglp1" "Syce1l" ...
# $ GOBP_SINGLE_STRAND_BREAK_REPAIR                                                                                                                                                             : chr [1:2081] "Ercc8" "Parp1" "Aplf" "Ercc6" ...

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/GeneOntology_C5")

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_MainEffectsModel_FGSEAResults_BasicC5gmt.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SD_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_MainEffectsModel_FGSEAResults_BasicC5gmt.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_EEinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_InteractionModel_FGSEAResults_BasicC5gmt.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SDinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_InteractionModel_FGSEAResults_BasicC5gmt.csv")

rm(temp1)

#Results seem more consistent this time around. *Phew*

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("BrainCellTypesAndFunction_CustomGeneSets_V1_Rat_GMT.gmt.txt")
str(GMT_ForRats)
# List of 549
# $ Astrocyte_All_Cahoy_JNeuro_2008                                           : chr [1:2082] "Gfap" "Aqp4" "Pla2g7" "Slc39a12" ...
# $ Astrocyte_All_Darmanis_PNAS_2015                                          : chr [1:2082] "Fgfr3" "Aqp4" "Gja1" "Agt" ...
# $ Astrocyte_All_Doyle_Cell_2008                                             : chr [1:2082] "Adora2b" "Cyp4f4" "Fzd2" "Nr2e1" ...

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/BrainCellTypeFunction")

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_MainEffectsModel_FGSEAResults_BrainCellTypeFunction.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SD_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_MainEffectsModel_FGSEAResults_BrainCellTypeFunction.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_EEinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_InteractionModel_FGSEAResults_BrainCellTypeFunction.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SDinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_InteractionModel_FGSEAResults_BrainCellTypeFunction.csv")

rm(temp1)



setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("c7.all.v7.3.symbols_RatOrtholog.gmt.txt")
str(GMT_ForRats)
# List of 5219
# $ KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP                                                                                                                     : chr [1:2091] "Rflnb" "Ampd3" "Nsg2" "Dusp6" ...
# $ KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_DN                                                                                                                     : chr [1:2091] "Hopx" "Id2" "Lgals1" "Haspin" ...
# $ KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_UP                                                                                                                    : chr [1:2091] "Eml5" "Klk8" "Galk1" "Mettl9" ...
# $ KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_DN                                                                                                                    : chr [1:2091] "Bcl2a1" "Id2" "Cox17" "Klrg1" ...
# $ KAECH_NAIVE_VS_MEMORY_CD8_TCELL_UP                                                                                                                       : chr [1:2091] "Mfhas1" "Prkcb" "Ramp1" "Retreg1" ...
# $ KAECH_NAIVE_VS_MEMORY_CD8_TCELL_DN                                                                                                                       : chr [1:2091] "Tspan31" "Ccr5" "Ifng" "Faslg" ...

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/Immune_C7")

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_MainEffectsModel_FGSEAResults_ImmuneC7.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SD_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_MainEffectsModel_FGSEAResults_ImmuneC7.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_EEinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_InteractionModel_FGSEAResults_ImmuneC7.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SDinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_InteractionModel_FGSEAResults_ImmuneC7.csv")

rm(temp1)

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("c3.all.v7.3.symbols_RatOrtholog.gmt.txt")
str(GMT_ForRats)
# List of 3731
# $ MIR153_5P                                                                                                  : chr [1:2099] "Mest" "Nudt4" "Agfg1" "Ipo5" ...
# $ MIR8485                                                                                                    : chr [1:2099] "Mboat2" "Ipo5" "Spata17" "Creb5" ...
# $ MIR3662                                                                                                    : chr [1:2099] "Vsig1" "Nudt4" "Agfg1" "Mbtps2" ...
# $ MIR607                                                                                                     : chr [1:2099] "Esrrg" "Pias2" "Agfg1" "Mbtps2" ...

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/Regulatory_C3")


temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_MainEffectsModel_FGSEAResults_Regulatory_C3.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SD_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_MainEffectsModel_FGSEAResults_Regulatory_C3.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_EEinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_InteractionModel_FGSEAResults_Regulatory_C3.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SDinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_InteractionModel_FGSEAResults_Regulatory_C3.csv")

rm(temp1)


setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("c2.all.v7.3.symbols_RatOrtholog.gmt.txt")
str(GMT_ForRats)
# List of 6255
# $ BENITEZ_GBM_PROTEASOME_INHIBITION_RESPONSE                                                                                                   : chr [1:2040] "Sox4" "Sesn2" "Sgta" "Sox11" ...
# $ BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_UP                                                                                               : chr [1:2040] "Cenpe" "Aspm" "Esm1" "Kif20b" ...
# $ BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_DN                                                                                               : chr [1:2040] "LOC683573" "Rgcc" "Fgfr4" "Bmp4" ...
# $ BLANCO_MELO_MERS_COV_INFECTION_MCR5_CELLS_UP                                                                                                 : chr [1:2040] "Hspa1a" "Hspa1b" "Dnaja4" "Lcn10" ...

#This GMT contains a lot of random crap (esp cancer), but it also has all of the reactome and kegg pathways - it might be worth it to whidle it down.
#Esp. since I already pulled out the brain-related custom gene sets.


setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/Curated_C2")


temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_MainEffectsModel_FGSEAResults_Curated_C2.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SD_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_MainEffectsModel_FGSEAResults_Curated_C2.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_EEinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_InteractionModel_FGSEAResults_Curated_C2.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SDinteract_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_InteractionModel_FGSEAResults_Curated_C2.csv")

rm(temp1)


#One more: Let's combine the traditional C5 (ontology) with the custom brain cell types and functions GMT so that we get an analysis that corrects for FDR across both.

library(fgsea)

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("c5withBrainCellTypesFunction.all.v7.3.symbols_RatOrtholog.gmt.txt")
str(GMT_ForRats)
# List of 15545
# $ GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE                                                                                                                                                       : chr [1:2082] "Akt3" "Ppargc1a" "Polg2" "Parp1" ...
# $ GOBP_REPRODUCTION                                                                                                                                                                           : chr [1:2082] "Ada" "Gnpda1" "Zglp1" "Syce1l" ...
# $ GOBP_SINGLE_STRAND_BREAK_REPAIR                                                                                                                                                             : chr [1:2082] "Ercc8" "Parp1" "Aplf" "Ercc6" ...

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/GeneOntologyC5_wBrainCellTypeFunction")

#Note - I also increased the # permutations for this, because the output from the others was looking like there might not be enough.

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_MainEffectsModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SD_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_MainEffectsModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_EEinteract_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_InteractionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SDinteract_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_InteractionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction.csv")

rm(temp1)

#There is a pretty clear negative correlation between the extremeness of the log2FC and A (average log2 expression) that I think might be biasing the results towards pathways with a large number of low-expressed genes.
#I'm going to try running a version using t-statistics to rank genes instead of log2FC:

NAcc_EE_Tstat_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$t.EnrichmentEE_MainEffectsModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_EE_Tstat_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

length(NAcc_EE_Tstat_forGSEA)
#[1] 16856

NAcc_EE_Tstat_forGSEARanked<-NAcc_EE_Tstat_forGSEA[order(NAcc_EE_Tstat_forGSEA)]
head(NAcc_EE_Tstat_forGSEARanked)
# Pcdhb6        Pcdhga2         Scn11a          Dhrs2         Pcdhb5 AABR07027137.1 
# -5.63          -5.55          -5.29          -5.12          -5.04          -4.74 


NAcc_SD_Tstat_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$t.SocialDefeatSD_MainEffectsModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_SD_Tstat_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

NAcc_SD_Tstat_forGSEARanked<-NAcc_SD_Tstat_forGSEA[order(NAcc_SD_Tstat_forGSEA)]
head(NAcc_SD_Tstat_forGSEARanked)
# Tmtc1        Reps2 LOC103690108        Mertk         Xkr4       Cyp4f4 
# -4.80        -4.24        -4.21        -4.07        -3.84        -3.80 

NAcc_EEinteract_Tstat_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$t.EnrichmentEE_InteractionModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_EEinteract_Tstat_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

NAcc_EEinteract_Tstat_forGSEARanked<-NAcc_EEinteract_Tstat_forGSEA[order(NAcc_EEinteract_Tstat_forGSEA)]
head(NAcc_EEinteract_Tstat_forGSEARanked)
# Pcdhb5    Pcdhga2     Pcdhb6      Dele1 AC120486.2    Pcdhgb6 
# -6.41      -6.20      -5.96      -5.62      -4.69      -4.36 

NAcc_SDinteract_Tstat_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$t.SocialDefeatSD_InteractionModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_SDinteract_Tstat_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

NAcc_SDinteract_Tstat_forGSEARanked<-NAcc_SDinteract_Tstat_forGSEA[order(NAcc_SDinteract_Tstat_forGSEA)]
head(NAcc_SDinteract_Tstat_forGSEARanked)
# AC120486.2     Frmpd1      Tmtc1     Nt5c1a     Pcdh17     Atp8a2 
# -5.60      -5.00      -4.37      -4.30      -4.25      -4.17 

NAcc_EEbySD_Tstat_forGSEA<-tapply(X=NACC_EESD_Combined_noNA$t.EnrichmentEE.SocialDefeatSD_InteractionModel, INDEX=NACC_EESD_Combined_noNA$gene_symbol, FUN=mean)
names(NAcc_EEbySD_Tstat_forGSEA)<-names(table(NACC_EESD_Combined_noNA$gene_symbol))

NAcc_EEbySD_Tstat_forGSEARanked<-NAcc_EEbySD_Tstat_forGSEA[order(NAcc_EEbySD_Tstat_forGSEA)]
head(NAcc_EEbySD_Tstat_forGSEARanked)
# Polr3g AABR07024869.1         Abca12         Rab33b       Zmpste24           Wwc2 
# -4.92          -4.60          -4.26          -4.13          -4.10          -3.91 

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/GeneOntologyC5_wBrainCellTypeFunction")

temp1<-fgsea(GMT_ForRats, NAcc_EE_Tstat_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_MainEffectsModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction_Tstat.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SD_Tstat_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_MainEffectsModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction_Tstat.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_EEinteract_Tstat_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_EE_InteractionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction_Tstat.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, NAcc_SDinteract_Tstat_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "NAcc_SD_InteractionModel_FGSEAResults_OntologyC5_wBrainCellTypeFunction_Tstat.csv")

rm(temp1)

#That didn't change the significance of the ependymal-related results, but did seem to expand the results to include other kinds of pathways.
#Let's take a look at the findings tomorrow - maybe it will clean up some of the HC results too.
#Ooh - even more pathways related to autism and social behavior highlighted now (???!!!)


#Let's try the collapse pathways function to see if we can simplify the results into something more digestible:

temp1<-fgsea(GMT_ForRats, NAcc_EE_Tstat_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)
str(temp1)
#Classes ‘data.table’ and 'data.frame':	10224 obs. of  8 variables:

temp1_JustFDR05<-temp1[temp1$padj<0.05,]
str(temp1_JustFDR05)
#Classes ‘data.table’ and 'data.frame':	489 obs. of  8 variables:

temp1_Collapsed<-collapsePathways(fgseaRes=temp1_JustFDR05, pathways=temp1_JustFDR05$pathway, stats=NAcc_EE_Tstat_forGSEARanked)
#This isn't working - it turns out that is because this function is only in newer version of FGSEA. I should probably update... nervous that I'm going to lose functionality in some way. I may run my other necessary analyses and come back to this...
rm(temp1)


##################

#Digging into the NACC social behavior results:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/GeneOntologyC5_wBrainCellTypeFunction")

LeadingEdge_SDEE_SocialBehaviorGeneSets<-read.csv("LeadingEdge_SDEE_SocialBehaviorGeneSets.csv", header=TRUE, stringsAsFactors = FALSE)
str(LeadingEdge_SDEE_SocialBehaviorGeneSets)

# 'data.frame':	181 obs. of  10 variables:
# $ SD_Gandal_2018_AutismSpectrumDisorder_Downregulated_Cortex: chr  "Tcf25" "Uba1" "Nap1l5" "Cpsf3" ...
# $ SD_Gandal_2018_AutismSpectrumDisorder_Upregulated_Cortex  : chr  "RT1-CE4" "RT1-CE5" "RT1-Da" "Cd44" ...
# $ SD_HP_AUTISTIC_BEHAVIOR                                   : chr  "Dcps" "Rora" "Dmxl2" "Scn8a" ...
# $ SD_HP_AUTISM                                              : chr  "Clip2" "Usf3" "Tsc2" "Nhs" ...
# $ SD_HP_IMPAIRED_SOCIAL_INTERACTIONS                        : chr  "Cox10" "Ireb2" "Alg11" "Uchl1" ...
# $ SD_HP_ABNORMAL_SOCIAL_BEHAVIOR                            : chr  "Ireb2" "Abca7" "Alg11" "Clip2" ...
# $ SD_HP_ABNORMAL_AGGRESSIVE_IMPULSIVE_OR_VIOLENT_BEHAVIOR   : chr  "Mab21l1" "Nfasc" "Brsk2" "Magel2" ...
# $ SD_HP_AGGRESSIVE_BEHAVIOR                                 : chr  "Mapk1" "Lingo1" "Cux2" "Pah" ...
# $ EE_Gandal_2018_AutismSpectrumDisorder_Downregulated_Cortex: chr  "Sncg" "Panx2" "Bles03" "Acbd6" ...
# $ EE_Gandal_2018_AutismSpectrumDisorder_Upregulated_Cortex  : chr  "Sh3glb1" "Tm4sf1" "Cmtm3" "Trip6" ...

#How much do they overlap with each other?

SimilarityBetweenSocialLeadingEdgeGenes<-matrix(0,10,10)
row.names(SimilarityBetweenSocialLeadingEdgeGenes)<-colnames(LeadingEdge_SDEE_SocialBehaviorGeneSets)
colnames(SimilarityBetweenSocialLeadingEdgeGenes)<-colnames(LeadingEdge_SDEE_SocialBehaviorGeneSets)

Length_SocialLeadingEdgeGenes<-matrix(0,10,1)
row.names(Length_SocialLeadingEdgeGenes)<-colnames(LeadingEdge_SDEE_SocialBehaviorGeneSets)

for(i in c(1:10)){
  GeneSetI<-LeadingEdge_SDEE_SocialBehaviorGeneSets[LeadingEdge_SDEE_SocialBehaviorGeneSets[,i]!="",i]
  Length_SocialLeadingEdgeGenes[i,1]<-length(GeneSetI)
  for(j in c(1:10)){
    GeneSetJ<-LeadingEdge_SDEE_SocialBehaviorGeneSets[LeadingEdge_SDEE_SocialBehaviorGeneSets[,j]!="",j]
    SimilarityBetweenSocialLeadingEdgeGenes[i,j]<-sum(GeneSetI%in%GeneSetJ)
  }
}

write.csv(Length_SocialLeadingEdgeGenes, "Length_SocialLeadingEdgeGenes.csv")
write.csv(SimilarityBetweenSocialLeadingEdgeGenes, "SimilarityBetweenSocialLeadingEdgeGenes.csv")


SimilarityBetweenSocialLeadingEdgeGenes_Proportion<-SimilarityBetweenSocialLeadingEdgeGenes
for(j in c(1:10)){
  SimilarityBetweenSocialLeadingEdgeGenes_Proportion[,j]<-SimilarityBetweenSocialLeadingEdgeGenes[,j]/Length_SocialLeadingEdgeGenes
}

write.csv(SimilarityBetweenSocialLeadingEdgeGenes_Proportion, "SimilarityBetweenSocialLeadingEdgeGenes_Proportion.csv")


#All right, let's try to loop it over the full gene set database (not just leading genes):

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("c5withBrainCellTypesFunction.all.v7.3.symbols_RatOrtholog.gmt.txt")
str(GMT_ForRats)
# List of 15545
# $ GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE                                                                                                                                                       : chr [1:2082] "Akt3" "Ppargc1a" "Polg2" "Parp1" ...
names(GMT_ForRats)
GMT_ForRats[[1]]

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/GeneOntologyC5_wBrainCellTypeFunction")

SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets<-matrix(0,15545,10)
row.names(SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets)<-names(GMT_ForRats)
colnames(SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets)<-colnames(LeadingEdge_SDEE_SocialBehaviorGeneSets)

for(i in c(1:10)){
  GeneSetI<-LeadingEdge_SDEE_SocialBehaviorGeneSets[LeadingEdge_SDEE_SocialBehaviorGeneSets[,i]!="",i]
  for(j in c(1:15545)){
    GeneSetJ<-GMT_ForRats[[j]][GMT_ForRats[[j]]!=""]
    SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets[j,i]<-sum(GeneSetI%in%GeneSetJ)
  }
}

write.csv(SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets, "SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets.csv")

SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets_Proportion<-SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets
for(j in c(1:15545)){
  SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets_Proportion[j,]<-SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets[j,]/Length_SocialLeadingEdgeGenes
}

write.csv(SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets_Proportion, "SimilarityBetweenSocialLeadingEdgeGenes_AndAllGeneSets_Proportion.csv")


#Pulling out gene sets containing the top genes for EE (NAcc):

GeneSetsWithEEGenes<-lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb5", "Pcdhga2", "Pcdhb6", "Pcdhb8", "Scn11a", "Dhrs2", "AABR07027137.1", "Megf8", "Mbd1", "Csf3r", "AABR07066700.1", "Ccl3")))
str(GeneSetsWithEEGenes)
write.csv(simplify2array(GeneSetsWithEEGenes), "GeneSetsWithEEGenes.csv")

GeneSetsWithEEGenes_wInteractionModel<-lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb5", "Pcdhga2", "Pcdhb6", "Pcdhb8", "Scn11a", "Dhrs2", "AABR07027137.1", "Megf8", "Mbd1", "Csf3r", "AABR07066700.1", "Ccl3", "Dele1", "Pcdhb7", "AC120486.2")))
str(GeneSetsWithEEGenes_wInteractionModel)
write.csv(simplify2array(GeneSetsWithEEGenes_wInteractionModel), "GeneSetsWithEEGenes_wInteractionModel.csv")

EE_Pcdhb5<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb5"))))
EE_Pcdhga2<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhga2"))))
EE_Pcdhb6<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb6"))))
EE_Pcdhb8<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb8"))))
EE_Scn11a<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Scn11a"))))
EE_Dhrs2<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Dhrs2"))))
EE_AABR07027137<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("AABR07027137.1"))))
EE_Megf8<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Megf8"))))
EE_Mbd1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Mbd1"))))
EE_Csf3r<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Csf3r"))))
EE_AABR07066700<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("AABR07066700.1"))))
EE_Ccl3<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Ccl3"))))
EE_Dele1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Dele1"))))
EE_Pcdhb7<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Pcdhb7"))))
EE_AC120486<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("AC120486.2"))))

head(data.frame(EE_Pcdhb5,EE_Pcdhga2,EE_Pcdhb6,EE_Pcdhb8,EE_Scn11a,EE_Dhrs2,EE_AABR07027137,EE_Megf8,EE_Mbd1,EE_Csf3r,EE_AABR07066700,EE_Ccl3, EE_Dele1,EE_Pcdhb7,EE_AC120486))

write.csv(data.frame(EE_Pcdhb5,EE_Pcdhga2,EE_Pcdhb6,EE_Pcdhb8,EE_Scn11a,EE_Dhrs2,EE_AABR07027137,EE_Megf8,EE_Mbd1,EE_Csf3r,EE_AABR07066700,EE_Ccl3, EE_Dele1,EE_Pcdhb7,EE_AC120486), "GeneSetsWithEEGenes_wInteractionModel_ByGene.csv")

#Pulling out gene sets containing the top genes for SD (NAcc):

GeneSetsWithSDGenes<-lapply(GMT_ForRats, function(y) sum(y%in%c("RT1-N2","RT1-CE4","RT1-CE5", "Slc26a8", "Tmtc1")))
str(GeneSetsWithSDGenes)
write.csv(simplify2array(GeneSetsWithSDGenes), "GeneSetsWithSDGenes.csv")

GeneSetsWithSDGenes_wInteractionModel<-lapply(GMT_ForRats, function(y) sum(y%in%c("RT1-N2","RT1-CE4","RT1-CE5", "Slc26a8", "Tmtc1", "AC120486.2", "Abca12", "Frmpd1", "LOC360919", "Polr3g")))
str(GeneSetsWithSDGenes_wInteractionModel)
write.csv(simplify2array(GeneSetsWithSDGenes_wInteractionModel), "GeneSetsWithSDGenes_wInteractionModel.csv")

SD_RT1N2<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("RT1-N2"))))
SD_RT1CE4<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("RT1-CE4"))))
SD_RT1CE5<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("RT1-CE5"))))
SD_Slc26a8<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Slc26a8"))))
SD_Tmtc1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Tmtc1"))))
SD_AC120486<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("AC120486.2"))))
SD_Abca12<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Abca12"))))
SD_Frmpd1<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Frmpd1"))))
SD_LOC360919<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("LOC360919"))))
SD_Polr3g<-simplify2array(lapply(GMT_ForRats, function(y) sum(y%in%c("Polr3g"))))

write.csv(data.frame(SD_RT1N2,SD_RT1CE4,SD_RT1CE5,SD_Slc26a8,SD_Tmtc1,SD_AC120486,SD_Abca12,SD_Frmpd1,SD_LOC360919,SD_Polr3g),"GeneSetsWithSDGenes_wInteractionModel_ByGene.csv")


##################################

#figure out which of the top DE results are also DE for HR/LR:
setwd("~/Documents/Microarray Gen/HRLR/ThesisMetaAnalysisOutput")

MetaAnalysesResults_AllAges<-read.csv("MetaAnalysesResults_AllAges.csv", header=T, stringsAsFactors = F)

str(MetaAnalysesResults_AllAges)
#Stupid date-gene problem. I'm going to ignore it for now - for a final version of these results we're going to have to fix that. Again.
colnames(MetaAnalysesResults_AllAges)[2]<-"gene_symbol"

colnames(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect)[2]<-"gene_symbol"

NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis<-join(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect, MetaAnalysesResults_AllAges, by="gene_symbol", type="left")

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers")

write.csv(NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis, "NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis.csv")

#Pcdhga1 is a HRvsLR DE gene (Adult FDR<0.01) that also has strong nominal effects for SD and EE
#Etv4 is a HRvsLR DE gene (Adult FDR<0.05) with nominal effects for both SD and EE
#Uhrf1 is a HRvsLR DE gene (Adult FDR<0.05) with nominal effects for SD and EE
#Tmem144 is an HRvsLR DE gene (Adult FDR<0.05 - also DE in many other models) with sig effect of SD and nominal effect of EE

#The top EE gene, Pcdhb5, has nominal effects for adult HR/LR (p=0.01) in the same direction, and similar direction of effect at P14.
#Same for Dhrs2 & Dhrs1, Fos, Sirt4, Synj2, Myh6, Foxred1, AABR07054319.1, Arfgap1, Pygm, Faxc, Rspo1, Zfp180, Tspan1, Adam23, Odf3, Unc45a, Car9, Hist1h2bk, Pex11a...
#Ifnlr1, Grhl3, Vars2, Nhlrc3, Barx2, Egfr, Slc13a5, Frmpd1, Atp8a2, Synj2, Nhlrc3, Cyp4f4, Mgarp, Vars2 show nominal effects for both SD and EE, and for adult HR/LR
#And even more overlap for just SD...

#This seems like a lot. Maybe we should do some sort of enrichment analysis/Venn diagram. (note - this was later done as part of fgsea,  now included above. There used to be code for a less formal/proper analysis too, but I cut it.)
#Also: I would like to know how often the effects are a reversal of HR/LR differences

#How well do the results from both variables correlate with the HR/LR direction of effect?

pdf("Plot_EEvsHRLRcoefficients_ModelwInteractionCovariates.pdf", width=5, height=5)
plot(EnrichmentEE~Adult_AdultMeta.estimate, data=NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis$Adult_AdultMeta.rawp<0.05 & NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis$EnrichmentEE.2<0.05,],  xlab="bHR/bLR Meta-Analysis Coefficients (genes with p<0.05)", ylab="EE Coefficients (genes with p<0.05)")

BestFitLine<-lm(EnrichmentEE~Adult_AdultMeta.estimate, data=NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis$Adult_AdultMeta.rawp<0.05 & NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis$EnrichmentEE.2<0.05,])
abline(BestFitLine)
dev.off()
#Definitely no correlation there...

pdf("Plot_SDvsHRLRcoefficients_ModelwInteractionCovariates.pdf", width=5, height=5)
plot(SocialDefeatSD~Adult_AdultMeta.estimate, data=NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis$Adult_AdultMeta.rawp<0.05 & NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis$SocialDefeatSD.2<0.05,],  xlab="bHR/bLR Meta-Analysis Coefficients (genes with p<0.05)", ylab="SD Coefficients (genes with p<0.05)")

BestFitLine<-lm(SocialDefeatSD~Adult_AdultMeta.estimate, data=NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis[NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis$Adult_AdultMeta.rawp<0.05 & NAcc_LimmaResults_EEbySD_RNAConcExtractDissect_vsMetaAnalysis$SocialDefeatSD.2<0.05,])
abline(BestFitLine)
dev.off()
#slight positive correlation, which is kind-of the opposite of what we would expect (because Negative means up in LRs)
summary.lm(BestFitLine)
#but very much not sig, so basically no correlation.

######################################

#Let's try performing the same analysis using the model output without an interaction term, since I am a little suspicious about the correlation that we are seeing between the EE and SD results.

colnames(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect)[2]<-"gene_symbol"

NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis<-join(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect, MetaAnalysesResults_AllAges, by="gene_symbol", type="left")

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers")

write.csv(NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis, "NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis.csv")

#Visually it seems like there is less overlap than before. 

#How well do the results from both variables correlate with the HR/LR direction of effect?

pdf("Plot_EEvsHRLRcoefficients_ModelNoInteractionCovariates.pdf", width=5, height=5)
plot(EnrichmentEE~Adult_AdultMeta.estimate, data=NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis$Adult_AdultMeta.rawp<0.05 & NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis$EnrichmentEE.2<0.05,],  xlab="bHR/bLR Meta-Analysis Coefficients (genes with p<0.05)", ylab="EE Coefficients (genes with p<0.05)")

BestFitLine<-lm(EnrichmentEE~Adult_AdultMeta.estimate, data=NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis$Adult_AdultMeta.rawp<0.05 & NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis$EnrichmentEE.2<0.05,])
abline(BestFitLine)
dev.off()
#Definitely no correlation there...

pdf("Plot_SDvsHRLRcoefficients_ModelNoInteractionCovariates.pdf", width=5, height=5)
plot(SocialDefeatSD~Adult_AdultMeta.estimate, data=NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis$Adult_AdultMeta.rawp<0.05 & NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis$SocialDefeatSD.2<0.05,],  xlab="bHR/bLR Meta-Analysis Coefficients (genes with p<0.05)", ylab="SD Coefficients (genes with p<0.05)")

BestFitLine<-lm(SocialDefeatSD~Adult_AdultMeta.estimate, data=NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis[NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis$Adult_AdultMeta.rawp<0.05 & NAcc_LimmaResults_onlyEESD_RNAConcExtractDissect_vsMetaAnalysis$SocialDefeatSD.2<0.05,])
abline(BestFitLine)
dev.off()
#slight positive correlation, which is kind-of the opposite of what we would expect (because Negative means up in LRs)
summary.lm(BestFitLine)
#but very much not sig, so basically no correlation.


################################################
