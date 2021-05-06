#Making a better .gmt file for analyzing our HC and STR results.
#Megan Hagenauer
#2021-03-26

################

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")


###############################################

#Reading in DropViz cell type markers:

#Dropviz.org accessed 03-25-2021
#Query->Parameters->Clusters->CompareClusters
#Ran two queries:
#All Target Clusters within the striatum (STR) vs. as comparison "Rest of Striatum"
#All Target Clusters within the hippocampus (HC) vs. as comparison "Rest of Hippocampus"
#Differential Expression Criteria: 
#Minimum Fold Ratio of 4 (log2 = difference of 2)
#Maximum P-value exponent: -30
#Min mean Log amount in Target: 0.5 (1.4 transcripts per cell)
#Max Mean Log amount in Comparison group: 6 (highest possible choice) - so basically we don't care if the gene is highly expressed, just that it is more highly expressed in the cell type of interest.

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/DropViz_HC")

DropViz_HC_Files<-list.files()

#To get started:
DropViz_HC_CellTypeClusterMarkers<-read.csv(DropViz_HC_Files[1], header=TRUE, stringsAsFactors = FALSE)
str(DropViz_HC_CellTypeClusterMarkers)
DropViz_HC_CellTypeClusterMarkers$CellType<-rep(paste(DropViz_HC_Files[1]), nrow(DropViz_HC_CellTypeClusterMarkers))

DropViz_HC_CellTypeClusterMarkers_Concise<-DropViz_HC_CellTypeClusterMarkers[DropViz_HC_CellTypeClusterMarkers$fc.disp>10,]

for(i in c(2:length(DropViz_HC_Files))){
Temp<-read.csv(DropViz_HC_Files[i], header=TRUE, stringsAsFactors = FALSE)
Temp$CellType<-rep(paste(DropViz_HC_Files[i]), nrow(Temp))

print(paste(DropViz_HC_Files[i]))
print(dim(Temp))
DropViz_HC_CellTypeClusterMarkers<-rbind(DropViz_HC_CellTypeClusterMarkers, Temp)

if(nrow(Temp)>50){
  print(sum(Temp$fc.disp>10))
  if(sum(Temp$fc.disp>10)>50){
    Temp_Concise<-Temp[Temp$fc.disp>10,]
  }else{Temp_Concise<-Temp[c(1:50),]}
  print(dim(Temp_Concise))
  DropViz_HC_CellTypeClusterMarkers_Concise<-rbind(DropViz_HC_CellTypeClusterMarkers_Concise, Temp_Concise)
  rm(Temp_Concise)
}
rm(Temp)
}

str(DropViz_HC_CellTypeClusterMarkers)
#'data.frame':	4780 obs. of  21 variables
str(DropViz_HC_CellTypeClusterMarkers_Concise)
#'data.frame':	1699 obs. of  21 variables

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")
write.csv(DropViz_HC_CellTypeClusterMarkers, "DropViz_HC_CellTypeClusterMarkers.csv")
write.csv(DropViz_HC_CellTypeClusterMarkers_Concise, "DropViz_HC_CellTypeClusterMarkers_Concise.csv")

##############

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/DropViz_STR")

DropViz_STR_Files<-list.files()

#To get started:
DropViz_STR_CellTypeClusterMarkers<-read.csv(DropViz_STR_Files[1], header=TRUE, stringsAsFactors = FALSE)
str(DropViz_STR_CellTypeClusterMarkers)
DropViz_STR_CellTypeClusterMarkers$CellType<-rep(paste(DropViz_STR_Files[1]), nrow(DropViz_STR_CellTypeClusterMarkers))

DropViz_STR_CellTypeClusterMarkers_Concise<-DropViz_STR_CellTypeClusterMarkers[DropViz_STR_CellTypeClusterMarkers$fc.disp>10,]

for(i in c(2:length(DropViz_STR_Files))){
  Temp<-read.csv(DropViz_STR_Files[i], header=TRUE, stringsAsFactors = FALSE)
  Temp$CellType<-rep(paste(DropViz_STR_Files[i]), nrow(Temp))
  
  print(paste(DropViz_STR_Files[i]))
  print(dim(Temp))
  DropViz_STR_CellTypeClusterMarkers<-rbind(DropViz_STR_CellTypeClusterMarkers, Temp)
  
  if(nrow(Temp)>50){
    print(sum(Temp$fc.disp>10))
    if(sum(Temp$fc.disp>10)>50){
      Temp_Concise<-Temp[Temp$fc.disp>10,]
    }else{Temp_Concise<-Temp[c(1:50),]}
    print(dim(Temp_Concise))
    DropViz_STR_CellTypeClusterMarkers_Concise<-rbind(DropViz_STR_CellTypeClusterMarkers_Concise, Temp_Concise)
    rm(Temp_Concise)
  }
  rm(Temp)
}

str(DropViz_STR_CellTypeClusterMarkers)
#'data.frame':	5021 obs. of  21 variables:
str(DropViz_STR_CellTypeClusterMarkers_Concise)
#'data.frame':	1769 obs. of  21 variables:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")
write.csv(DropViz_STR_CellTypeClusterMarkers, "DropViz_STR_CellTypeClusterMarkers.csv")
write.csv(DropViz_STR_CellTypeClusterMarkers_Concise, "DropViz_STR_CellTypeClusterMarkers_Concise.csv")

#################
#reading in some of the other files:

Zeisel_CA1PyramidalIndex<-read.csv("Zeisel CA1 Pyramidal Index.csv", header=TRUE, stringsAsFactors=FALSE)
str(Zeisel_CA1PyramidalIndex)
#'data.frame':	409 obs. of  1 variable:
#$ Zeisel.CA1.Pyramidal.Index: chr  "0610010B08RIK_LOC4" "1500015A07RIK" "2610018G03RIK" "4831440E17RIK" ...

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_toMakeRat")
GMTfiles_ToMakeRat<-list.files()


#Original read-in code - it turned out that this was actually problematic (see notes below):
# c2.all.v7.3.symbols<-read.delim("c2.all.v7.3.symbols.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# c3.all.v7.3.symbols<-read.delim("c3.all.v7.3.symbols.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# c5.all.v7.3.symbols<-read.delim("c5.all.v7.3.symbols.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# c5.go.v7.3.symbols<-read.delim("c5.go.v7.3.symbols.gmt.txt" , sep="\t", header=FALSE, stringsAsFactors = FALSE)
# c7.all.v7.3.symbols<-read.delim("c7.all.v7.3.symbols.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# c7.immunesigdb.v7.3.symbols<-read.delim("c7.immunesigdb.v7.3.symbols.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# c8.all.v7.3.symbols<-read.delim("c8.all.v7.3.symbols.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# c8.JustBrainBloodwOurs.v7.3.symbols<-read.delim("c8.JustBrainBloodwOurs.v7.3.symbols.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# h.all.v7.3.symbols<-read.delim("h.all.v7.3.symbols.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# msigdb.v7.3.symbols<-read.delim("msigdb.v7.3.symbols.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# CombinedHCspecific<-read.delim("TableS1_CombinedHCspecific.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)

#Alright, it turns out that the original read-in process wasn't working properly - it appeared to occassionally be misinterpreting parts of the GMT file for new line symbols (I think - I never did quite get to the heart of it - but my guess is that there are \n's mixed into the file text somewhere), which was causing there to artificially be more "pathways" (rows) than expected in the file, with the added pathways just being gene names (because they were actually part of the listed genes for the previous pathway).
#This function seemed to work to fix the issue:

gmtPathways_WithLink<-function (gmt.file) 
{
  pathwayLines <- strsplit(readLines(gmt.file), "\\t")
  pathways<-t(as.data.frame(stri_list2matrix(pathwayLines, byrow = FALSE), stringsAsFactors=FALSE))
}

c8.all.v7.3.symbol_Temp<-gmtPathways_WithLink("c8.all.v7.3.symbols.gmt.txt")
str(c8.all.v7.3.symbol_Temp)
#List of 673
# chr [1:673, 1:1772] "BUSSLINGER_ESOPHAGEAL_QUIESCENT_BASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_PROLIFERATING_BASAL_CELLS" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:673] "V1" "V2" "V3" "V4" ...
# ..$ : NULL
head(c8.all.v7.3.symbol_Temp)
c8.all.v7.3.symbol_Temp[c(1:10), c(1:10)]
#Got it.
rm(c8.all.v7.3.symbol_Temp)


c2.all.v7.3.symbols<-gmtPathways_WithLink("c2.all.v7.3.symbols.gmt.txt")
str(c2.all.v7.3.symbols)
# chr [1:6255, 1:1943] "BENITEZ_GBM_PROTEASOME_INHIBITION_RESPONSE" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_UP" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:6255] "V1" "V2" "V3" "V4" ...
# ..$ : NULL

c3.all.v7.3.symbols<-gmtPathways_WithLink("c3.all.v7.3.symbols.gmt.txt")
str(c3.all.v7.3.symbols)
# chr [1:3731, 1:2002] "MIR153_5P" "MIR8485" "MIR3662" "MIR607" "MIR616_5P" "MIR371B_5P" "MIR373_5P" "MIR6867_5P" "MIR12136" "MIR548AJ_3P_MIR548X_3P" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:3731] "V1" "V2" "V3" "V4" ...
# ..$ : NULL

c5.all.v7.3.symbols<-gmtPathways_WithLink("c5.all.v7.3.symbols.gmt.txt")
str(c5.all.v7.3.symbols)
# chr [1:14996, 1:1984] "GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE" "GOBP_REPRODUCTION" "GOBP_SINGLE_STRAND_BREAK_REPAIR" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14996] "V1" "V2" "V3" "V4" ...
# ..$ : NULL

c5.go.v7.3.symbols<-gmtPathways_WithLink("c5.go.v7.3.symbols.gmt.txt")
str(c5.go.v7.3.symbols)
# chr [1:10183, 1:1984] "GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE" "GOBP_REPRODUCTION" "GOBP_SINGLE_STRAND_BREAK_REPAIR" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:10183] "V1" "V2" "V3" "V4" ...
# ..$ : NULL

c7.all.v7.3.symbols<-gmtPathways_WithLink("c7.all.v7.3.symbols.gmt.txt")
str(c7.all.v7.3.symbols)
# chr [1:5219, 1:1994] "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP" "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_DN" "KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_UP" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:5219] "V1" "V2" "V3" "V4" ...
# ..$ : NULL

c7.immunesigdb.v7.3.symbols<-gmtPathways_WithLink("c7.immunesigdb.v7.3.symbols.gmt.txt")
str(c7.immunesigdb.v7.3.symbols)
# chr [1:4872, 1:202] "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP" "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_DN" "KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_UP" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:4872] "V1" "V2" "V3" "V4" ...
# ..$ : NULL

c8.all.v7.3.symbols<-gmtPathways_WithLink("c8.all.v7.3.symbols.gmt.txt")
str(c8.all.v7.3.symbols)
# chr [1:673, 1:1772] "BUSSLINGER_ESOPHAGEAL_QUIESCENT_BASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_PROLIFERATING_BASAL_CELLS" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:673] "V1" "V2" "V3" "V4" ...
# ..$ : NULL

c8.JustBrainBloodwOurs.v7.3.symbols<-gmtPathways_WithLink("c8.JustBrainBloodwOurs.v7.3.symbols.gmt.txt")
str(c8.JustBrainBloodwOurs.v7.3.symbols)
# chr [1:242, 1:2822] "FAN_EMBRYONIC_CTX_BIG_GROUPS_CAJAL_RETZIUS" "FAN_EMBRYONIC_CTX_BIG_GROUPS_BRAIN_ENDOTHELIAL" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:242] "V1" "V2" "V3" "V4" ...
# ..$ : NULL

h.all.v7.3.symbols<-gmtPathways_WithLink("h.all.v7.3.symbols.gmt.txt")
str(h.all.v7.3.symbols)
# chr [1:50, 1:202] "HALLMARK_TNFA_SIGNALING_VIA_NFKB" "HALLMARK_HYPOXIA" "HALLMARK_CHOLESTEROL_HOMEOSTASIS" "HALLMARK_MITOTIC_SPINDLE" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:50] "V1" "V2" "V3" "V4" ...
# ..$ : NULL

#I didn't end up re-doing the read-in process for this file because it didn't turn out to have the same problem (artificial line breaks) that the other MsigDb generated files had:
# CombinedHCspecific<-gmtPathways_WithLink("TableS1_CombinedHCspecific.gmt.txt")
# str(CombinedHCspecific)
#Original read in:
CombinedHCspecific<-read.delim("TableS1_CombinedHCspecific.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)



#While I'm at it, let's add in the gene sets for the top genes associated with chronic stress, internalizing behavior, and psychiatric illness:

#I already compiled some of these in another analysis:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/08_ComparisonWNestler")

HC_NACC_StressVsCtrl_Mice<-read.csv("HC_NACC_StressVsCtrl_Mice_Joined.csv", header=TRUE, stringsAsFactors = FALSE)
str(HC_NACC_StressVsCtrl_Mice)
#This document containts the genes identified as affected by chronic stress in Bagot et al. 2016 (NACC and HC), Bagot et al. 2017 (NACC and HC), and Pena et al. 2019 (just NACC)

#Hmmm... honestly, it may be easier to build a .gmt using the original suppl tables. Let me try that instead.
#Also might be useful to divide up the genes into upregulated and down-regulated effects.
#And maybe make the p-value cut-offs more stringent. P<0.05 nominal is pretty loosey goosey - probably contains a lot of noise.

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/StressStudies")

list.files()

Bagot_2016_SuscVsCTRL_NACC<-read.csv("Bagot_2016_SuscVsCTRL_NACC.csv", header=TRUE, stringsAsFactors = FALSE)
Bagot_2017_SuscVsCTRL_NACC<-read.csv("Bagot_2017_SuscVsCTRL_NACC.csv", header=TRUE, stringsAsFactors = FALSE)
Bagot_2016_SuscVsCTRL_HC<-read.csv("Bagot_2016_SuscVsCTRL_vHC.csv", header=TRUE, stringsAsFactors = FALSE)
Bagot_2017_SuscVsCTRL_HC<-read.csv("Bagot_2017_SuscVsCTRL_HC.csv", header=TRUE, stringsAsFactors = FALSE)
Pena_2019_SuscVsCTRL_NACC<-read.csv("Pena_2019_ELSSDvsCTRL_NACC.csv", header=TRUE, stringsAsFactors = FALSE)
Gray_2014_S2_CRS_HC<-read.csv("Gray_2014_S2_CRS_HC.csv", header=TRUE, stringsAsFactors = FALSE)
Gray_2014_S3_CRSwFST_HC<-read.csv("Gray_2014_S3_CRSwFST_HC.csv", header=TRUE, stringsAsFactors = FALSE)
Gray_2014_Suppl_FSTvsCtrl_HC<-read.csv("Gray_2014_Suppl_FSTvsCtrl_HC.csv", header=TRUE, stringsAsFactors = FALSE)
Gray_2014_S7_CortVsVeh_HC<-read.csv("Gray_2014_S7_CortVsVeh_HC.csv", header=TRUE, stringsAsFactors = FALSE)
Gray_2014_Suppl6_StressCortANOVA_HC<-read.csv("Gray_2014_Suppl6_StressCortANOVA_HC.csv", header=TRUE, stringsAsFactors = FALSE)

StressGeneSets<-data.frame(GeneSymbol_Mouse=character(), DatasetName=character(), stringsAsFactors = FALSE)
str(StressGeneSets)

str(Bagot_2016_SuscVsCTRL_NACC)
#'data.frame':	54 obs. of  12 variables:
sum(Bagot_2016_SuscVsCTRL_NACC$p_value<0.005)
#[1] 18
#Ouch.
sum(Bagot_2016_SuscVsCTRL_NACC$q_value<0.10)
#[1] 6
str(Bagot_2016_SuscVsCTRL_HC)
#'data.frame':	100 obs. of  12 variables:
sum(Bagot_2016_SuscVsCTRL_HC$p_value<0.005)
#[1] 34
sum(Bagot_2016_SuscVsCTRL_HC$q_value<0.10)
#[1] 14
#Ouch. So should I split it into up vs. down regulation or make the p-value more strict?  I don't think I can do both.

#Oh - also, to make proper gene sents, we need to remove duplicate symbols:

Bagot_2016_SuscVsCTRL_NACC_Genes<-Bagot_2016_SuscVsCTRL_NACC$gene[Bagot_2016_SuscVsCTRL_NACC$p_value<0.005]
sum(duplicated(Bagot_2016_SuscVsCTRL_NACC_Genes))
#[1] 0

Bagot_2016_SuscVsCTRL_HC_Genes<-Bagot_2016_SuscVsCTRL_HC$gene[Bagot_2016_SuscVsCTRL_HC$p_value<0.005]
sum(duplicated(Bagot_2016_SuscVsCTRL_HC_Genes))
#[1] 0

StressGeneSets<-rbind(StressGeneSets, data.frame(GeneSymbol_Mouse=Bagot_2016_SuscVsCTRL_NACC_Genes, DatasetName=rep("Bagot_2016_StressSusceptibleVsCTRL_NACC", length(Bagot_2016_SuscVsCTRL_NACC_Genes)), stringsAsFactors=FALSE), data.frame(GeneSymbol_Mouse=Bagot_2016_SuscVsCTRL_HC_Genes, DatasetName=rep("Bagot_2016_StressSusceptibleVsCTRL_HC",length(Bagot_2016_SuscVsCTRL_HC_Genes)), stringsAsFactors=FALSE))

str(StressGeneSets)


str(Bagot_2017_SuscVsCTRL_NACC)
#'data.frame':	314 obs. of  10 variables:
sum(Bagot_2017_SuscVsCTRL_NACC$P.Value<0.005)
#[1] 20
#Ouch
sum(Bagot_2017_SuscVsCTRL_NACC$adj.P.Val<0.10)
#[1] 0
str(Bagot_2017_SuscVsCTRL_HC)
#'data.frame':	164 obs. of  10 variables:
sum(Bagot_2017_SuscVsCTRL_HC$P.Value<0.005)
#[1] 16
sum(Bagot_2017_SuscVsCTRL_HC$adj.P.Val<0.10)
#[1] 0

Bagot_2017_SuscVsCTRL_NACC_Genes<-Bagot_2017_SuscVsCTRL_NACC$GeneID[Bagot_2017_SuscVsCTRL_NACC$P.Value<0.005]
sum(duplicated(Bagot_2017_SuscVsCTRL_NACC_Genes))
#[1] 0

Bagot_2017_SuscVsCTRL_HC_Genes<-Bagot_2017_SuscVsCTRL_HC$GeneID[Bagot_2017_SuscVsCTRL_HC$P.Value<0.005]
sum(duplicated(Bagot_2017_SuscVsCTRL_HC_Genes))
#[1] 0

StressGeneSets<-rbind(StressGeneSets, data.frame(GeneSymbol_Mouse=Bagot_2017_SuscVsCTRL_NACC_Genes, DatasetName=rep("Bagot_2017_StressSusceptibleVsCTRL_NACC", length(Bagot_2017_SuscVsCTRL_NACC_Genes)), stringsAsFactors=FALSE), data.frame(GeneSymbol_Mouse=Bagot_2017_SuscVsCTRL_HC_Genes, DatasetName=rep("Bagot_2017_StressSusceptibleVsCTRL_HC",length(Bagot_2017_SuscVsCTRL_HC_Genes)), stringsAsFactors=FALSE))

str(StressGeneSets)


str(Pena_2019_SuscVsCTRL_NACC)
#'data.frame':	75 obs. of  23 variables:
#It may be easiest just to grab all of these (they represent three stress conditions and two sexes)
sum(Pena_2019_SuscVsCTRL_NACC$pval_2Stresses_vs_ELS2Ctl<0.005)
#[1] 10
sum(Pena_2019_SuscVsCTRL_NACC$pval_2Stresses_vs_StdCtl<0.005)
#[1] 19
sum(Pena_2019_SuscVsCTRL_NACC$pval_StdAdult_vs_StdCtl<0.005)
#[1] 15

Pena_2019_SuscVsCTRL_NACC_Genes<-Pena_2019_SuscVsCTRL_NACC$gene_name[(Pena_2019_SuscVsCTRL_NACC$pval_2Stresses_vs_ELS2Ctl<0.005)|(Pena_2019_SuscVsCTRL_NACC$pval_2Stresses_vs_StdCtl<0.005)|(Pena_2019_SuscVsCTRL_NACC$pval_StdAdult_vs_StdCtl<0.005)]
sum(duplicated(Pena_2019_SuscVsCTRL_NACC_Genes))
#[1] 5

Pena_2019_SuscVsCTRL_NACC_Genes<-unique(Pena_2019_SuscVsCTRL_NACC_Genes)

StressGeneSets<-rbind(StressGeneSets, data.frame(GeneSymbol_Mouse=Pena_2019_SuscVsCTRL_NACC_Genes, DatasetName=rep("Pena_2019_StressVsCTRL_NACC", length(Pena_2019_SuscVsCTRL_NACC_Genes)), stringsAsFactors=FALSE))
                      

str(Gray_2014_S2_CRS_HC)
#'data.frame':	791 obs. of  12 variables:
sum(Gray_2014_S2_CRS_HC$p_CRS<0.005)
#[1] 84
#big difference, but not terrible.
sum(Gray_2014_S2_CRS_HC$p_CRS<0.005 & Gray_2014_S2_CRS_HC$Regulation_CRS=="up")
#[1] 34
sum(Gray_2014_S2_CRS_HC$p_CRS<0.005 & Gray_2014_S2_CRS_HC$Regulation_CRS=="down")
#[1] 50

Gray_2014_S2_CRS_HC_Genes_Up<-Gray_2014_S2_CRS_HC$Gene.Symbol[which(Gray_2014_S2_CRS_HC$p_CRS<0.005 & Gray_2014_S2_CRS_HC$Regulation_CRS=="up")]
Gray_2014_S2_CRS_HC_Genes_Down<-Gray_2014_S2_CRS_HC$Gene.Symbol[which(Gray_2014_S2_CRS_HC$p_CRS<0.005 & Gray_2014_S2_CRS_HC$Regulation_CRS=="down")]

Gray_2014_S2_CRS_HC_Genes_Up<-unique(Gray_2014_S2_CRS_HC_Genes_Up)
Gray_2014_S2_CRS_HC_Genes_Down<-unique(Gray_2014_S2_CRS_HC_Genes_Down)

str(Gray_2014_S3_CRSwFST_HC)
#'data.frame':	4448 obs. of  13 variables
sum(Gray_2014_S3_CRSwFST_HC$p._CRS_FST<0.005)
#[1] 1223 
#!!!
#Maybe trim it further so that it is more likely to make GSEA cut-offs:
sum(Gray_2014_S3_CRSwFST_HC$p._CRS_FST<0.0005)
#[1] 280

Gray_2014_S3_CRSwFST_HC_Genes_Up<-Gray_2014_S3_CRSwFST_HC$Symbol[which(Gray_2014_S3_CRSwFST_HC$p._CRS_FST<0.005 & Gray_2014_S3_CRSwFST_HC$Regulation_CRS_FST=="up")]
Gray_2014_S3_CRSwFST_HC_Genes_Down<-Gray_2014_S3_CRSwFST_HC$Symbol[which(Gray_2014_S3_CRSwFST_HC$p._CRS_FST<0.005 & Gray_2014_S3_CRSwFST_HC$Regulation_CRS_FST=="down")]

Gray_2014_S3_CRSwFST_HC_Genes_Up<-unique(Gray_2014_S3_CRSwFST_HC_Genes_Up)
Gray_2014_S3_CRSwFST_HC_Genes_Down<-unique(Gray_2014_S3_CRSwFST_HC_Genes_Down)

str(Gray_2014_Suppl_FSTvsCtrl_HC)
#'data.frame':	1369 obs. of  25 variables:
sum(Gray_2014_Suppl_FSTvsCtrl_HC$p<0.005)
#[1] 193
sum(Gray_2014_Suppl_FSTvsCtrl_HC$p<0.005 & Gray_2014_Suppl_FSTvsCtrl_HC$Regulation=="up")
#[1] 94

Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Up<-Gray_2014_Suppl_FSTvsCtrl_HC$Gene.Symbol[which(Gray_2014_Suppl_FSTvsCtrl_HC$p<0.005 & Gray_2014_Suppl_FSTvsCtrl_HC$Regulation=="up")]
Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Down<-Gray_2014_Suppl_FSTvsCtrl_HC$Gene.Symbol[which(Gray_2014_Suppl_FSTvsCtrl_HC$p<0.005 & Gray_2014_Suppl_FSTvsCtrl_HC$Regulation=="down")]  

Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Up<-unique(Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Up)
Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Down<-unique(Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Down)

str(Gray_2014_S7_CortVsVeh_HC)
#'data.frame':	1157 obs. of  13 variables:
sum(Gray_2014_S7_CortVsVeh_HC$p_AcuteCORT<0.005)
#[1] 102
sum(Gray_2014_S7_CortVsVeh_HC$p_AcuteCORT<0.005 & Gray_2014_S7_CortVsVeh_HC$Regulation_AcuteCORT=="up")
#[1] 57

Gray_2014_S7_CortVsVeh_HC_Genes_Up<-Gray_2014_S7_CortVsVeh_HC$Gene.Symbol[which(Gray_2014_S7_CortVsVeh_HC$p_AcuteCORT<0.005 & Gray_2014_S7_CortVsVeh_HC$Regulation_AcuteCORT=="up")]
Gray_2014_S7_CortVsVeh_HC_Genes_Down<-Gray_2014_S7_CortVsVeh_HC$Gene.Symbol[which(Gray_2014_S7_CortVsVeh_HC$p_AcuteCORT<0.005 & Gray_2014_S7_CortVsVeh_HC$Regulation_AcuteCORT=="down")]

Gray_2014_S7_CortVsVeh_HC_Genes_Up<-unique(Gray_2014_S7_CortVsVeh_HC_Genes_Up)
Gray_2014_S7_CortVsVeh_HC_Genes_Down<-unique(Gray_2014_S7_CortVsVeh_HC_Genes_Down)

str(Gray_2014_Suppl6_StressCortANOVA_HC)
#'data.frame':	8269 obs. of  10 variables:
sum(Gray_2014_Suppl6_StressCortANOVA_HC$p_fromANOVA<0.005)
#[1] 3841
sum(Gray_2014_Suppl6_StressCortANOVA_HC$p_fromANOVA<0.0005)
#[1] 1799
sum(Gray_2014_Suppl6_StressCortANOVA_HC$p_fromANOVA<0.00005)
#[1] 848
sum(Gray_2014_Suppl6_StressCortANOVA_HC$p_fromANOVA<0.000005)
#[1] 400
#I would need to make the cut off for this quite strict to have it make it into GSEA analyses, especially since it can't be split by directionality.

Gray_2014_Suppl6_StressCortANOVA_HC_Genes<-Gray_2014_Suppl6_StressCortANOVA_HC$Symbol[which(Gray_2014_Suppl6_StressCortANOVA_HC$p_fromANOVA<0.00005)]

Gray_2014_Suppl6_StressCortANOVA_HC_Genes<-unique(Gray_2014_Suppl6_StressCortANOVA_HC_Genes)

StressGeneSets<-rbind(StressGeneSets, 
                      data.frame(GeneSymbol_Mouse=Gray_2014_S2_CRS_HC_Genes_Up, DatasetName=rep("Gray_2014_ChronicRestraintStress_Upregulated_HC", length(Gray_2014_S2_CRS_HC_Genes_Up)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Mouse=Gray_2014_S2_CRS_HC_Genes_Down, DatasetName=rep("Gray_2014_ChronicRestraintStress_Downregulated_HC", length(Gray_2014_S2_CRS_HC_Genes_Down)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Mouse=Gray_2014_S3_CRSwFST_HC_Genes_Up, DatasetName=rep("Gray_2014_ChronicRestraintAndForcedSwimStress_Upregulated_HC", length(Gray_2014_S3_CRSwFST_HC_Genes_Up)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Mouse=Gray_2014_S3_CRSwFST_HC_Genes_Down, DatasetName=rep("Gray_2014_ChronicRestraintAndForcedSwimStress_Downregulated_HC", length(Gray_2014_S3_CRSwFST_HC_Genes_Down)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Mouse=Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Up, DatasetName=rep("Gray_2014_ForcedSwimStress_Upregulated_HC", length(Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Up)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Mouse=Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Down, DatasetName=rep("Gray_2014_ForcedSwimStress_Downregulated_HC", length(Gray_2014_Suppl_FSTvsCtrl_HC_Genes_Down)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Mouse=Gray_2014_S7_CortVsVeh_HC_Genes_Up, DatasetName=rep("Gray_2014_AcuteCorticosterone_Upregulated_HC", length(Gray_2014_S7_CortVsVeh_HC_Genes_Up)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Mouse=Gray_2014_S7_CortVsVeh_HC_Genes_Down, DatasetName=rep("Gray_2014_AcuteCorticosterone_Downregulated_HC", length(Gray_2014_S7_CortVsVeh_HC_Genes_Down)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Mouse=Gray_2014_Suppl6_StressCortANOVA_HC_Genes, DatasetName=rep("Gray_2014_AllStressConditions_HC", length(Gray_2014_Suppl6_StressCortANOVA_HC_Genes)), stringsAsFactors=FALSE)
                      )

write.csv(StressGeneSets, "StressGeneSets.csv")

#Double-checking whether any NAs snuck in there:
sum(is.na(StressGeneSets$GeneSymbol_Mouse)==TRUE)
#[1] 0

####################

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")

HippocampalDE_AcrossDatasets<-read.csv("TableS4_HippocampalDE_AcrossDatasets_wAVEDirection.csv", header=TRUE, stringsAsFactors = FALSE)

Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$bLR.vs..bHR==1)])
Birt_Hagenauer_2020_BredHighRespondersToNovelty_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$bLR.vs..bHR==-1)])

#Not enough DE genes in this dataset to split them by direction of effect:

Garafola_Hen_2014_CongenitallyLearnedHelpless_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$Congenitally.Helpless..cLH..and.Helpless.Resistant..cNLH.==1|HippocampalDE_AcrossDatasets$Congenitally.Helpless..cLH..and.Helpless.Resistant..cNLH.==-1)])

Blaveri_2010_FlindersSensitiveLine_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$Flinder.Sensitive.vs..Flinders.Resistant==1)])
Blaveri_2010_FlindersSensitiveLine_Downregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$Flinder.Sensitive.vs..Flinders.Resistant==-1)])

Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$Flinders.Sensitive.vs..Sprague.Dawley==1)])
Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$Flinders.Sensitive.vs..Sprague.Dawley==-1)])

DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$NIH.HS.High.Anxiety.vs..NIH.HS.Low.Anxiety==1)])
DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$NIH.HS.High.Anxiety.vs..NIH.HS.Low.Anxiety==-1)])

Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$RLA.vs..RHA==1)])
#Sabariego_2013_RomanLowAvoidanceRats_Downregulated_HC<-HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$RLA.vs..RHA==-1)]
#character(0)

Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$SLA.vs..SHA==1)])
Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$SLA.vs..SHA==-1)])

Meckes_2018_WKYvsF344_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$WKY.vs..F344==1)])
Meckes_2018_WKYvsF344_Downregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$WKY.vs..F344==-1)])

Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$WMI.vs..WLI..females.==1)])
Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$WMI.vs..WLI..females.==-1)])

Andrus_2012_WistarMoreImmobileRats_Upregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$WMI.vs..WLI..males.==1)])
Andrus_2012_WistarMoreImmobileRats_Downregulated_HC<-unique(HippocampalDE_AcrossDatasets$X[which(HippocampalDE_AcrossDatasets$WMI.vs..WLI..males.==-1)])


InternalizingBehaviorGeneSets<-data.frame(GeneSymbol_Rat=character(), DatasetName=character(), stringsAsFactors = FALSE)
str(InternalizingBehaviorGeneSets)

InternalizingBehaviorGeneSets<-rbind(InternalizingBehaviorGeneSets, 
                      data.frame(GeneSymbol_Rat=Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC, DatasetName=rep("Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC", length(Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Birt_Hagenauer_2020_BredHighRespondersToNovelty_Upregulated_HC, DatasetName=rep("Birt_Hagenauer_2020_BredHighRespondersToNovelty_Upregulated_HC", length(Birt_Hagenauer_2020_BredHighRespondersToNovelty_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Garafola_Hen_2014_CongenitallyLearnedHelpless_HC, DatasetName=rep("Garafola_Hen_2014_CongenitallyLearnedHelpless_HC", length(Garafola_Hen_2014_CongenitallyLearnedHelpless_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Blaveri_2010_FlindersSensitiveLine_Upregulated_HC, DatasetName=rep("Blaveri_2010_FlindersSensitiveLine_Upregulated_HC", length(Blaveri_2010_FlindersSensitiveLine_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Blaveri_2010_FlindersSensitiveLine_Downregulated_HC, DatasetName=rep("Blaveri_2010_FlindersSensitiveLine_Downregulated_HC", length(Blaveri_2010_FlindersSensitiveLine_Downregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC, DatasetName=rep("Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC", length(Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC, DatasetName=rep("Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC", length(Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC, DatasetName=rep("DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC", length(DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC, DatasetName=rep("DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC", length(DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC, DatasetName=rep("Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC", length(Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC, DatasetName=rep("Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC", length(Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC, DatasetName=rep("Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC", length(Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Meckes_2018_WKYvsF344_Upregulated_HC, DatasetName=rep("Meckes_2018_WKYvsF344_Upregulated_HC", length(Meckes_2018_WKYvsF344_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Meckes_2018_WKYvsF344_Downregulated_HC, DatasetName=rep("Meckes_2018_WKYvsF344_Downregulated_HC", length(Meckes_2018_WKYvsF344_Downregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC, DatasetName=rep("Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC", length(Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC, DatasetName=rep("Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC", length(Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Andrus_2012_WistarMoreImmobileRats_Upregulated_HC, DatasetName=rep("Andrus_2012_WistarMoreImmobileRats_Upregulated_HC", length(Andrus_2012_WistarMoreImmobileRats_Upregulated_HC)), stringsAsFactors=FALSE),
                      data.frame(GeneSymbol_Rat=Andrus_2012_WistarMoreImmobileRats_Downregulated_HC, DatasetName=rep("Andrus_2012_WistarMoreImmobileRats_Downregulated_HC", length(Andrus_2012_WistarMoreImmobileRats_Downregulated_HC)), stringsAsFactors=FALSE)
                      )
                      
Replicated_InternalizingBehavior_DifferentiallyExpressed_HC<-unique(InternalizingBehaviorGeneSets$GeneSymbol_Rat[duplicated(InternalizingBehaviorGeneSets$GeneSymbol_Rat)])
# [1] "Chd1l"        "Pcdhga1"      "Ak1"          "Drg1"         "Ghdc"         "Tmem144"      "Trhr"         "Zfp90"        "Trim45"       "Acox3"       
# [11] "Apln"         "Cav1"         "Eif2b1"       "Fcrl2"        "Htra1"        "Mfge8"        "Mknk1"        "Pdlim5"       "Pgm1"         "Rarres2"     
# [21] "Rbm3"         "Rhpn2"        "Robo3"        "Sh3bgr"       "Slc19a3"      "Slc25a18"     "Slc27a1"      "Sparcl1"      "Sun2"         "Tek"         
# [31] "Tmem176a"     "Tmem2"        "Zfp612"       "Slc35d3"      "Larp6"        "Cmc1"         "Abcg2"        "Ccdc107"      "Comt"         "Dcn"         
# [41] "Fam111a"      "Mgst1"        "P2rx4"        "Ppfibp2"      "Rab3b"        "Slc1a3"       "Sst"          "St5"          "Thy1"         "Tradd"       
# [51] "Tufm"         "C1qa"         "Fmo5"         "Fxyd7"        "Prss35"       "Rtn4ip1"      "Tes"          "Birc2"        "Camkk2"       "Ccnd1"       
# [61] "Cd48"         "Clic2"        "Commd6"       "Cyp4f4"       "F2r"          "Fkbp9"        "Fxyd5"        "Hapln1"       "Mfap3"        "Otub1"       
# [71] "Phf20l1"      "Pik3c2a"      "Pxdn"         "RGD1359508"   "Rnase4"       "Tmem100"      "Tubgcp2"      "Vwa5b2"       "Col6a1"       "Tmco5a"      
# [81] "Csad"         "Sncg"         "Frem3"        "Abcc9"        "Cd74"         "Gda"          "Ndufv2"       "Bmp4"         "Ccdc137"      "Ist1"        
# [91] "Naa50"        "Nudt4"        "Prrt2"        "Zfp821"       "Abhd13"       "Arhgap32"     "Arl16"        "Bmp3"         "Ccdc174"      "Ccndbp1"     
# [101] "Csdc2"        "Cyp2j3"       "Cyp4v3"       "Endod1"       "Eva1a"        "Fn1"          "Grm5"         "Hnrnpa1"      "Irgq"         "Jund"        
# [111] "Lcp1"         "Letm2"        "LOC102546433" "LOC103691688" "LOC688778"    "Lrrc61"       "Lsm14b"       "Mfap3l"       "Mfsd10"       "Nat8l"       
# [121] "Nqo1"         "Pfas"         "Pfkfb2"       "Pggt1b"       "Phip"         "Pkia"         "Ppp4c"        "Ptgr1"        "Ripk2"        "S100a10"     
# [131] "Sec23b"       "Sec23ip"      "Spint1"       "Sppl2a"       "Thrsp"        "Uggt2"        "Vps13c"       "Wipi2"        "Adam1a"       "Aqp9"        
# [141] "Arhgap20"     "Blvrb"        "Bmpr2"        "Bphl"         "Ccdc77"       "Cisd2"        "Grem2"        "Hcrt"         "Hspb1"        "Lrrc8b"      
# [151] "Magi1"        "Mis12"        "Mtmr1"        "Mx2"          "Pnpo"         "Ptprj"        "Rdx"          "Retsat"       "Serpinb9"     "Syt12"       
# [161] "Pvrl1"        "Spp1"         "Tgfbi"        "Mdh1"         "Ptpro"        "Slc3a1"       "Acad11"       "Akr1cl"       "C1qc"         "Gal"         
# [171] "Lcat"         "Medag"        "Tdg"          "Tnnt1"        "Trappc6a"     "Ucp2"         "Adpgk"        "Atg12"        "Bhlhe22"      "Cep95"       
# [181] "Chi3l1"       "Creg1"        "Cyp7b1"       "Deptor"       "Dnajc30"      "Dut"          "Ephx2"        "Gaa"          "Gad1"         "Gfra1"       
# [191] "Gldc"         "Golm1"        "Ilf3"         "Kidins220"    "Lancl1"       "LOC100911253" "LOC102552880" "LOC103692950" "Mdm4"         "Mpeg1"       
# [201] "Opa1"         "Pcdh8"        "Prelp"        "Prepl"        "Pwp2"         "Rnf6"         "Rsrp1"        "Stxbp1"       "Syt1"         "Zfp180"      
# [211] "Acat1"        "Akap5"        "Blnk"         "C3"           "Cib2"         "Dpyd"         "Dsp"          "F11r"         "Fbln2"        "Fntb"        
# [221] "Itga7"        "LOC100174910" "Lum"          "Nqo2"         "Perp"         "Plekhh1"      "RGD1309104"   "Pcdh15"       "Pqlc2"        "RT1-Da"      
# [231] "Sostdc1"      "Id3"          "Prlr"         "Vdac1"        "Celsr3"       "Sec63"        "Clec9a"       "Traf6"        "Slc22a7"      "Tspan7"      
# [241] "Mppe1"        "Acbd4"        "Capg"         "Ccdc62"       "Chst2"        "Kdm5b"        "Ltk"          "Numbl"        "Slc24a3"      "Crtac1"      
# [251] "Fxyd3"        "RGD1305680"   "Abcb10"       "Acvr1c"       "Dtnb"         "Kcnab2"       "Ophn1"        "Pik3r1"       "Serinc2"      "Strbp"       
# [261] "Tef"          "Pigh"         "Atp1b2"       "Pdzd2"        "Bcat1"        "Cav2"         "Ehd4"         "Etfa"         "Igfbp5"       "Inpp5f"      
# [271] "P2ry14"       "Slc35f1"      "Tcn2"         "Thoc1"        "Inadl"        "Psme4"        "RGD1308772"   "Agtpbp1"      "Calml4"       "Clic6"       
# [281] "Creb1"        "Dab2"         "F5"           "Gls"          "Kalrn"        "Krt8"         "Mcat"         "Mfrp"         "Nedd4"        "Pak1"        
# [291] "Sez6l2"       "Slc31a1"      "Slc6a20"      "Ttr"          "Yipf3"    

InternalizingBehaviorGeneSets<-rbind(InternalizingBehaviorGeneSets, 
                                     data.frame(GeneSymbol_Rat=Replicated_InternalizingBehavior_DifferentiallyExpressed_HC, DatasetName=rep("Replicated_InternalizingBehavior_DifferentiallyExpressed_HC", length(Replicated_InternalizingBehavior_DifferentiallyExpressed_HC)), stringsAsFactors=FALSE))

write.csv(InternalizingBehaviorGeneSets, "InternalizingBehaviorGeneSets.csv")


#############

#Genes with strong evidence supporting their differential expression in psychiatric disorder (meta-analysis):

Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis<-read.csv("Gandal_2018_PsychTranscriptomeVsGenetics_Suppl_TableS1_MicroarrayMetaAnalysis.csv", header=TRUE, stringsAsFactors = FALSE)

Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis<-read.csv("Gandal_2018_PsychTranscriptomeVsGenetics_Suppl_TableS1_RNASeqMetaAnalysis.csv", header=TRUE, stringsAsFactors = FALSE)

str(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis)

#these two analyses produced such large gene lists I needed to add a p-value cut-off to make an effective gene set:
length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.beta_log2FC>0)])
#[1] 1096


length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.P.value<0.001)])
#[1] 532

Gandal_2018_AutismSpectrumDisorder_Upregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.P.value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.P.value<0.001)])
#[1] 616

Gandal_2018_AutismSpectrumDisorder_Downregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$ASD.P.value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.P.value<0.001)])
#[1] 461

Gandal_2018_Schizophrenia_Upregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.P.value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.P.value<0.001)])
#[1] 636

Gandal_2018_Schizophrenia_Downregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$SCZ.P.value<0.001)])
  
length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.P.value<0.001)])
#[1] 178

Gandal_2018_BipolarDisorder_Upregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.P.value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.P.value<0.001)])
#[1] 191

Gandal_2018_BipolarDisorder_Downregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$BD.P.value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.P.value<0.001)])
#[1] 183

Gandal_2018_MajorDepressiveDisorder_Upregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.P.value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.P.value<0.001)])
#[1] 145

Gandal_2018_MajorDepressiveDisorder_Downregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$MDD.P.value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.P.value<0.001)])
#[1] 289

Gandal_2018_AlcoholAbuseDisorder_Upregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.beta_log2FC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.P.value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.P.value<0.001)])
#[1] 312

Gandal_2018_AlcoholAbuseDisorder_Downregulated_Cortex<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$hgnc_symbol[which(Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.FDR<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.beta_log2FC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_MicroarrayMetaAnalysis$AAD.P.value<0.001)])

#Hmmm... I forgot that the RNA-Seq data is only given in Ensembl notation. I'll need to decipher that to gene symbols to use it for a gene set:
length(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis$X)

library(org.Hs.eg.db)

x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

xx[1]
# $`1`
# [1] "A1BG"
names(xx[1])
#[1] "1"
names(xx[800])
#[1] "100128946"
#Hmm...
#Apparently this package only maps to Entrez, not Ensembl. :(

columns(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
# [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT" 

#Ah - I think I need to map both Symbol and Ensembl back to Entrez to align them.
#Or maybe not?

keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
# [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT"     

uniKeys <- head(keys(org.Hs.eg.db, keytype="ENSEMBL"))
cols <- c("SYMBOL")
select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="ENSEMBL")
# 'select()' returned 1:1 mapping between keys and columns
# ENSEMBL   SYMBOL
# 1 ENSG00000121410     A1BG
# 2 ENSG00000175899      A2M
# 3 ENSG00000256069    A2MP1
# 4 ENSG00000171428     NAT1
# 5 ENSG00000156006     NAT2
# 6 ENSG00000196136 SERPINA3

#Lovely.

colnames(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis)[1]<-"ENSEMBL"
str(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis)
#'data.frame':	15823 obs. of  9 variables:

EnsemblVsGeneSymbol_Human<-select(org.Hs.eg.db, keys=Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis$ENSEMBL, columns=cols, keytype="ENSEMBL")
str(EnsemblVsGeneSymbol_Human)
# 'data.frame':	16187 obs. of  2 variables:
# $ ENSEMBL: chr  "ENSG00000000003" "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" ...
# $ SYMBOL : chr  "TSPAN6" "DPM1" "SCYL3" "C1orf112" ...

Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols<-join(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis, EnsemblVsGeneSymbol_Human, by="ENSEMBL", type="left", match="all")
str(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols)
#'data.frame':	15823 obs. of  10 variables:


length(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SYMBOL[which(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.adj.P.Val<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.logFC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.P.Value<0.001)])
#[1] 344

Gandal_2018_Schizophrenia_Upregulated_PFC<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SYMBOL[which(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.adj.P.Val<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.logFC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.P.Value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SYMBOL[which(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.adj.P.Val<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.logFC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.P.Value<0.001)])
#[1] 535

Gandal_2018_Schizophrenia_Downregulated_PFC<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SYMBOL[which(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.adj.P.Val<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.logFC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SCZ.P.Value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SYMBOL[which(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.adj.P.Val<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.logFC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.P.Value<0.001)])
#[1] 77

Gandal_2018_BipolarDisorder_Upregulated_PFC<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SYMBOL[which(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.adj.P.Val<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.logFC>0 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.P.Value<0.001)])

length(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SYMBOL[which(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.adj.P.Val<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.logFC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.P.Value<0.001)])  
#[1] 107

Gandal_2018_BipolarDisorder_Downregulated_PFC<-unique(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$SYMBOL[which(Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.adj.P.Val<0.05 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.logFC<0 & Gandal_2018_PsychTranscriptomeVsGenetics_RNASeqMetaAnalysis_wSymbols$BD.P.Value<0.001)]) 


#Can I make the combining gene sets code a little smoother?
Gandal_2018_GeneSets<-list(Gandal_2018_Schizophrenia_Downregulated_PFC,Gandal_2018_Schizophrenia_Upregulated_PFC, Gandal_2018_Schizophrenia_Downregulated_Cortex, Gandal_2018_Schizophrenia_Upregulated_Cortex,Gandal_2018_BipolarDisorder_Downregulated_Cortex,Gandal_2018_BipolarDisorder_Downregulated_PFC, Gandal_2018_BipolarDisorder_Upregulated_Cortex, Gandal_2018_BipolarDisorder_Upregulated_PFC, Gandal_2018_MajorDepressiveDisorder_Downregulated_Cortex, Gandal_2018_MajorDepressiveDisorder_Upregulated_Cortex, Gandal_2018_AlcoholAbuseDisorder_Downregulated_Cortex,Gandal_2018_AlcoholAbuseDisorder_Upregulated_Cortex, Gandal_2018_AutismSpectrumDisorder_Downregulated_Cortex, Gandal_2018_AutismSpectrumDisorder_Upregulated_Cortex)

Gandal_2018_GeneSet_Names<-list("Gandal_2018_Schizophrenia_Downregulated_PFC","Gandal_2018_Schizophrenia_Upregulated_PFC", "Gandal_2018_Schizophrenia_Downregulated_Cortex", "Gandal_2018_Schizophrenia_Upregulated_Cortex","Gandal_2018_BipolarDisorder_Downregulated_Cortex","Gandal_2018_BipolarDisorder_Downregulated_PFC", "Gandal_2018_BipolarDisorder_Upregulated_Cortex", "Gandal_2018_BipolarDisorder_Upregulated_PFC", "Gandal_2018_MajorDepressiveDisorder_Downregulated_Cortex", "Gandal_2018_MajorDepressiveDisorder_Upregulated_Cortex", "Gandal_2018_AlcoholAbuseDisorder_Downregulated_Cortex","Gandal_2018_AlcoholAbuseDisorder_Upregulated_Cortex", "Gandal_2018_AutismSpectrumDisorder_Downregulated_Cortex", "Gandal_2018_AutismSpectrumDisorder_Upregulated_Cortex")

HumanPsychiatricGeneSets<-data.frame(GeneSymbol_Human=character(), DatasetName=character(), stringsAsFactors = FALSE)
str(HumanPsychiatricGeneSets)

Gandal_2018_GeneSets[[1]]

for(i in c(1:length(Gandal_2018_GeneSet_Names))){
HumanPsychiatricGeneSets<-rbind(HumanPsychiatricGeneSets, 
                                data.frame(GeneSymbol_Human=as.character(Gandal_2018_GeneSets[[i]]), DatasetName=as.character(rep(Gandal_2018_GeneSet_Names[[i]], length(Gandal_2018_GeneSets[[i]]))), stringsAsFactors=FALSE))
}

str(HumanPsychiatricGeneSets)
# 'data.frame':	4420 obs. of  2 variables:
#   $ GeneSymbol_Human: chr  "BAD" "DBNDD1" "MSL3" "TSPOAP1" ...
# $ DatasetName     : chr  "Gandal_2018_Schizophrenia_Downregulated_PFC" "Gandal_2018_Schizophrenia_Downregulated_PFC" "Gandal_2018_Schizophrenia_Downregulated_PFC" "Gandal_2018_Schizophrenia_Downregulated_PFC" ...

#Note: There are definitely some NAs in here that need to be removed - maybe after adding species-specific annotation.

HumanPsychiatricGeneSets<-HumanPsychiatricGeneSets[is.na(HumanPsychiatricGeneSets$GeneSymbol_Human)==FALSE,]
str(HumanPsychiatricGeneSets)
#'data.frame':	4416 obs. of  2 variables:
# $ GeneSymbol_Human: chr  "BAD" "DBNDD1" "MSL3" "TSPOAP1" ...
# $ DatasetName     : chr  "Gandal_2018_Schizophrenia_Downregulated_PFC" "Gandal_2018_Schizophrenia_Downregulated_PFC" "Gandal_2018_Schizophrenia_Downregulated_PFC" "Gandal_2018_Schizophrenia_Downregulated_PFC" ...

setwd("/Users/mhh/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")
write.csv(HumanPsychiatricGeneSets, "HumanPsychiatricGeneSets.csv")


#It occurs to me that I should probably make a better ortholog version of our original cell type marker database too:

BrainInABlenderDatabase<-read.csv("BrainInABlenderDatabase_DateGenesFixed.csv", header=TRUE, stringsAsFactors = FALSE)
str(BrainInABlenderDatabase)

BrainInABlenderDatabase_Human<-BrainInABlenderDatabase[BrainInABlenderDatabase$Species=="Human",]
BrainInABlenderDatabase_Mouse<-BrainInABlenderDatabase[BrainInABlenderDatabase$Species=="Mice",]


#################

#Fixing lurking issues with genes with date names (since that problem in the gene symbol annotation was only changed mid-2020):

DateGeneDecoder<-read.delim("DateGeneDecoder.txt", sep=",", header=TRUE, stringsAsFactors = FALSE)
str(DateGeneDecoder)
# 'data.frame':	28 obs. of  3 variables:
# $ OldSymbol: chr  "Sept1" "Sept2" "Sept3" "Sept4" ...
# $ NewSymbol: chr  "Septin1" "Septin2" "Septin3" "Septin4" ...
# $ ExcelDate: chr  "1-Sep" "2-Sep" "3-Sep" "4-Sep" ...

DateGeneDecoder_Human<-DateGeneDecoder
DateGeneDecoder_Human$OldSymbol<-toupper(DateGeneDecoder$OldSymbol)
DateGeneDecoder_Human$NewSymbol<-toupper(DateGeneDecoder$NewSymbol)


StressGeneSets[StressGeneSets$GeneSymbol_Mouse%in%DateGeneDecoder$OldSymbol|StressGeneSets$GeneSymbol_Mouse%in%DateGeneDecoder$ExcelDate,]
# GeneSymbol_Mouse                                 DatasetName
# 1492            5-Sep Gray_2014_ForcedSwimStress_Downregulated_HC

StressGeneSets$GeneSymbol_Mouse[StressGeneSets$GeneSymbol_Mouse=="5-Sep"]<-"Septin5"

#Well that one was easy...

write.csv(StressGeneSets, "StressGeneSets_DateGenesFixed.csv")

HumanPsychiatricGeneSets[HumanPsychiatricGeneSets$GeneSymbol_Human%in%DateGeneDecoder_Human$OldSymbol|HumanPsychiatricGeneSets$GeneSymbol_Human%in%DateGeneDecoder_Human$ExcelDate,]
#     GeneSymbol_Human                                           DatasetName
# 67             MARCH2           Gandal_2018_Schizophrenia_Downregulated_PFC
# 240            SEPT11           Gandal_2018_Schizophrenia_Downregulated_PFC
# 775            MARCH3             Gandal_2018_Schizophrenia_Upregulated_PFC
# 2275            2-Mar        Gandal_2018_BipolarDisorder_Upregulated_Cortex
# 2359            3-Mar        Gandal_2018_BipolarDisorder_Upregulated_Cortex
# 3235            5-Sep   Gandal_2018_AlcoholAbuseDisorder_Upregulated_Cortex
# 4315            3-Mar Gandal_2018_AutismSpectrumDisorder_Upregulated_Cortex

HumanPsychiatricGeneSets$GeneSymbol_Human[HumanPsychiatricGeneSets$GeneSymbol_Human=="MARCH2"]<-"MARCHF2"
HumanPsychiatricGeneSets$GeneSymbol_Human[HumanPsychiatricGeneSets$GeneSymbol_Human=="SEPT11"]<-"SEPTIN11"
HumanPsychiatricGeneSets$GeneSymbol_Human[HumanPsychiatricGeneSets$GeneSymbol_Human=="2-Mar"]<-"MARCHF2"
HumanPsychiatricGeneSets$GeneSymbol_Human[HumanPsychiatricGeneSets$GeneSymbol_Human=="3-Mar"]<-"MARCHF3"
HumanPsychiatricGeneSets$GeneSymbol_Human[HumanPsychiatricGeneSets$GeneSymbol_Human=="5-Sep"]<-"SEPTIN5"

#That was also easy.

write.csv(HumanPsychiatricGeneSets, "HumanPsychiatricGeneSets_DateGenesFixed.csv")

#I should probably double-check that there aren't any more lurking in the DropViz or Zeisel data:

DropViz_HC_CellTypeClusterMarkers$gene[DropViz_HC_CellTypeClusterMarkers$gene%in%DateGeneDecoder$OldSymbol|DropViz_HC_CellTypeClusterMarkers$gene%in%DateGeneDecoder$ExcelDate]
#character(0)
DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene%in%DateGeneDecoder$OldSymbol|DropViz_STR_CellTypeClusterMarkers$gene%in%DateGeneDecoder$ExcelDate]
#character(0)

#We got them all. :)

BrainInABlenderDatabase_Human$Gene.Symbol..Human.[BrainInABlenderDatabase_Human$Gene.Symbol..Human.%in%DateGeneDecoder_Human$OldSymbol|BrainInABlenderDatabase_Human$Gene.Symbol..Human.%in%DateGeneDecoder_Human$ExcelDate]
#character(0)
BrainInABlenderDatabase_Human$Gene.Symbol..Mouse.[BrainInABlenderDatabase_Human$Gene.Symbol..Mouse.%in%DateGeneDecoder$OldSymbol|BrainInABlenderDatabase_Human$Gene.Symbol..Mouse.%in%DateGeneDecoder$ExcelDate]
#character(0)

BrainInABlenderDatabase_Mouse$Gene.Symbol..Mouse.[BrainInABlenderDatabase_Mouse$Gene.Symbol..Mouse.%in%DateGeneDecoder$OldSymbol|BrainInABlenderDatabase_Mouse$Gene.Symbol..Mouse.%in%DateGeneDecoder$ExcelDate]
#[1] "Nov" "Nov" "Nov"
BrainInABlenderDatabase_Mouse$Gene.Symbol..Mouse.[BrainInABlenderDatabase_Mouse$Gene.Symbol..Mouse.=="Nov"]<-"Ccn3"

BrainInABlenderDatabase_Mouse$Gene.Symbol..Human.[BrainInABlenderDatabase_Mouse$Gene.Symbol..Human.%in%DateGeneDecoder_Human$OldSymbol|BrainInABlenderDatabase_Mouse$Gene.Symbol..Human.%in%DateGeneDecoder_Human$ExcelDate]
[1] "NOV" "NOV" "NOV"
BrainInABlenderDatabase_Mouse$Gene.Symbol..Human.[BrainInABlenderDatabase_Mouse$Gene.Symbol..Human.=="NOV"]<-"CCN3"

#It is... a little shocking to me that there aren't more. Skeptic that I am, it makes me nervous that the date gene symbols actually converted to dates and then were forced into character format so they turned into numbers. Because that's the stupid sort of shit that Excel does.
#So I went back to Excel and took a peek and that is exactly what happened. 
#So I converted the numbers back to dates (hoping that this version of excel is using the same numeric date references as whichever version created the database...). The dates all line up with March and Sept genes, so I fixed them by hand hoping for the best and reloaded the data.


str(InternalizingBehaviorGeneSets)
InternalizingBehaviorGeneSets[InternalizingBehaviorGeneSets$GeneSymbol_Rat%in%DateGeneDecoder$OldSymbol|InternalizingBehaviorGeneSets$GeneSymbol_Rat%in%DateGeneDecoder$ExcelDate,]
# GeneSymbol_Rat                                         DatasetName
# 365             Nov   Blaveri_2010_FlindersSensitiveLine_Upregulated_HC
# 431           Sept6   Blaveri_2010_FlindersSensitiveLine_Upregulated_HC
# 836          March8 Blaveri_2010_FlindersSensitiveLine_Downregulated_HC
# 999           Sept2 Blaveri_2010_FlindersSensitiveLine_Downregulated_HC
# 1000          Sept5 Blaveri_2010_FlindersSensitiveLine_Downregulated_HC
# 1001          Sept8 Blaveri_2010_FlindersSensitiveLine_Downregulated_HC
# 1905         Sept10                Meckes_2018_WKYvsF344_Upregulated_HC

InternalizingBehaviorGeneSets$GeneSymbol_Rat[InternalizingBehaviorGeneSets$GeneSymbol_Rat=="Sept6"]<-"Septin6"
InternalizingBehaviorGeneSets$GeneSymbol_Rat[InternalizingBehaviorGeneSets$GeneSymbol_Rat=="Sept2"]<-"Septin2"
InternalizingBehaviorGeneSets$GeneSymbol_Rat[InternalizingBehaviorGeneSets$GeneSymbol_Rat=="Sept5"]<-"Septin5"
InternalizingBehaviorGeneSets$GeneSymbol_Rat[InternalizingBehaviorGeneSets$GeneSymbol_Rat=="Sept8"]<-"Septin8"
InternalizingBehaviorGeneSets$GeneSymbol_Rat[InternalizingBehaviorGeneSets$GeneSymbol_Rat=="Sept10"]<-"Septin10"
InternalizingBehaviorGeneSets$GeneSymbol_Rat[InternalizingBehaviorGeneSets$GeneSymbol_Rat=="March8"]<-"Marchf8"
InternalizingBehaviorGeneSets$GeneSymbol_Rat[InternalizingBehaviorGeneSets$GeneSymbol_Rat=="Nov"]<-"Ccn3"

write.csv(InternalizingBehaviorGeneSets, "InternalizingBehaviorGeneSets_DateGenesFixed.csv")


#How about MsigDb:

str(msigdb.v7.3.symbols)
#'data.frame':	49589 obs. of  142 variables:
#'That's weird - that data.frame seems to be truncated at 142 genes per category.
#I just opened it up in Excel and there are definitely SEPTIN genes - so hypothetically the gene symbol annotation in these files has the date problem fixed.

#I'm pretty sure my old home-brewed HC specific GMT has date problems. Let's look:
#Looks like I may have already gone in and fixed it in Excel. Excellent.

str(CombinedHCspecific)
sum("CDS1"%in%CombinedHCspecific)#Sanity check - apparently this code doesn't work unless the column is specified, which is awkward.
sum("CDS1"%in%as.matrix(CombinedHCspecific[,c(3:ncol(CombinedHCspecific))]))
#[1] 1
#That code works. :)

TempMatrix<-as.matrix(CombinedHCspecific[,c(3:ncol(CombinedHCspecific))])

TempMatrix[TempMatrix%in%DateGeneDecoder_Human$OldSymbol|TempMatrix%in%DateGeneDecoder_Human$ExcelDate]
#[1] "NOV"   "NOV"   "SEPT9"

TempMatrix[TempMatrix=="NOV"]<-"CCN3"
TempMatrix[TempMatrix=="SEPT9"]<-"SEPTIN9"

CombinedHCspecific[,c(3:ncol(CombinedHCspecific))]<-TempMatrix
str(CombinedHCspecific)

write.table(CombinedHCspecific, "TableS1_CombinedHCspecific_DateGenesFixed.gmt.txt", sep="\t", col.names = FALSE)


#how about some of the other gmts?

c2.all.v7.3.symbols
TempMatrix<-as.matrix(c2.all.v7.3.symbols[,c(3:ncol(c2.all.v7.3.symbols))])
TempMatrix[TempMatrix%in%DateGeneDecoder_Human$OldSymbol|TempMatrix%in%DateGeneDecoder_Human$ExcelDate]
#character(0)

#Out of laziness, I just reran the same code for all of the others and they all look fine - up to date annotation.
TempMatrix<-as.matrix(c3.all.v7.3.symbols[,c(3:ncol(c3.all.v7.3.symbols))])
TempMatrix<-as.matrix(c5.all.v7.3.symbols[,c(3:ncol(c5.all.v7.3.symbols))])
TempMatrix<-as.matrix(c7.all.v7.3.symbols[,c(3:ncol(c7.all.v7.3.symbols))])
TempMatrix<-as.matrix(c8.all.v7.3.symbols[,c(3:ncol(c8.all.v7.3.symbols))])
TempMatrix<-as.matrix(h.all.v7.3.symbols[,c(3:ncol(h.all.v7.3.symbols))])
#character(0) for all after running the code above


#Zeisel Pyramidal neurons:

str(Zeisel_CA1PyramidalIndex)
Zeisel_CA1PyramidalIndex[Zeisel_CA1PyramidalIndex$Symbol_Mouse_AsUpperCase%in%DateGeneDecoder_Human$OldSymbol|Zeisel_CA1PyramidalIndex$Symbol_Mouse_AsUpperCase%in%DateGeneDecoder$ExcelDate,]
# Symbol_Mouse_AsUpperCase                       Dataset
# 223                   MARCH3 Zeisel_CA1PyramidalNeurons_HC
# 327                    SEPT9 Zeisel_CA1PyramidalNeurons_HC

Zeisel_CA1PyramidalIndex$Symbol_Mouse_AsUpperCase[Zeisel_CA1PyramidalIndex$Symbol_Mouse_AsUpperCase=="MARCH3"]<-"MARCHF3"
Zeisel_CA1PyramidalIndex$Symbol_Mouse_AsUpperCase[Zeisel_CA1PyramidalIndex$Symbol_Mouse_AsUpperCase=="SEPT9"]<-"SEPTIN9"


str(Zeisel_CA1PyramidalIndex)
setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")
write.csv(Zeisel_CA1PyramidalIndex, "Zeisel_CA1PyramidalIndex_FixedDateGenes.csv")

str(DropViz_HC_CellTypeClusterMarkers)

DropViz_HC_CellTypeClusterMarkers$gene[DropViz_HC_CellTypeClusterMarkers$gene%in%DateGeneDecoder$OldSymbol|DropViz_HC_CellTypeClusterMarkers$gene%in%DateGeneDecoder$ExcelDate]
#[1] "Sept10"  "Sept10"  "March10" "Sept4"   "Sept11"  "Nov"     "Sept4"   "March8"  "Sept7"  

DropViz_HC_CellTypeClusterMarkers$gene[DropViz_HC_CellTypeClusterMarkers$gene=="Sept10"]<-"Septin10"
DropViz_HC_CellTypeClusterMarkers$gene[DropViz_HC_CellTypeClusterMarkers$gene=="Sept4"]<-"Septin4"
DropViz_HC_CellTypeClusterMarkers$gene[DropViz_HC_CellTypeClusterMarkers$gene=="Sept11"]<-"Septin11"
DropViz_HC_CellTypeClusterMarkers$gene[DropViz_HC_CellTypeClusterMarkers$gene=="Sept7"]<-"Septin7"
DropViz_HC_CellTypeClusterMarkers$gene[DropViz_HC_CellTypeClusterMarkers$gene=="March10"]<-"Marchf10"
DropViz_HC_CellTypeClusterMarkers$gene[DropViz_HC_CellTypeClusterMarkers$gene=="March8"]<-"Marchf8"
DropViz_HC_CellTypeClusterMarkers$gene[DropViz_HC_CellTypeClusterMarkers$gene=="Nov"]<-"Ccn3"

DropViz_HC_CellTypeClusterMarkers_Concise$gene[DropViz_HC_CellTypeClusterMarkers_Concise$gene%in%DateGeneDecoder$OldSymbol|DropViz_HC_CellTypeClusterMarkers_Concise$gene%in%DateGeneDecoder$ExcelDate]
#[1] "Sept4"
DropViz_HC_CellTypeClusterMarkers_Concise$gene[DropViz_HC_CellTypeClusterMarkers_Concise$gene=="Sept4"]<-"Septin4"

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")

write.csv(DropViz_HC_CellTypeClusterMarkers, "DropViz_HC_CellTypeClusterMarkers_FixedDateGenes.csv")
write.csv(DropViz_HC_CellTypeClusterMarkers_Concise, "DropViz_HC_CellTypeClusterMarkers_Concise_FixedDateGenes.csv")

str(DropViz_STR_CellTypeClusterMarkers)

DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene%in%DateGeneDecoder$OldSymbol|DropViz_STR_CellTypeClusterMarkers$gene%in%DateGeneDecoder$ExcelDate]
#[1] "Sept10"  "March10" "March1"  "Sept4"   "Sept10"  "Sept4"   "March8"  "Sept8"   "Sept7"   "Sept2"  

DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene=="Sept10"]<-"Septin10"
DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene=="Sept4"]<-"Septin4"
DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene=="Sept7"]<-"Septin7"
DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene=="Sept8"]<-"Septin8"
DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene=="Sept2"]<-"Septin2"

DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene=="March10"]<-"Marchf10"
DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene=="March1"]<-"Marchf1"
DropViz_STR_CellTypeClusterMarkers$gene[DropViz_STR_CellTypeClusterMarkers$gene=="March8"]<-"Marchf8"


DropViz_STR_CellTypeClusterMarkers_Concise$gene[DropViz_STR_CellTypeClusterMarkers_Concise$gene%in%DateGeneDecoder$OldSymbol|DropViz_STR_CellTypeClusterMarkers_Concise$gene%in%DateGeneDecoder$ExcelDate]
#[1] "Sept4"
DropViz_STR_CellTypeClusterMarkers_Concise$gene[DropViz_STR_CellTypeClusterMarkers_Concise$gene=="Sept4"]<-"Septin4"

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")

write.csv(DropViz_STR_CellTypeClusterMarkers, "DropViz_STR_CellTypeClusterMarkers_FixedDateGenes.csv")
write.csv(DropViz_STR_CellTypeClusterMarkers_Concise, "DropViz_STR_CellTypeClusterMarkers_Concise_FixedDateGenes.csv")


###############################################

#Determining formal orthologs between the data from different species:

setwd("~/Documents/Microarray Gen/GeneOrthology_HumanRatMouse")

#I downloaded formal ortholog information from: http://www.informatics.jax.org/homology.shtml
#Let's reorganize it to make it easier to work with:
Orthologs<-read.delim("InformaticsJax_OrthologDatabase_20210228.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
str(Orthologs)
Orthologs_Mice<-Orthologs[Orthologs$Common.Organism.Name=="mouse, laboratory",]
str(Orthologs_Mice)
#'data.frame':	20891 obs. of  12 variables:

colnames(Orthologs_Mice)[c(2:12)]<-paste(colnames(Orthologs_Mice)[c(2:12)], "_Mouse", sep="")

Orthologs_Rat<-Orthologs[Orthologs$Common.Organism.Name=="rat",]
str(Orthologs_Rat)
#"'data.frame':	20606 obs. of  12 variables:"

colnames(Orthologs_Rat)[c(2:12)]<-paste(colnames(Orthologs_Rat)[c(2:12)], "_Rat", sep="")

Orthologs_Mice_Rats<-join(Orthologs_Mice, Orthologs_Rat, by="HomoloGene.ID", type="left", match="all")
str(Orthologs_Mice_Rats)
#'data.frame':	22814 obs. of  23 variables:


Orthologs_Human<-Orthologs[Orthologs$Common.Organism.Name=="human",]
str(Orthologs_Human)
#"'data.frame':	19124 obs. of  12 variables:"

colnames(Orthologs_Human)[c(2:12)]<-paste(colnames(Orthologs_Human)[c(2:12)], "_Human", sep="")

Orthologs_Mice_Humans<-join(Orthologs_Mice, Orthologs_Human, by="HomoloGene.ID", type="left", match="all")
str(Orthologs_Mice_Humans)
#'data.frame':	21021 obs. of  23 variables:

Orthologs_Rats_Humans<-join(Orthologs_Rat, Orthologs_Human, by="HomoloGene.ID", type="left", match="all")
str(Orthologs_Rats_Humans)
#'data.frame':	20728 obs. of  23 variables:

#Out of curiousity, how much does this actually gain me over just using gene symbol equivalency?
sum(Orthologs_Mice_Rats$Symbol_Mouse==Orthologs_Mice_Rats$Symbol_Rat, na.rm=TRUE)
#[1] 15728
sum(is.na(Orthologs_Mice_Rats$Symbol_Rat)==TRUE|is.na(Orthologs_Mice_Rats$Symbol_Mouse)==TRUE)
#[1] 2176
15728/(22814-2176)
#[1] 0.7620893
#Interesting - symbols are only equivalent for ~3/4 of the genes (according to this database)

Orthologs_Mice_Rats_Symbols<-data.frame(Orthologs_Mice_Rats$Symbol_Mouse, Orthologs_Mice_Rats$Symbol_Rat, stringsAsFactors = FALSE)
colnames(Orthologs_Mice_Rats_Symbols)<-c("Symbol_Mouse", "Symbol_Rat")
sum(is.na(Orthologs_Mice_Rats_Symbols$Symbol_Rat))
#[1] 2176

Orthologs_Mice_Rats_Symbols_NoNA<-Orthologs_Mice_Rats_Symbols[is.na(Orthologs_Mice_Rats_Symbols$Symbol_Rat)==FALSE & is.na(Orthologs_Mice_Rats$Symbol_Mouse)==FALSE,]
str(Orthologs_Mice_Rats_Symbols_NoNA)
# 'data.frame':	20638 obs. of  2 variables:
# $ Symbol_Mouse: chr  "Acadm" "Acadvl" "Acat1" "Acvr1" ...
# $ Symbol_Rat  : chr  "Acadm" "Acadvl" "Acat1" "Acvr1" ...

Orthologs_Mice_Humans_Symbols<-data.frame(Orthologs_Mice_Humans$Symbol_Mouse, Orthologs_Mice_Humans$Symbol_Human, stringsAsFactors = FALSE)
colnames(Orthologs_Mice_Humans_Symbols)<-c("Symbol_Mouse", "Symbol_Human")
sum(is.na(Orthologs_Mice_Humans_Symbols$Symbol_Human)==TRUE | is.na(Orthologs_Mice_Humans$Symbol_Mouse)==TRUE)
#[1] 3705

Orthologs_Mice_Humans_Symbols_NoNA<-Orthologs_Mice_Humans_Symbols[is.na(Orthologs_Mice_Humans_Symbols$Symbol_Human)==FALSE & is.na(Orthologs_Mice_Humans$Symbol_Mouse)==FALSE,]
str(Orthologs_Mice_Humans_Symbols_NoNA)
#'data.frame':	17316 obs. of  2 variables:
# $ Symbol_Mouse: chr  "Acadm" "Acadvl" "Acat1" "Acvr1" ...
# $ Symbol_Human: chr  "ACADM" "ACADVL" "ACAT1" "ACVR1" ...

#hmmm... let's double check whether the date genes are problematic in this dataset:
Orthologs_Mice_Humans_Symbols_NoNA$Symbol_Mouse[Orthologs_Mice_Humans_Symbols_NoNA$Symbol_Mouse=="Septin9"]
#[1] "Septin9"
#Nope, looks like they have been updated. *phew*

Orthologs_Rats_Humans_Symbols<-data.frame(Orthologs_Rats_Humans$Symbol_Rat, Orthologs_Rats_Humans$Symbol_Human, stringsAsFactors = FALSE)
colnames(Orthologs_Rats_Humans_Symbols)<-c("Symbol_Rat", "Symbol_Human")
sum(is.na(Orthologs_Rats_Humans_Symbols$Symbol_Human)==TRUE | is.na(Orthologs_Rats_Humans_Symbols$Symbol_Rat)==TRUE)
#[1] 3746

Orthologs_Rats_Humans_Symbols_NoNA<-Orthologs_Rats_Humans_Symbols[is.na(Orthologs_Rats_Humans_Symbols$Symbol_Human)==FALSE & is.na(Orthologs_Rats_Humans_Symbols$Symbol_Rat)==FALSE,]
str(Orthologs_Rats_Humans_Symbols_NoNA)
# 'data.frame':	16982 obs. of  2 variables:
# $ Symbol_Rat  : chr  "Acadm" "Acadvl" "Acat1" "Acvr1" ...
# $ Symbol_Human: chr  "ACADM" "ACADVL" "ACAT1" "ACVR1" ...

write.csv(Orthologs_Mice_Rats,"Orthologs_Mice_Rats_InformaticsJax_OrthologDatabase_20210228.txt")
write.csv(Orthologs_Rats_Humans, "Orthologs_Rats_Humans_InformaticsJax_OrthologDatabase_20210228.txt")
write.csv(Orthologs_Mice_Humans, "Orthologs_Mice_Humans_InformaticsJax_OrthologDatabase_20210228.txt")

write.csv(Orthologs_Mice_Rats_Symbols_NoNA,"Orthologs_Mice_Rats_Symbols_NoNA_InformaticsJax_OrthologDatabase_20210228.txt")
write.csv(Orthologs_Rats_Humans_Symbols_NoNA, "Orthologs_Rats_Humans_Symbols_NoNA_InformaticsJax_OrthologDatabase_20210228.txt")
write.csv(Orthologs_Mice_Humans_Symbols_NoNA, "Orthologs_Mice_Humans_Symbols_NoNA_InformaticsJax_OrthologDatabase_20210228.txt")


str(DropViz_HC_CellTypeClusterMarkers)
colnames(DropViz_HC_CellTypeClusterMarkers)[1]<-"Symbol_Mouse"
colnames(DropViz_HC_CellTypeClusterMarkers_Concise)[1]<-"Symbol_Mouse"
colnames(DropViz_STR_CellTypeClusterMarkers)[1]<-"Symbol_Mouse"
colnames(DropViz_STR_CellTypeClusterMarkers_Concise)[1]<-"Symbol_Mouse"

str(DropViz_HC_CellTypeClusterMarkers)
#'data.frame':	4780 obs. of  21 variables:
#'
DropViz_HC_CellTypeClusterMarkers_RatOrthologs<-join(DropViz_HC_CellTypeClusterMarkers, Orthologs_Mice_Rats_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(DropViz_HC_CellTypeClusterMarkers_RatOrthologs)
#'data.frame':	4450 obs. of  43 variables:

write.csv(DropViz_HC_CellTypeClusterMarkers_RatOrthologs, "DropViz_HC_CellTypeClusterMarkers_RatOrthologs.csv")

DropViz_HC_CellTypeClusterMarkers_HumanOrthologs<-join(DropViz_HC_CellTypeClusterMarkers, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(DropViz_HC_CellTypeClusterMarkers_HumanOrthologs)
#'data.frame':	4271 obs. of  22 variables:

write.csv(DropViz_HC_CellTypeClusterMarkers_HumanOrthologs, "DropViz_HC_CellTypeClusterMarkers_HumanOrthologs.csv")


DropViz_HC_CellTypeClusterMarkers_Concise_RatOrthologs<-join(DropViz_HC_CellTypeClusterMarkers_Concise, Orthologs_Mice_Rats_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(DropViz_HC_CellTypeClusterMarkers_Concise_RatOrthologs)
#'data.frame':	1597 obs. of  22 variables:

write.csv(DropViz_HC_CellTypeClusterMarkers_Concise_RatOrthologs, "DropViz_HC_CellTypeClusterMarkers_Concise_RatOrthologs.csv")


DropViz_HC_CellTypeClusterMarkers_Concise_HumanOrthologs<-join(DropViz_HC_CellTypeClusterMarkers_Concise, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(DropViz_HC_CellTypeClusterMarkers_Concise_HumanOrthologs)
#'data.frame':	1518 obs. of  22 variables:

write.csv(DropViz_HC_CellTypeClusterMarkers_Concise_HumanOrthologs, "DropViz_HC_CellTypeClusterMarkers_Concise_HumanOrthologs.csv")


DropViz_STR_CellTypeClusterMarkers_RatOrthologs<-join(DropViz_STR_CellTypeClusterMarkers, Orthologs_Mice_Rats_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(DropViz_STR_CellTypeClusterMarkers_RatOrthologs)
#'data.frame':	4729 obs. of  22 variables:

write.csv(DropViz_STR_CellTypeClusterMarkers_RatOrthologs, "DropViz_STR_CellTypeClusterMarkers_RatOrthologs.csv")


DropViz_STR_CellTypeClusterMarkers_HumanOrthologs<-join(DropViz_STR_CellTypeClusterMarkers, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(DropViz_STR_CellTypeClusterMarkers_HumanOrthologs)
#'data.frame':	4476 obs. of  22 variables:

write.csv(DropViz_STR_CellTypeClusterMarkers_HumanOrthologs, "DropViz_STR_CellTypeClusterMarkers_HumanOrthologs.csv")


DropViz_STR_CellTypeClusterMarkers_Concise_RatOrthologs<-join(DropViz_STR_CellTypeClusterMarkers_Concise, Orthologs_Mice_Rats_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(DropViz_STR_CellTypeClusterMarkers_Concise_RatOrthologs)
#'data.frame':	1667 obs. of  22 variables:

write.csv(DropViz_STR_CellTypeClusterMarkers_Concise_RatOrthologs, "DropViz_STR_CellTypeClusterMarkers_Concise_RatOrthologs.csv")


DropViz_STR_CellTypeClusterMarkers_Concise_HumanOrthologs<-join(DropViz_STR_CellTypeClusterMarkers_Concise, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(DropViz_STR_CellTypeClusterMarkers_Concise_HumanOrthologs)
#'data.frame':	1579 obs. of  22 variables:
#'
write.csv(DropViz_STR_CellTypeClusterMarkers_Concise_HumanOrthologs, "DropViz_STR_CellTypeClusterMarkers_Concise_HumanOrthologs.csv")

str(Zeisel_CA1PyramidalIndex)
#'data.frame':	409 obs. of  1 variable:

colnames(Zeisel_CA1PyramidalIndex)[1]<-"Symbol_Mouse_AsUpperCase"
str(Zeisel_CA1PyramidalIndex)
Zeisel_CA1PyramidalIndex$Dataset<-"Zeisel_CA1PyramidalNeurons_HC"
str(Zeisel_CA1PyramidalIndex)
Orthologs_Mice_Rats_Symbols_NoNA$Symbol_Mouse_AsUpperCase<-toupper(Orthologs_Mice_Rats_Symbols_NoNA$Symbol_Mouse)

Zeisel_CA1PyramidalIndex_RatOrthologs<-join(Zeisel_CA1PyramidalIndex, Orthologs_Mice_Rats_Symbols_NoNA, by="Symbol_Mouse_AsUpperCase", type="inner", match="all")
str(Zeisel_CA1PyramidalIndex_RatOrthologs)
#'data.frame':	335 obs. of  4 variables:

write.csv(Zeisel_CA1PyramidalIndex_RatOrthologs, "Zeisel_CA1PyramidalIndex_RatOrthologs.csv")

Orthologs_Mice_Humans_Symbols_NoNA$Symbol_Mouse_AsUpperCase<-toupper(Orthologs_Mice_Humans_Symbols_NoNA$Symbol_Mouse)

Zeisel_CA1PyramidalIndex_HumanOrthologs<-join(Zeisel_CA1PyramidalIndex, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Mouse_AsUpperCase", type="inner", match="all")
str(Zeisel_CA1PyramidalIndex_HumanOrthologs)
#'data.frame':	333 obs. of  4 variables:
write.csv(Zeisel_CA1PyramidalIndex_HumanOrthologs, "Zeisel_CA1PyramidalIndex_HumanOrthologs.csv")


str(StressGeneSets)
"'data.frame':	2469 obs. of  2 variables:"

colnames(StressGeneSets)<-"Symbol_Mouse"


StressGeneSets_RatOrthologs<-join(StressGeneSets, Orthologs_Mice_Rats_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(StressGeneSets_RatOrthologs)
#'data.frame':	1998 obs. of  3 variables:
#That's a bigger hit - many of these symbols must be old.

write.csv(StressGeneSets_RatOrthologs, "StressGeneSets_RatOrthologs.csv")

StressGeneSets_HumanOrthologs<-join(StressGeneSets, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(StressGeneSets_HumanOrthologs)
#'data.frame':	1919 obs. of  3 variables:
#'#Ouch! almost half of the gene symbols from the original datasets don't have human orthologs. Gah.

write.csv(StressGeneSets_HumanOrthologs, "StressGeneSets_HumanOrthologs.csv")

StressGeneSets_MouseOrthologs_Symbols_NoNA<-unique(StressGeneSets[,c(1,2)])
str(StressGeneSets_MouseOrthologs_Symbols_NoNA)
#'data.frame':	2469 obs. of  2 variables:
sum(is.na(StressGeneSets_MouseOrthologs_Symbols_NoNA$Symbol_Mouse))
#[1] 0
write.csv(StressGeneSets_MouseOrthologs_Symbols_NoNA, "StressGeneSets_MouseOrthologs_Symbols_NoNA.csv")

str(StressGeneSets_RatOrthologs)
#'data.frame':	1998 obs. of  4 variables:
StressGeneSets_RatOrthologs_Symbols_NoNA<-unique(StressGeneSets_RatOrthologs[,c(3,2)])
str(StressGeneSets_RatOrthologs_Symbols_NoNA)
#'data.frame':	1996 obs. of  2 variables:
sum(is.na(StressGeneSets_RatOrthologs_Symbols_NoNA$Symbol_Rat))
#[1] 0
write.csv(StressGeneSets_RatOrthologs_Symbols_NoNA, "StressGeneSets_RatOrthologs_Symbols_NoNA.csv")

str(StressGeneSets_HumanOrthologs)
#'data.frame':	1919 obs. of  4 variables:
StressGeneSets_HumanOrthologs_Symbols_NoNA<-unique(StressGeneSets_HumanOrthologs[,c(3,2)])
str(StressGeneSets_HumanOrthologs_Symbols_NoNA)
#'data.frame':	1915 obs. of  2 variables:
sum(is.na(StressGeneSets_HumanOrthologs_Symbols_NoNA$Symbol_Human))
#[1] 0
write.csv(StressGeneSets_HumanOrthologs_Symbols_NoNA, "StressGeneSets_HumanOrthologs_Symbols_NoNA.csv")

str(HumanPsychiatricGeneSets)
#'data.frame':	4416 obs. of  2 variables:
colnames(HumanPsychiatricGeneSets)[1]<-"Symbol_Human"

HumanPsychiatricGeneSets_MouseOrthologs<-join(HumanPsychiatricGeneSets, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Human", type="inner", match="all")
str(HumanPsychiatricGeneSets_MouseOrthologs)
#'data.frame':	3971 obs. of  3 variables:

write.csv(HumanPsychiatricGeneSets_MouseOrthologs, "HumanPsychiatricGeneSets_MouseOrthologs.csv")

HumanPsychiatricGeneSets_RatOrthologs<-join(HumanPsychiatricGeneSets, Orthologs_Rats_Humans_Symbols_NoNA, by="Symbol_Human", type="inner", match="all")
str(HumanPsychiatricGeneSets_RatOrthologs)
#'data.frame':	3866 obs. of  3 variables:

write.csv(HumanPsychiatricGeneSets_RatOrthologs, "HumanPsychiatricGeneSets_RatOrthologs.csv")

str(HumanPsychiatricGeneSets)
#'data.frame':	4416 obs. of  2 variables:
HumanPsychiatricGeneSets_HumanOrthologs_JustSymbols_NoNA<-unique(HumanPsychiatricGeneSets[,c(1:2)])
str(HumanPsychiatricGeneSets_HumanOrthologs_JustSymbols_NoNA)
#'data.frame':	4416 obs. of  2 variables:
sum(is.na(HumanPsychiatricGeneSets_HumanOrthologs_JustSymbols_NoNA$Symbol_Human))
#[1] 0
write.csv(HumanPsychiatricGeneSets_HumanOrthologs_JustSymbols_NoNA, "HumanPsychiatricGeneSets_HumanOrthologs_JustSymbols_NoNA.csv")

str(HumanPsychiatricGeneSets_MouseOrthologs)
#'data.frame':	3971 obs. of  3 variables:
HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA<-unique(HumanPsychiatricGeneSets_MouseOrthologs[,c(3:2)])
str(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA)
#'data.frame':	3970 obs. of  2 variables:
sum(is.na(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$Symbol_Mouse))
#[1] 0
write.csv(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA, "HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA.csv")

str(HumanPsychiatricGeneSets_RatOrthologs)
#'data.frame':	3866 obs. of  3 variables:
HumanPsychiatricGeneSets_RatOrthologs_JustSymbols_NoNA<-unique(HumanPsychiatricGeneSets_RatOrthologs[,c(3:2)])
str(HumanPsychiatricGeneSets_RatOrthologs_JustSymbols_NoNA)
#'data.frame':	3863 obs. of  2 variables:
sum(is.na(HumanPsychiatricGeneSets_RatOrthologs_JustSymbols_NoNA$Symbol_Rat))
#[1] 0
write.csv(HumanPsychiatricGeneSets_RatOrthologs_JustSymbols_NoNA, "HumanPsychiatricGeneSets_RatOrthologs_JustSymbols_NoNA.csv")



colnames(BrainInABlenderDatabase_Human)[4]<-"Symbol_Human"
colnames(BrainInABlenderDatabase_Human)[5]<-"Symbol_Mouse"
colnames(BrainInABlenderDatabase_Mouse)[4]<-"Symbol_Human"
colnames(BrainInABlenderDatabase_Mouse)[5]<-"Symbol_Mouse"

BrainInABlenderDatabase_Mouse_RatOrthologs<-join(BrainInABlenderDatabase_Mouse, Orthologs_Mice_Rats_Symbols_NoNA, by="Symbol_Mouse", type="inner", match="all")
str(BrainInABlenderDatabase_Mouse_RatOrthologs)

BrainInABlenderDatabase_Human_RatOrthologs<-join(BrainInABlenderDatabase_Human, Orthologs_Rats_Humans_Symbols_NoNA, by="Symbol_Human", type="inner", match="all")
str(BrainInABlenderDatabase_Human_RatOrthologs)

write.csv(BrainInABlenderDatabase_Human_RatOrthologs, "BrainInABlenderDatabase_Human_RatOrthologs.csv")
write.csv(BrainInABlenderDatabase_Mouse_RatOrthologs, "BrainInABlenderDatabase_Mouse_RatOrthologs.csv")

#hmm... I'm pretty sure that still includes NAs and duplicates though. Let's simplify it down:

BrainInABlenderDatabase_Human_JustSymbols<-BrainInABlenderDatabase_Human[,c(4,5,13)] 
BrainInABlenderDatabase_Mouse_JustSymbols<-BrainInABlenderDatabase_Mouse[,c(4,5,13)] 

BrainInABlenderDatabase_Human_RatOrthologs_JustSymbols<-BrainInABlenderDatabase_Human_RatOrthologs[,c(15,13)] 
BrainInABlenderDatabase_Mouse_RatOrthologs_JustSymbols<-BrainInABlenderDatabase_Mouse_RatOrthologs[,c(15,13)] 


str(BrainInABlenderDatabase_Mouse_JustSymbols)
#'data.frame':	3240 obs. of  3 variables:
BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA<-unique(BrainInABlenderDatabase_Mouse_JustSymbols[,c(1,3)])
head(BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA)
str(BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA)
#'data.frame':	2853 obs. of  2 variables:
sum(is.na(BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA$Symbol_Human))
#[1] 29
BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA<-BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA[is.na(BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA$Symbol_Human)==FALSE,]
str(BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA)
#'data.frame':	2824 obs. of  2 variables:

str(BrainInABlenderDatabase_Mouse_JustSymbols)
#'data.frame':	3240 obs. of  3 variables:
BrainInABlenderDatabase_Mouse_JustSymbols_MouseOrthologs_NoNA<-unique(BrainInABlenderDatabase_Mouse_JustSymbols[,c(2,3)])
head(BrainInABlenderDatabase_Mouse_JustSymbols_MouseOrthologs_NoNA)
str(BrainInABlenderDatabase_Mouse_JustSymbols_MouseOrthologs_NoNA)
#'data.frame':	3218 obs. of  2 variables:
sum(is.na(BrainInABlenderDatabase_Mouse_JustSymbols_MouseOrthologs_NoNA$Symbol_Mouse))
#[1] 0

str(BrainInABlenderDatabase_Mouse_RatOrthologs_JustSymbols)
#'data.frame':	2735 obs. of  2 variables:
BrainInABlenderDatabase_Mouse_RatOrthologs_JustSymbols_RatOrthologs_NoNA<-unique(BrainInABlenderDatabase_Mouse_RatOrthologs_JustSymbols[,c(1,2)])
head(BrainInABlenderDatabase_Mouse_RatOrthologs_JustSymbols_RatOrthologs_NoNA)
str(BrainInABlenderDatabase_Mouse_RatOrthologs_JustSymbols_RatOrthologs_NoNA)
#'data.frame':	2690 obs. of  2 variables:
sum(is.na(BrainInABlenderDatabase_Mouse_RatOrthologs_JustSymbols_RatOrthologs_NoNA$Symbol_Rat))
#[1] 0

str(BrainInABlenderDatabase_Human_JustSymbols)
#'data.frame':	143 obs. of  3 variables:
BrainInABlenderDatabase_Human_JustSymbols_HumanOrthologs_NoNA<-unique(BrainInABlenderDatabase_Human_JustSymbols[,c(1,3)])
head(BrainInABlenderDatabase_Human_JustSymbols_HumanOrthologs_NoNA)
str(BrainInABlenderDatabase_Human_JustSymbols_HumanOrthologs_NoNA)
#'data.frame':	143 obs. of  2 variables:
sum(is.na(BrainInABlenderDatabase_Human_JustSymbols_HumanOrthologs_NoNA$Symbol_Human))
#[1] 0

str(BrainInABlenderDatabase_Human_JustSymbols)
#'data.frame':	143 obs. of  3 variables:
BrainInABlenderDatabase_Human_JustSymbols_MouseOrthologs_NoNA<-unique(BrainInABlenderDatabase_Human_JustSymbols[,c(2,3)])
head(BrainInABlenderDatabase_Human_JustSymbols_MouseOrthologs_NoNA)
str(BrainInABlenderDatabase_Human_JustSymbols_MouseOrthologs_NoNA)
#'data.frame':	139 obs. of  2 variables:
sum(is.na(BrainInABlenderDatabase_Human_JustSymbols_MouseOrthologs_NoNA$Symbol_Mouse))
#[1] 5
BrainInABlenderDatabase_Human_JustSymbols_MouseOrthologs_NoNA<-BrainInABlenderDatabase_Human_JustSymbols_MouseOrthologs_NoNA[is.na(BrainInABlenderDatabase_Human_JustSymbols_MouseOrthologs_NoNA$Symbol_Mouse)==FALSE,]
str(BrainInABlenderDatabase_Human_JustSymbols_MouseOrthologs_NoNA)
#'data.frame':	134 obs. of  2 variables:

str(BrainInABlenderDatabase_Human_RatOrthologs_JustSymbols)
#'data.frame':	129 obs. of  2 variables:
BrainInABlenderDatabase_Human_RatOrthologs_JustSymbols_RatOrthologs_NoNA<-unique(BrainInABlenderDatabase_Human_RatOrthologs_JustSymbols[,c(1,2)])
head(BrainInABlenderDatabase_Human_RatOrthologs_JustSymbols_RatOrthologs_NoNA)
str(BrainInABlenderDatabase_Human_RatOrthologs_JustSymbols_RatOrthologs_NoNA)
#'data.frame':	128 obs. of  2 variables:
sum(is.na(BrainInABlenderDatabase_Human_RatOrthologs_JustSymbols_RatOrthologs_NoNA$Symbol_Rat))
#[1] 0

#Can I just combine these together now to make life easier?

str(BrainInABlenderDatabase_Human_JustSymbols_HumanOrthologs_NoNA)
str(BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA)

BrainInABlenderDatabase_All_JustSymbols_HumanOrthologs_NoNA<-rbind.data.frame(BrainInABlenderDatabase_Mouse_JustSymbols_HumanOrthologs_NoNA, BrainInABlenderDatabase_Human_JustSymbols_HumanOrthologs_NoNA)

BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA<-rbind.data.frame(BrainInABlenderDatabase_Mouse_JustSymbols_MouseOrthologs_NoNA, BrainInABlenderDatabase_Human_JustSymbols_MouseOrthologs_NoNA)

BrainInABlenderDatabase_All_JustSymbols_RatOrthologs_NoNA<-rbind.data.frame(BrainInABlenderDatabase_Mouse_RatOrthologs_JustSymbols_RatOrthologs_NoNA, BrainInABlenderDatabase_Human_RatOrthologs_JustSymbols_RatOrthologs_NoNA)

write.csv(BrainInABlenderDatabase_All_JustSymbols_HumanOrthologs_NoNA, "BrainInABlenderDatabase_All_JustSymbols_HumanOrthologs_NoNA.csv")
write.csv(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA, "BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA.csv")
write.csv(BrainInABlenderDatabase_All_JustSymbols_RatOrthologs_NoNA, "BrainInABlenderDatabase_All_JustSymbols_RatOrthologs_NoNA.csv")



str(InternalizingBehaviorGeneSets)

str(InternalizingBehaviorGeneSets)
# 'data.frame':	3219 obs. of  2 variables:
#   $ GeneSymbol_Rat: chr  "Abhd10" "Acbd4" "Acss2" "Adam19" ...
# $ DatasetName   : chr  "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" ...

colnames(InternalizingBehaviorGeneSets)[1]<-"Symbol_Rat"

InternalizingBehaviorGeneSets_MouseOrthologs<-join(InternalizingBehaviorGeneSets, Orthologs_Mice_Rats_Symbols_NoNA, by="Symbol_Rat", type="inner", match="all")
str(InternalizingBehaviorGeneSets_MouseOrthologs)
# 'data.frame':	2815 obs. of  4 variables:
#   $ Symbol_Rat              : chr  "Abhd10" "Acbd4" "Acss2" "Adam19" ...
# $ DatasetName             : chr  "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" ...
# $ Symbol_Mouse            : chr  "Abhd10" "Acbd4" "Acss2" "Adam19" ...
# $ Symbol_Mouse_AsUpperCase: chr  "ABHD10" "ACBD4" "ACSS2" "ADAM19" ...

write.csv(InternalizingBehaviorGeneSets_MouseOrthologs, "InternalizingBehaviorGeneSets_MouseOrthologs.csv")

InternalizingBehaviorGeneSets_HumanOrthologs<-join(InternalizingBehaviorGeneSets, Orthologs_Rats_Humans_Symbols_NoNA, by="Symbol_Rat", type="inner", match="all")
str(InternalizingBehaviorGeneSets_HumanOrthologs)
# 'data.frame':	2690 obs. of  3 variables:
#   $ Symbol_Rat  : chr  "Abhd10" "Acbd4" "Acss2" "Adam19" ...
# $ DatasetName : chr  "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" ...
# $ Symbol_Human: chr  "ABHD10" "ACBD4" "ACSS2" "ADAM19" ...

write.csv(InternalizingBehaviorGeneSets_HumanOrthologs, "InternalizingBehaviorGeneSets_HumanOrthologs.csv")

#Let's make a version with duplicates and NAs removed:

str(InternalizingBehaviorGeneSets)
#'data.frame':	3219 obs. of  2 variables:
InternalizingBehaviorGeneSets_RatOrthologs_NoNA<-unique(InternalizingBehaviorGeneSets[,c(1:2)])
str(InternalizingBehaviorGeneSets_RatOrthologs_NoNA)
#'data.frame':	3219 obs. of  2 variables:
InternalizingBehaviorGeneSets_RatOrthologs_NoNA<-InternalizingBehaviorGeneSets_RatOrthologs_NoNA[is.na(InternalizingBehaviorGeneSets_RatOrthologs_NoNA$Symbol_Rat)==FALSE,]
str(InternalizingBehaviorGeneSets_RatOrthologs_NoNA)
#'data.frame':	3219 obs. of  2 variables:

str(InternalizingBehaviorGeneSets_MouseOrthologs)
#'data.frame':	2815 obs. of  4 variables:
InternalizingBehaviorGeneSets_MouseOrthologs_NoNA<-unique(InternalizingBehaviorGeneSets_MouseOrthologs[,c(3,2)])
str(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA)
#'data.frame':	2815 obs. of  2 variables:
InternalizingBehaviorGeneSets_MouseOrthologs_NoNA<-InternalizingBehaviorGeneSets_MouseOrthologs_NoNA[is.na(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$Symbol_Mouse)==FALSE,]
str(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA)
#'data.frame':	2815 obs. of  2 variables:

str(InternalizingBehaviorGeneSets_HumanOrthologs)
#'data.frame':	2690 obs. of  3 variables:
InternalizingBehaviorGeneSets_HumanOrthologs_NoNA<-unique(InternalizingBehaviorGeneSets_HumanOrthologs[,c(3,2)])
str(InternalizingBehaviorGeneSets_HumanOrthologs_NoNA)
#'data.frame':	2690 obs. of  2 variables:
InternalizingBehaviorGeneSets_HumanOrthologs_NoNA<-InternalizingBehaviorGeneSets_HumanOrthologs_NoNA[is.na(InternalizingBehaviorGeneSets_HumanOrthologs_NoNA$Symbol_Human)==FALSE,]
str(InternalizingBehaviorGeneSets_HumanOrthologs_NoNA)
#'data.frame':	2690 obs. of  2 variables:


#Ok, now tackling the GMTs.
#This is going to be substantially more of a pain in the butt.
#Probably easiest to transpose the data.frames, run a for loop to replace gene symbols with orthologs and trim out missing values for each column, and transpose back.

#Old Code:

# str(c2.all.v7.3.symbols)
# str(t(c2.all.v7.3.symbols))
# #seems to transpose into a character matrix? I guess I could have used that shortcut earlier, LOL.
# 
# #Making some matrices to hold my results:
# 
# TempMatrix<-t(c2.all.v7.3.symbols)
# TempMatrix_MouseOrtholog<-matrix(data=as.character(), ncol=ncol(TempMatrix), nrow=(nrow(TempMatrix)+100))
# TempMatrix_RatOrtholog<-matrix(data=as.character(), ncol=ncol(TempMatrix), nrow=(nrow(TempMatrix)+100))
# #I had to add some rows because not all of the orthologs are a 1:1 match.
# 
# for(i in c(3:ncol(TempMatrix))){
#   
#   TempMatrixColumn<-data.frame("Symbol_Human"=TempMatrix[,i])
#   
#   TempColumnOrthologRats<-join(TempMatrixColumn, Orthologs_Rats_Humans_Symbols_NoNA, by="Symbol_Human", type="inner", match="all")
#   if(length(unique(TempColumnOrthologRats$Symbol_Rat))>1){
#   TempMatrix_RatOrtholog[c(3:(length(unique(TempColumnOrthologRats$Symbol_Rat))+2)),i]<-unique(TempColumnOrthologRats$Symbol_Rat)}
#   else{}
#   
#   TempColumnOrthologMice<-join(TempMatrixColumn, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Human", type="inner", match="all")
#   if(length(unique(TempColumnOrthologMice$Symbol_Mouse))>1){
#   TempMatrix_MouseOrtholog[c(3:(length(unique(TempColumnOrthologMice$Symbol_Mouse))+2)),i]<-unique(TempColumnOrthologMice$Symbol_Mouse)}
#   else{}
#   
#   rm(TempMatrixColumn, TempColumnOrthologRats,TempColumnOrthologMice)
#   
# }

# 
# TempMatrix_MouseOrtholog[c(1:2),]<-TempMatrix[c(1:2),]
# TempMatrix_RatOrtholog[c(1:2),]<-TempMatrix[c(1:2),]
# 
#   c2.all.v7.3.symbols_MouseOrtholog<-t(TempMatrix_MouseOrtholog)
#   c2.all.v7.3.symbols_RatOrtholog<-t(TempMatrix_RatOrtholog)
#   
#   write.table(c2.all.v7.3.symbols_MouseOrtholog, "c2.all.v7.3.symbols_MouseOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
#   write.table(c2.all.v7.3.symbols_RatOrtholog, "c2.all.v7.3.symbols_RatOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
#   
#   #I'm having trouble opening those files.
#   
#   write.csv(c2.all.v7.3.symbols_MouseOrtholog, "c2.all.v7.3.symbols_MouseOrtholog.csv", na="")
#   #That didn't work either - still can't open the file.
#   #let's test whether it is just a problem with writing out:
#   write.csv(TempMatrix, "TempMatrix.csv")
#   #That didn't work either.
#   #Let's try something that wrote out fine earlier.
#   write.csv(HumanPsychiatricGeneSets_RatOrthologs, "TestingWritingOut.csv")
#   #That also isn't working. I think it may be time to save the workspace and code and restart the computer. :(
#   #restarted the computer and now everything is opening up fine. :)


#New Code:


setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")

FindOrthologsForGMTs<-function(GMTFileAsCharacterMatrix){
  
  TempMatrix<-t(GMTFileAsCharacterMatrix)
  TempMatrix_MouseOrtholog<-matrix(data=as.character(), ncol=ncol(TempMatrix), nrow=(nrow(TempMatrix)+100))
  TempMatrix_RatOrtholog<-matrix(data=as.character(), ncol=ncol(TempMatrix), nrow=(nrow(TempMatrix)+100))
  #I had to add some rows because not all of the orthologs are a 1:1 match.
  
  for(i in c(1:ncol(TempMatrix))){
    
    TempMatrixColumn<-data.frame("Symbol_Human"=TempMatrix[,i])
    
    TempColumnOrthologRats<-join(TempMatrixColumn, Orthologs_Rats_Humans_Symbols_NoNA, by="Symbol_Human", type="inner", match="all")
    if(length(unique(TempColumnOrthologRats$Symbol_Rat))>1){
      TempMatrix_RatOrtholog[c(3:(length(unique(TempColumnOrthologRats$Symbol_Rat))+2)),i]<-unique(TempColumnOrthologRats$Symbol_Rat)}
    else{}
    
    TempColumnOrthologMice<-join(TempMatrixColumn, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Human", type="inner", match="all")
    if(length(unique(TempColumnOrthologMice$Symbol_Mouse))>1){
      TempMatrix_MouseOrtholog[c(3:(length(unique(TempColumnOrthologMice$Symbol_Mouse))+2)),i]<-unique(TempColumnOrthologMice$Symbol_Mouse)}
    else{}
    
    rm(TempMatrixColumn, TempColumnOrthologRats,TempColumnOrthologMice)
    
  }
  TempMatrix_MouseOrtholog[c(1:2),]<-TempMatrix[c(1:2),]
  TempMatrix_RatOrtholog[c(1:2),]<-TempMatrix[c(1:2),]
  return(list(t(TempMatrix_MouseOrtholog), t(TempMatrix_RatOrtholog)))
}


Temp<-FindOrthologsForGMTs(h.all.v7.3.symbols)
str(Temp)
# List of 2
# $ : chr [1:50, 1:302] "HALLMARK_TNFA_SIGNALING_VIA_NFKB" "HALLMARK_HYPOXIA" "HALLMARK_CHOLESTEROL_HOMEOSTASIS" "HALLMARK_MITOTIC_SPINDLE" ...
# $ : chr [1:50, 1:302] "HALLMARK_TNFA_SIGNALING_VIA_NFKB" "HALLMARK_HYPOXIA" "HALLMARK_CHOLESTEROL_HOMEOSTASIS" "HALLMARK_MITOTIC_SPINDLE" ...

h.all.v7.3.symbols_MouseOrtholog<-Temp[[1]]
h.all.v7.3.symbols_RatOrtholog<-Temp[[2]]

write.table(h.all.v7.3.symbols_MouseOrtholog, "h.all.v7.3.symbols_MouseOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
write.table(h.all.v7.3.symbols_RatOrtholog, "h.all.v7.3.symbols_RatOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

rm(Temp)



  Temp<-FindOrthologsForGMTs(c5.all.v7.3.symbols)
  str(Temp)

  c5.all.v7.3.symbols_MouseOrtholog<-Temp[[1]]
  c5.all.v7.3.symbols_RatOrtholog<-Temp[[2]]
  # List of 2
  # $ : chr [1:14996, 1:2084] "GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE" "GOBP_REPRODUCTION" "GOBP_SINGLE_STRAND_BREAK_REPAIR" "GOBP_REGULATION_OF_DNA_RECOMBINATION" ...
  # $ : chr [1:14996, 1:2084] "GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE" "GOBP_REPRODUCTION" "GOBP_SINGLE_STRAND_BREAK_REPAIR" "GOBP_REGULATION_OF_DNA_RECOMBINATION" ...
  
  write.table(c5.all.v7.3.symbols_MouseOrtholog, "c5.all.v7.3.symbols_MouseOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
  write.table(c5.all.v7.3.symbols_RatOrtholog, "c5.all.v7.3.symbols_RatOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
  
 rm(Temp)
  
  
 Temp<-FindOrthologsForGMTs(c8.all.v7.3.symbols)
  str(Temp)
  List of 2
  # $ : chr [1:673, 1:1872] "BUSSLINGER_ESOPHAGEAL_QUIESCENT_BASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_PROLIFERATING_BASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_EARLY_SUPRABASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_LATE_SUPRABASAL_CELLS" ...
  # $ : chr [1:673, 1:1872] "BUSSLINGER_ESOPHAGEAL_QUIESCENT_BASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_PROLIFERATING_BASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_EARLY_SUPRABASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_LATE_SUPRABASAL_CELLS" ...
  
  c8.all.v7.3.symbols_MouseOrtholog<-Temp[[1]]
  c8.all.v7.3.symbols_RatOrtholog<-Temp[[2]]
  
  write.table(c8.all.v7.3.symbols_MouseOrtholog, "c8.all.v7.3.symbols_MouseOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
  write.table(c8.all.v7.3.symbols_RatOrtholog, "c8.all.v7.3.symbols_RatOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
  
  rm(Temp)
  
 
  Temp<-FindOrthologsForGMTs(c3.all.v7.3.symbols)
  str(Temp)
  # List of 2
  # $ : chr [1:3731, 1:2102] "MIR153_5P" "MIR8485" "MIR3662" "MIR607" ...
  # $ : chr [1:3731, 1:2102] "MIR153_5P" "MIR8485" "MIR3662" "MIR607" ...
  
 c3.all.v7.3.symbols_MouseOrtholog<-Temp[[1]]
 c3.all.v7.3.symbols_RatOrtholog<-Temp[[2]]
 
write.table(c3.all.v7.3.symbols_MouseOrtholog, "c3.all.v7.3.symbols_MouseOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
write.table(c3.all.v7.3.symbols_RatOrtholog, "c3.all.v7.3.symbols_RatOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
 
rm(Temp)

Temp<-FindOrthologsForGMTs(c7.all.v7.3.symbols)
str(Temp)
# List of 2
# $ : chr [1:5219, 1:2094] "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP" "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_DN" "KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_UP" "KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_DN" ...
# $ : chr [1:5219, 1:2094] "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP" "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_DN" "KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_UP" "KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_DN" ...

c7.all.v7.3.symbols_MouseOrtholog<-Temp[[1]]
c7.all.v7.3.symbols_RatOrtholog<-Temp[[2]]

write.table(c7.all.v7.3.symbols_MouseOrtholog, "c7.all.v7.3.symbols_MouseOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
write.table(c7.all.v7.3.symbols_RatOrtholog, "c7.all.v7.3.symbols_RatOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

rm(Temp)


Temp<-FindOrthologsForGMTs(c2.all.v7.3.symbols)
str(Temp)
# List of 2
# $ : chr [1:6255, 1:2043] "BENITEZ_GBM_PROTEASOME_INHIBITION_RESPONSE" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_UP" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_DN" "BLANCO_MELO_MERS_COV_INFECTION_MCR5_CELLS_UP" ...
# $ : chr [1:6255, 1:2043] "BENITEZ_GBM_PROTEASOME_INHIBITION_RESPONSE" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_UP" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_DN" "BLANCO_MELO_MERS_COV_INFECTION_MCR5_CELLS_UP" ...

c2.all.v7.3.symbols_MouseOrtholog<-Temp[[1]]
c2.all.v7.3.symbols_RatOrtholog<-Temp[[2]]

write.table(c2.all.v7.3.symbols_MouseOrtholog, "c2.all.v7.3.symbols_MouseOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
write.table(c2.all.v7.3.symbols_RatOrtholog, "c2.all.v7.3.symbols_RatOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

rm(Temp)



CombinedHCspecific
CombinedHCspecific[,1]
#This one was originally based on a mixture of mouse and human data - ouch
#Park=Mouse coexpression data (but in all caps)
#Johnson=Human coexpression data
#Cembrowski=mouse (but in all caps)

CombinedHCspecific_Mouse<-CombinedHCspecific[c(1:14,39:69),]
CombinedHCspecific_Human<-CombinedHCspecific[c(15:38),]

TempMatrix<-t(CombinedHCspecific_Human)
#This should with the old code above:
TempMatrix_MouseOrtholog<-matrix(data=as.character(), ncol=ncol(TempMatrix), nrow=(nrow(TempMatrix)+100))
TempMatrix_RatOrtholog<-matrix(data=as.character(), ncol=ncol(TempMatrix), nrow=(nrow(TempMatrix)+100))
# #I had to add some rows because not all of the orthologs are a 1:1 match.
# 
for(i in c(1:ncol(TempMatrix))){

  TempMatrixColumn<-data.frame("Symbol_Human"=TempMatrix[,i])

  TempColumnOrthologRats<-join(TempMatrixColumn, Orthologs_Rats_Humans_Symbols_NoNA, by="Symbol_Human", type="inner", match="all")
  if(length(unique(TempColumnOrthologRats$Symbol_Rat))>1){
  TempMatrix_RatOrtholog[c(3:(length(unique(TempColumnOrthologRats$Symbol_Rat))+2)),i]<-unique(TempColumnOrthologRats$Symbol_Rat)}
  else{}

  TempColumnOrthologMice<-join(TempMatrixColumn, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Human", type="inner", match="all")
  if(length(unique(TempColumnOrthologMice$Symbol_Mouse))>1){
  TempMatrix_MouseOrtholog[c(3:(length(unique(TempColumnOrthologMice$Symbol_Mouse))+2)),i]<-unique(TempColumnOrthologMice$Symbol_Mouse)}
  else{}

  rm(TempMatrixColumn, TempColumnOrthologRats,TempColumnOrthologMice)

}


TempMatrix_MouseOrtholog[c(1:2),]<-TempMatrix[c(1:2),]
TempMatrix_RatOrtholog[c(1:2),]<-TempMatrix[c(1:2),]

CombinedHCspecific_Human_MouseOrtholog<-t(TempMatrix_MouseOrtholog)
CombinedHCspecific_Human_RatOrtholog<-t(TempMatrix_RatOrtholog)

write.table(CombinedHCspecific_Human_MouseOrtholog, "CombinedHCspecific_Human_MouseOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
write.table(CombinedHCspecific_Human_RatOrtholog, "CombinedHCspecific_Human_RatOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

write.table(CombinedHCspecific_Human, "CombinedHCspecific_Human.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

TempMatrix<-t(CombinedHCspecific_Mouse)

TempMatrix_HumanOrtholog<-matrix(data=as.character(), ncol=ncol(TempMatrix), nrow=(nrow(TempMatrix)+100))
TempMatrix_RatOrtholog<-matrix(data=as.character(), ncol=ncol(TempMatrix), nrow=(nrow(TempMatrix)+100))
#I had to add some rows because not all of the orthologs are a 1:1 match.

for(i in c(1:ncol(TempMatrix))){
  
  TempMatrixColumn<-data.frame("Symbol_Mouse_AsUpperCase"=TempMatrix[,i])
  
  TempColumnOrthologRats<-join(TempMatrixColumn, Orthologs_Mice_Rats_Symbols_NoNA, by="Symbol_Mouse_AsUpperCase", type="inner", match="all")
  if(length(unique(TempColumnOrthologRats$Symbol_Rat))>1){
    TempMatrix_RatOrtholog[c(3:(length(unique(TempColumnOrthologRats$Symbol_Rat))+2)),i]<-unique(TempColumnOrthologRats$Symbol_Rat)}
  else{}
  
  TempColumnOrthologHuman<-join(TempMatrixColumn, Orthologs_Mice_Humans_Symbols_NoNA, by="Symbol_Mouse_AsUpperCase", type="inner", match="all")
  if(length(unique(TempColumnOrthologHuman$Symbol_Human))>1){
    TempMatrix_HumanOrtholog[c(3:(length(unique(TempColumnOrthologHuman$Symbol_Human))+2)),i]<-unique(TempColumnOrthologHuman$Symbol_Human)}
  else{}
  
  rm(TempMatrixColumn, TempColumnOrthologRats,TempColumnOrthologHuman)
  
}

TempMatrix_HumanOrtholog[c(1:2),]<-TempMatrix[c(1:2),]
TempMatrix_RatOrtholog[c(1:2),]<-TempMatrix[c(1:2),]

CombinedHCspecific_Mouse_HumanOrtholog<-t(TempMatrix_HumanOrtholog)
CombinedHCspecific_Mouse_RatOrtholog<-t(TempMatrix_RatOrtholog)

write.table(CombinedHCspecific_Mouse_HumanOrtholog, "CombinedHCspecific_Mouse_HumanOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
write.table(CombinedHCspecific_Mouse_RatOrtholog, "CombinedHCspecific_Mouse_RatOrtholog.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

write.table(CombinedHCspecific_Mouse, "CombinedHCspecific_Mouse.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

###################################

#Test run before going any further with things (just to make sure I didn't break anything essential while messing around with the .gmt files):

library(fgsea)


#First, for comparison, the old gmt:
setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_toMakeRat")

OldGMT_ForRats<-gmtPathways("c5.all.v7.3.symbols.gmt.txt")
str(OldGMT_ForRats)
# List of 14996
# $ GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE                                                                                                                                                       : chr [1:21] "AKT3" "PPARGC1A" "POLG2" "PARP1" ...
# $ GOBP_REPRODUCTION                                                                                                                                                                           : chr [1:1442] "ADA" "GNPDA1" "ZGLP1" "SYCE1L" ...

names(OldGMT_ForRats)
# [2] "GOBP_REPRODUCTION"                                                                                                                                   
# [3] "GOBP_SINGLE_STRAND_BREAK_REPAIR"    
OldGMT_ForRats[[1]]
# [1] "AKT3"     "PPARGC1A" "POLG2"    "PARP1"    "DNA2"     "TYMP"     "PRIMPOL"  "STOX1"    "SLC25A4"  "LIG3"     "MEF2A"    "MPV17"    "OPA1"     "SLC25A36" "TOP3A"   
# [16] "TP53"     "PIF1"     "SESN2"    "SLC25A33" "MGME1"    "LONP1"  

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")
GMT_ForRats<-gmtPathways("c5.all.v7.3.symbols_RatOrtholog.gmt.txt")
str(GMT_ForRats)
# List of 14996
# $ GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE                                                                                                                                                       : chr [1:2081] "Akt3" "Ppargc1a" "Polg2" "Parp1" ...
# $ GOBP_REPRODUCTION                                                                                                                                                                           : chr [1:2081] "Ada" "Gnpda1" "Zglp1" "Syce1l" ...
# $ GOBP_SINGLE_STRAND_BREAK_REPAIR                                                                                                                                                             : chr [1:2081] "Ercc8" "Parp1" "Aplf" "Ercc6" ...
# $ GOBP_REGULATION_OF_DNA_RECOMBINATION                                                                                                                                                        : chr [1:2081] "Parp3" "Actr2" "Rad50" "Alyref" ...
names(GMT_ForRats)
GMT_ForRats[[1]]
#That definitely is not reading in properly. I tried opening it up in textedit and it doesn't look right either. The delimiter must be incorrect.
#Although my code should have outputted it as tab-delimited text (which is the correct format).
#Ah - I needed quote=FALSE in the output code.  Now it looks right.
#before going back and re-outputting everything though, let's see if the version with quote=FALSE now works.

#Extracting out the results for each comparison:

NACC_EESD_Combined<-read.csv("NACC_EESD_Combined_DateGenesFixed.csv", header=TRUE, stringsAsFactors = FALSE)

colnames(NACC_EESD_Combined)

NAcc_EE_Betas_forGSEA<-NACC_EESD_Combined$Coef.EnrichmentEE_MainEffectsModel
str(NAcc_EE_Betas_forGSEA)
#names(NAcc_EE_Betas_forGSEA)<-toupper(NACC_EESD_Combined$gene_symbol)
names(NAcc_EE_Betas_forGSEA)<-NACC_EESD_Combined$gene_symbol
sum(duplicated(names(NAcc_EE_Betas_forGSEA)))
#[1] 908
sum(is.na(names(NAcc_EE_Betas_forGSEA)))
#[1] 766
length(NAcc_EE_Betas_forGSEA)
#[1] 17765

NAcc_EE_Betas_forGSEA<-NAcc_EE_Betas_forGSEA[is.na(names(NAcc_EE_Betas_forGSEA))==F]
length(NAcc_EE_Betas_forGSEA)
#[1] 16999
sum(is.na(names(NAcc_EE_Betas_forGSEA)))
#[1] 0

NAcc_EE_Betas_forGSEARanked<-NAcc_EE_Betas_forGSEA[order(NAcc_EE_Betas_forGSEA)]
head(NAcc_EE_Betas_forGSEARanked)
# Scn11a AABR07028009.1            Ak7     AC128059.4   LOC103693430        Tmem212 
# -2.236         -2.169         -2.026         -1.862         -1.845         -1.836 

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "TestRun_NAcc_EE_MainEffectsModel_FGSEAResults_BasicC5gmt.csv")
#It's working! Yessss...
#Now I just need to make our brain .gmt, because most of these pathways are completely irrelevant...

##################################

#Notes for when I come back to this:

#Reformat Dropviz (+Zeisel), StressGenes,  and Psychiatric genes into .gmt format and combine them.
#When doing BrainInABlender, remove NAs and duplicates
#Clean up c3gmt to be only gene sets related to the brain/neural function? (maybe clean up other gene sets too?)


DropViz_HC_Mouse_Temp<-list()
DropViz_HC_Rat_Temp<-list()
DropViz_HC_Human_Temp<-list()

names(table(DropViz_HC_CellTypeClusterMarkers_Concise_RatOrthologs$CellType))

for(i in c(1:length(names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))))){
  
  DropViz_HC_Mouse_Temp[[i]]<-DropViz_HC_CellTypeClusterMarkers_Concise$Symbol_Mouse[DropViz_HC_CellTypeClusterMarkers_Concise$CellType==names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))[i]]
  
  DropViz_HC_Rat_Temp[[i]]<-DropViz_HC_CellTypeClusterMarkers_Concise_RatOrthologs$Symbol_Rat[DropViz_HC_CellTypeClusterMarkers_Concise_RatOrthologs$CellType==names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))[i]]
  
  DropViz_HC_Human_Temp[[i]]<-DropViz_HC_CellTypeClusterMarkers_Concise_HumanOrthologs$Symbol_Human[DropViz_HC_CellTypeClusterMarkers_Concise_HumanOrthologs$CellType==names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))[i]]
  
}

names(DropViz_HC_Mouse_Temp)<-names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))
names(DropViz_HC_Rat_Temp)<-names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))
names(DropViz_HC_Human_Temp)<-names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))

str(DropViz_HC_Mouse_Temp)

library(stringi)

DropViz_HC_Mouse_GMT_Temp <- as.data.frame(stri_list2matrix(DropViz_HC_Mouse_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(DropViz_HC_Mouse_GMT_Temp)<-names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))
str(DropViz_HC_Mouse_GMT_Temp)

str(t(DropViz_HC_Mouse_GMT_Temp))
head(t(DropViz_HC_Mouse_GMT_Temp))
DropViz_HC_Mouse_GMT<-t(DropViz_HC_Mouse_GMT_Temp)
write.table(DropViz_HC_Mouse_GMT, "DropViz_HC_Mouse_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

DropViz_HC_Rat_GMT_Temp <- as.data.frame(stri_list2matrix(DropViz_HC_Rat_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(DropViz_HC_Rat_GMT_Temp)<-names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))
str(DropViz_HC_Rat_GMT_Temp)

str(t(DropViz_HC_Rat_GMT_Temp))
head(t(DropViz_HC_Rat_GMT_Temp))
DropViz_HC_Rat_GMT<-t(DropViz_HC_Rat_GMT_Temp)
write.table(DropViz_HC_Rat_GMT, "DropViz_HC_Rat_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

DropViz_HC_Human_GMT_Temp <- as.data.frame(stri_list2matrix(DropViz_HC_Human_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(DropViz_HC_Human_GMT_Temp)<-names(table(DropViz_HC_CellTypeClusterMarkers_Concise$CellType))
str(DropViz_HC_Human_GMT_Temp)

str(t(DropViz_HC_Human_GMT_Temp))
head(t(DropViz_HC_Human_GMT_Temp))
DropViz_HC_Human_GMT<-t(DropViz_HC_Human_GMT_Temp)
write.table(DropViz_HC_Human_GMT, "DropViz_HC_Human_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)


#Let's test it:
GMT_ForRats<-gmtPathways("DropViz_HC_Rat_GMT.gmt.txt")
str(GMT_ForRats)

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "TestRun_NAcc_EE_MainEffectsModel_FGSEAResults_DropViz_HC_CellTypes.csv")
#It's working! Yessssss


DropViz_STR_Mouse_Temp<-list()
DropViz_STR_Rat_Temp<-list()
DropViz_STR_Human_Temp<-list()

names(table(DropViz_STR_CellTypeClusterMarkers_Concise_RatOrthologs$CellType))

for(i in c(1:length(names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))))){
  
  DropViz_STR_Mouse_Temp[[i]]<-DropViz_STR_CellTypeClusterMarkers_Concise$Symbol_Mouse[DropViz_STR_CellTypeClusterMarkers_Concise$CellType==names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))[i]]
  
  DropViz_STR_Rat_Temp[[i]]<-DropViz_STR_CellTypeClusterMarkers_Concise_RatOrthologs$Symbol_Rat[DropViz_STR_CellTypeClusterMarkers_Concise_RatOrthologs$CellType==names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))[i]]
  
  DropViz_STR_Human_Temp[[i]]<-DropViz_STR_CellTypeClusterMarkers_Concise_HumanOrthologs$Symbol_Human[DropViz_STR_CellTypeClusterMarkers_Concise_HumanOrthologs$CellType==names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))[i]]
  
}

names(DropViz_STR_Mouse_Temp)<-names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))
names(DropViz_STR_Rat_Temp)<-names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))
names(DropViz_STR_Human_Temp)<-names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))

str(DropViz_STR_Mouse_Temp)

library(stringi)

DropViz_STR_Mouse_GMT_Temp <- as.data.frame(stri_list2matrix(DropViz_STR_Mouse_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(DropViz_STR_Mouse_GMT_Temp)<-names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))
str(DropViz_STR_Mouse_GMT_Temp)

str(t(DropViz_STR_Mouse_GMT_Temp))
head(t(DropViz_STR_Mouse_GMT_Temp))
DropViz_STR_Mouse_GMT<-t(DropViz_STR_Mouse_GMT_Temp)
write.table(DropViz_STR_Mouse_GMT, "DropViz_STR_Mouse_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

DropViz_STR_Rat_GMT_Temp <- as.data.frame(stri_list2matrix(DropViz_STR_Rat_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(DropViz_STR_Rat_GMT_Temp)<-names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))
str(DropViz_STR_Rat_GMT_Temp)

str(t(DropViz_STR_Rat_GMT_Temp))
head(t(DropViz_STR_Rat_GMT_Temp))
DropViz_STR_Rat_GMT<-t(DropViz_STR_Rat_GMT_Temp)
write.table(DropViz_STR_Rat_GMT, "DropViz_STR_Rat_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

DropViz_STR_Human_GMT_Temp <- as.data.frame(stri_list2matrix(DropViz_STR_Human_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(DropViz_STR_Human_GMT_Temp)<-names(table(DropViz_STR_CellTypeClusterMarkers_Concise$CellType))
str(DropViz_STR_Human_GMT_Temp)

str(t(DropViz_STR_Human_GMT_Temp))
head(t(DropViz_STR_Human_GMT_Temp))
DropViz_STR_Human_GMT<-t(DropViz_STR_Human_GMT_Temp)
write.table(DropViz_STR_Human_GMT, "DropViz_STR_Human_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

#Another Test:
GMT_ForRats<-gmtPathways("DropViz_STR_Rat_GMT.gmt.txt")
str(GMT_ForRats)

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "TestRun_NAcc_EE_MainEffectsModel_FGSEAResults_DropViz_STR_CellTypes.csv")
#Very pretty. :)

BrainInABlenderDatabase_All_JustSymbols_HumanOrthologs_NoNA
BrainInABlenderDatabase_All_JustSymbols_RatOrthologs_NoNA
BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA

BrainInABlender_Mouse_Temp<-list()
BrainInABlender_Rat_Temp<-list()
BrainInABlender_Human_Temp<-list()

names(table(BrainInABlenderDatabase_All_JustSymbols_RatOrthologs_NoNA$Tag))

for(i in c(1:length(names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))))){
  
  BrainInABlender_Mouse_Temp[[i]]<-BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Symbol_Mouse[BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag==names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))[i]]
  
  BrainInABlender_Rat_Temp[[i]]<-BrainInABlenderDatabase_All_JustSymbols_RatOrthologs_NoNA$Symbol_Rat[BrainInABlenderDatabase_All_JustSymbols_RatOrthologs_NoNA$Tag==names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))[i]]
  
  BrainInABlender_Human_Temp[[i]]<-BrainInABlenderDatabase_All_JustSymbols_HumanOrthologs_NoNA$Symbol_Human[BrainInABlenderDatabase_All_JustSymbols_HumanOrthologs_NoNA$Tag==names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))[i]]
  
}

names(BrainInABlender_Mouse_Temp)<-names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))
names(BrainInABlender_Rat_Temp)<-names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))
names(BrainInABlender_Human_Temp)<-names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))

str(BrainInABlender_Mouse_Temp)

library(stringi)

BrainInABlender_Mouse_GMT_Temp <- as.data.frame(stri_list2matrix(BrainInABlender_Mouse_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(BrainInABlender_Mouse_GMT_Temp)<-names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))
str(BrainInABlender_Mouse_GMT_Temp)

str(t(BrainInABlender_Mouse_GMT_Temp))
head(t(BrainInABlender_Mouse_GMT_Temp))
BrainInABlender_Mouse_GMT<-t(BrainInABlender_Mouse_GMT_Temp)
write.table(BrainInABlender_Mouse_GMT, "BrainInABlender_Mouse_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

BrainInABlender_Rat_GMT_Temp <- as.data.frame(stri_list2matrix(BrainInABlender_Rat_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(BrainInABlender_Rat_GMT_Temp)<-names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))
str(BrainInABlender_Rat_GMT_Temp)

str(t(BrainInABlender_Rat_GMT_Temp))
head(t(BrainInABlender_Rat_GMT_Temp))
BrainInABlender_Rat_GMT<-t(BrainInABlender_Rat_GMT_Temp)
write.table(BrainInABlender_Rat_GMT, "BrainInABlender_Rat_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

BrainInABlender_Human_GMT_Temp <- as.data.frame(stri_list2matrix(BrainInABlender_Human_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(BrainInABlender_Human_GMT_Temp)<-names(table(BrainInABlenderDatabase_All_JustSymbols_MouseOrthologs_NoNA$Tag))
str(BrainInABlender_Human_GMT_Temp)

str(t(BrainInABlender_Human_GMT_Temp))
head(t(BrainInABlender_Human_GMT_Temp))
BrainInABlender_Human_GMT<-t(BrainInABlender_Human_GMT_Temp)
write.table(BrainInABlender_Human_GMT, "BrainInABlender_Human_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)



InternalizingBehavior_Mouse_Temp<-list()
InternalizingBehavior_Rat_Temp<-list()
InternalizingBehavior_Human_Temp<-list()

names(table(InternalizingBehaviorGeneSets_RatOrthologs_NoNA$DatasetName))

for(i in c(1:length(names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))))){
  
  InternalizingBehavior_Mouse_Temp[[i]]<-InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$Symbol_Mouse[InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName==names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))[i]]
  
  InternalizingBehavior_Rat_Temp[[i]]<-InternalizingBehaviorGeneSets_RatOrthologs_NoNA$Symbol_Rat[InternalizingBehaviorGeneSets_RatOrthologs_NoNA$DatasetName==names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))[i]]
  
  InternalizingBehavior_Human_Temp[[i]]<-InternalizingBehaviorGeneSets_HumanOrthologs_NoNA$Symbol_Human[InternalizingBehaviorGeneSets_HumanOrthologs_NoNA$DatasetName==names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))[i]]
  
}

names(InternalizingBehavior_Mouse_Temp)<-names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))
names(InternalizingBehavior_Rat_Temp)<-names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))
names(InternalizingBehavior_Human_Temp)<-names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))

str(InternalizingBehavior_Mouse_Temp)

library(stringi)

InternalizingBehavior_Mouse_GMT_Temp <- as.data.frame(stri_list2matrix(InternalizingBehavior_Mouse_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(InternalizingBehavior_Mouse_GMT_Temp)<-names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))
str(InternalizingBehavior_Mouse_GMT_Temp)

str(t(InternalizingBehavior_Mouse_GMT_Temp))
head(t(InternalizingBehavior_Mouse_GMT_Temp))
InternalizingBehavior_Mouse_GMT<-t(InternalizingBehavior_Mouse_GMT_Temp)
write.table(InternalizingBehavior_Mouse_GMT, "InternalizingBehavior_Mouse_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

InternalizingBehavior_Rat_GMT_Temp <- as.data.frame(stri_list2matrix(InternalizingBehavior_Rat_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(InternalizingBehavior_Rat_GMT_Temp)<-names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))
str(InternalizingBehavior_Rat_GMT_Temp)

str(t(InternalizingBehavior_Rat_GMT_Temp))
head(t(InternalizingBehavior_Rat_GMT_Temp))
InternalizingBehavior_Rat_GMT<-t(InternalizingBehavior_Rat_GMT_Temp)
write.table(InternalizingBehavior_Rat_GMT, "InternalizingBehavior_Rat_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

InternalizingBehavior_Human_GMT_Temp <- as.data.frame(stri_list2matrix(InternalizingBehavior_Human_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(InternalizingBehavior_Human_GMT_Temp)<-names(table(InternalizingBehaviorGeneSets_MouseOrthologs_NoNA$DatasetName))
str(InternalizingBehavior_Human_GMT_Temp)

str(t(InternalizingBehavior_Human_GMT_Temp))
head(t(InternalizingBehavior_Human_GMT_Temp))
InternalizingBehavior_Human_GMT<-t(InternalizingBehavior_Human_GMT_Temp)
write.table(InternalizingBehavior_Human_GMT, "InternalizingBehavior_Human_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)


StressGeneSets_Mouse_Temp<-list()
StressGeneSets_Rat_Temp<-list()
StressGeneSets_Human_Temp<-list()

names(table(StressGeneSets_RatOrthologs_Symbols_NoNA$Dataset))

for(i in c(1:length(names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))))){
  
  StressGeneSets_Mouse_Temp[[i]]<-StressGeneSets_MouseOrthologs_Symbols_NoNA$Symbol_Mouse[StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset==names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))[i]]
  
  StressGeneSets_Rat_Temp[[i]]<-StressGeneSets_RatOrthologs_Symbols_NoNA$Symbol_Rat[StressGeneSets_RatOrthologs_Symbols_NoNA$Dataset==names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))[i]]
  
  StressGeneSets_Human_Temp[[i]]<-StressGeneSets_HumanOrthologs_Symbols_NoNA$Symbol_Human[StressGeneSets_HumanOrthologs_Symbols_NoNA$Dataset==names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))[i]]
  
}

names(StressGeneSets_Mouse_Temp)<-names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))
names(StressGeneSets_Rat_Temp)<-names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))
names(StressGeneSets_Human_Temp)<-names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))

str(StressGeneSets_Mouse_Temp)

library(stringi)

StressGeneSets_Mouse_GMT_Temp <- as.data.frame(stri_list2matrix(StressGeneSets_Mouse_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(StressGeneSets_Mouse_GMT_Temp)<-names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))
str(StressGeneSets_Mouse_GMT_Temp)

str(t(StressGeneSets_Mouse_GMT_Temp))
head(t(StressGeneSets_Mouse_GMT_Temp))
StressGeneSets_Mouse_GMT<-t(StressGeneSets_Mouse_GMT_Temp)
write.table(StressGeneSets_Mouse_GMT, "StressGeneSets_Mouse_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

StressGeneSets_Rat_GMT_Temp <- as.data.frame(stri_list2matrix(StressGeneSets_Rat_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(StressGeneSets_Rat_GMT_Temp)<-names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))
str(StressGeneSets_Rat_GMT_Temp)

str(t(StressGeneSets_Rat_GMT_Temp))
head(t(StressGeneSets_Rat_GMT_Temp))
StressGeneSets_Rat_GMT<-t(StressGeneSets_Rat_GMT_Temp)
write.table(StressGeneSets_Rat_GMT, "StressGeneSets_Rat_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

StressGeneSets_Human_GMT_Temp <- as.data.frame(stri_list2matrix(StressGeneSets_Human_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(StressGeneSets_Human_GMT_Temp)<-names(table(StressGeneSets_MouseOrthologs_Symbols_NoNA$Dataset))
str(StressGeneSets_Human_GMT_Temp)

str(t(StressGeneSets_Human_GMT_Temp))
head(t(StressGeneSets_Human_GMT_Temp))
StressGeneSets_Human_GMT<-t(StressGeneSets_Human_GMT_Temp)
write.table(StressGeneSets_Human_GMT, "StressGeneSets_Human_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)




HumanPsychiatricGeneSets_HumanOrthologs_JustSymbols_NoNA
HumanPsychiatricGeneSets_RatOrthologs_JustSymbols_NoNA
HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA

HumanPsychiatricGeneSets_Mouse_Temp<-list()
HumanPsychiatricGeneSets_Rat_Temp<-list()
HumanPsychiatricGeneSets_Human_Temp<-list()

names(table(HumanPsychiatricGeneSets_RatOrthologs_JustSymbols_NoNA$DatasetName))

for(i in c(1:length(names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))))){
  
  HumanPsychiatricGeneSets_Mouse_Temp[[i]]<-HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$Symbol_Mouse[HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName==names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))[i]]
  
  HumanPsychiatricGeneSets_Rat_Temp[[i]]<-HumanPsychiatricGeneSets_RatOrthologs_JustSymbols_NoNA$Symbol_Rat[HumanPsychiatricGeneSets_RatOrthologs_JustSymbols_NoNA$DatasetName==names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))[i]]
  
  HumanPsychiatricGeneSets_Human_Temp[[i]]<-HumanPsychiatricGeneSets_HumanOrthologs_JustSymbols_NoNA$Symbol_Human[HumanPsychiatricGeneSets_HumanOrthologs_JustSymbols_NoNA$DatasetName==names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))[i]]
  
}

names(HumanPsychiatricGeneSets_Mouse_Temp)<-names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))
names(HumanPsychiatricGeneSets_Rat_Temp)<-names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))
names(HumanPsychiatricGeneSets_Human_Temp)<-names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))

str(HumanPsychiatricGeneSets_Mouse_Temp)

library(stringi)

HumanPsychiatricGeneSets_Mouse_GMT_Temp <- as.data.frame(stri_list2matrix(HumanPsychiatricGeneSets_Mouse_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(HumanPsychiatricGeneSets_Mouse_GMT_Temp)<-names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))
str(HumanPsychiatricGeneSets_Mouse_GMT_Temp)

str(t(HumanPsychiatricGeneSets_Mouse_GMT_Temp))
head(t(HumanPsychiatricGeneSets_Mouse_GMT_Temp))
HumanPsychiatricGeneSets_Mouse_GMT<-t(HumanPsychiatricGeneSets_Mouse_GMT_Temp)
write.table(HumanPsychiatricGeneSets_Mouse_GMT, "HumanPsychiatricGeneSets_Mouse_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

HumanPsychiatricGeneSets_Rat_GMT_Temp <- as.data.frame(stri_list2matrix(HumanPsychiatricGeneSets_Rat_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(HumanPsychiatricGeneSets_Rat_GMT_Temp)<-names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))
str(HumanPsychiatricGeneSets_Rat_GMT_Temp)

str(t(HumanPsychiatricGeneSets_Rat_GMT_Temp))
head(t(HumanPsychiatricGeneSets_Rat_GMT_Temp))
HumanPsychiatricGeneSets_Rat_GMT<-t(HumanPsychiatricGeneSets_Rat_GMT_Temp)
write.table(HumanPsychiatricGeneSets_Rat_GMT, "HumanPsychiatricGeneSets_Rat_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = TRUE)

HumanPsychiatricGeneSets_Human_GMT_Temp <- as.data.frame(stri_list2matrix(HumanPsychiatricGeneSets_Human_Temp, byrow = FALSE), stringsAsFactors=FALSE)
names(HumanPsychiatricGeneSets_Human_GMT_Temp)<-names(table(HumanPsychiatricGeneSets_MouseOrthologs_JustSymbols_NoNA$DatasetName))
str(HumanPsychiatricGeneSets_Human_GMT_Temp)

str(t(HumanPsychiatricGeneSets_Human_GMT_Temp))
head(t(HumanPsychiatricGeneSets_Human_GMT_Temp))
HumanPsychiatricGeneSets_Human_GMT<-t(HumanPsychiatricGeneSets_Human_GMT_Temp)
write.table(HumanPsychiatricGeneSets_Human_GMT, "HumanPsychiatricGeneSets_Human_GMT.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)


#Note:

#To make a brain cell type GMT, for each species (mouse, rat, human) I will need to combine together DropViz_HC, DropViz_STR, BrainInABlender, Zeisel, and any relevant cell types in C8, and regional patterns from CombinedHCspecific_mouse
#Some of these things may be easier to do in Excel (e.g., adding Zeisel and the Cembrowski HC regional gene sets - or maybe just save those for function?).

DropViz_HC_Human_GMT, DropViz_STR_Human_GMT, BrainInABlender_Human_GMT 

DropViz_HC_Rat_GMT, DropViz_STR_Rat_GMT, BrainInABlender_Rat_GMT

DropViz_HC_Mouse_GMT,  DropViz_STR_Mouse_GMT, BrainInABlender_Mouse_GMT

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_toMakeRat")

c8_WhichAreRelevant<-read.delim("c8_WhichAreRelevant.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
str(c8_WhichAreRelevant)
#'data.frame':	673 obs. of  2 variables:
table(c8_WhichAreRelevant$V2)
# 0   1 
# 462 211

str(c8.all.v7.3.symbols)
# chr [1:673, 1:1772] "BUSSLINGER_ESOPHAGEAL_QUIESCENT_BASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_PROLIFERATING_BASAL_CELLS" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:673] "V1" "V2" "V3" "V4" ...
# ..$ : NULL
str(c8.all.v7.3.symbols_RatOrtholog)
#chr [1:673, 1:1872] "BUSSLINGER_ESOPHAGEAL_QUIESCENT_BASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_PROLIFERATING_BASAL_CELLS" ...
str(c8.all.v7.3.symbols_MouseOrtholog)
#chr [1:673, 1:1872] "BUSSLINGER_ESOPHAGEAL_QUIESCENT_BASAL_CELLS" "BUSSLINGER_ESOPHAGEAL_PROLIFERATING_BASAL_CELLS" ...

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")

c8.all.v7.3.symbols_BrainRelevant<-c8.all.v7.3.symbols[c8_WhichAreRelevant$V2==1,]
str(c8.all.v7.3.symbols_BrainRelevant)
# chr [1:211, 1:1772] "FAN_EMBRYONIC_CTX_BIG_GROUPS_CAJAL_RETZIUS" "FAN_EMBRYONIC_CTX_BIG_GROUPS_BRAIN_ENDOTHELIAL" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:211] "V39" "V40" "V41" "V42" ...
# ..$ : NULL

c8.all.v7.3.symbols_MouseOrtholog_BrainRelevant<-c8.all.v7.3.symbols_MouseOrtholog[c8_WhichAreRelevant$V2==1,]
str(c8.all.v7.3.symbols_MouseOrtholog_BrainRelevant)

c8.all.v7.3.symbols_RatOrtholog_BrainRelevant<-c8.all.v7.3.symbols_RatOrtholog[c8_WhichAreRelevant$V2==1,]
str(c8.all.v7.3.symbols_RatOrtholog_BrainRelevant)

write.table(c8.all.v7.3.symbols_BrainRelevant, "c8.all.v7.3.symbols_BrainRelevant.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

write.table(c8.all.v7.3.symbols_MouseOrtholog_BrainRelevant, "c8.all.v7.3.symbols_MouseOrtholog_BrainRelevant.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

write.table(c8.all.v7.3.symbols_RatOrtholog_BrainRelevant, "c8.all.v7.3.symbols_RatOrtholog_BrainRelevant.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

#To make a brain cell type GMT, for each species (mouse, rat, human) I will need to combine together DropViz_HC, DropViz_STR, BrainInABlender, Zeisel, and any relevant cell types in C8, and regional patterns from CombinedHCspecific_mouse
#Some of these things may be easier to do in Excel (e.g., adding Zeisel and the Cembrowski HC regional gene sets - or maybe just save those for function?).

#Objects of interest for Brain Cell Type GMT:
DropViz_HC_Human_GMT, DropViz_STR_Human_GMT, BrainInABlender_Human_GMT, c8.all.v7.3.symbols_BrainRelevant, Zeisel_CA1PyramidalIndex_HumanOrthologs, CombinedHCspecific_Mouse_HumanOrtholog 

DropViz_HC_Rat_GMT, DropViz_STR_Rat_GMT, BrainInABlender_Rat_GMT, c8.all.v7.3.symbols_RatOrtholog_BrainRelevant, Zeisel_CA1PyramidalIndex_RatOrthologs, CombinedHCspecific_Mouse_RatOrtholog

DropViz_HC_Mouse_GMT,  DropViz_STR_Mouse_GMT, BrainInABlender_Mouse_GMT, c8.all.v7.3.symbols_MouseOrtholog_BrainRelevant, Zeisel_CA1PyramidalIndex, CombinedHCspecific_Mouse

#Notes about structure:
#The Dropviz and BrainInABlender based GMT files have the actual name as the row.names, no link. I wonder if my test-run with these in fgsea skipped the first two columns of data...
#C8 and CombinedHCSpecific: first 2 columns are the name and the link.
#CombinedHCSpecific: First 14 rows are Hipposeq
#Zeisel_CA1PyramidalIndex_HumanOrthologs needs editing - extra columns of crap.

#It seems like this might actually be easier to do in Excel. I just sorted all of the GMTs by species and am going to try messing around there.

#Let's see if it worked:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("BrainCellTypes_V1_Rat_GMT.gmt.txt")
str(GMT_ForRats)
# List of 289
# $ Astrocyte_All_Cahoy_JNeuro_2008                                           : chr [1:1213] "Gfap" "Aqp4" "Pla2g7" "Slc39a12" ...
# $ Astrocyte_All_Darmanis_PNAS_2015                                          : chr [1:1213] "Fgfr3" "Aqp4" "Gja1" "Agt" ...
# $ Astrocyte_All_Doyle_Cell_2008                                             : chr [1:1213] "Adora2b" "Cyp4f4" "Fzd2" "Nr2e1" ...
#Looks good.

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")
write.csv(temp1, "TestRun_NAcc_EE_MainEffectsModel_FGSEAResults_BrainCellTypesGMT.csv")
#Works beautifully. :)


setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_toMakeRat")

c2_WhichAreRelevant<-read.csv("c2_WhichAreRelevant_NonCuratedBrain.csv", header=FALSE, stringsAsFactors = FALSE)
str(c2_WhichAreRelevant)
# 'data.frame':	6255 obs. of  2 variables:
#   $ V1: chr  "BENITEZ_GBM_PROTEASOME_INHIBITION_RESPONSE" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_UP" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_DN" "BLANCO_MELO_MERS_COV_INFECTION_MCR5_CELLS_UP" ...
# $ V2: int  0 0 0 0 0 0 0 0 0 0 ...
table(c2_WhichAreRelevant$V2)
# 0    1 
# 6097  158 

str(c2.all.v7.3.symbols)
# chr [1:6255, 1:1943] "BENITEZ_GBM_PROTEASOME_INHIBITION_RESPONSE" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_UP" ...
# - attr(*, "dimnames")=List of 2
str(c2.all.v7.3.symbols_MouseOrtholog)
# chr [1:6255, 1:2043] "BENITEZ_GBM_PROTEASOME_INHIBITION_RESPONSE" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_UP" ...
str(c2.all.v7.3.symbols_RatOrtholog)
#chr [1:6255, 1:2043] "BENITEZ_GBM_PROTEASOME_INHIBITION_RESPONSE" "BLANCO_MELO_SARS_COV_1_INFECTION_MCR5_CELLS_UP" ...


c2.all.v7.3.symbols_BrainRelevant<-c2.all.v7.3.symbols[c2_WhichAreRelevant$V2==1,]
str(c2.all.v7.3.symbols_BrainRelevant)
# chr [1:158, 1:1943] "RHEIN_ALL_GLUCOCORTICOID_THERAPY_UP" "RHEIN_ALL_GLUCOCORTICOID_THERAPY_DN" "CHEBOTAEV_GR_TARGETS_UP" "CHEBOTAEV_GR_TARGETS_DN" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:158] "V205" "V206" "V329" "V330" ...
# ..$ : NULL

c2.all.v7.3.symbols_MouseOrtholog_BrainRelevant<-c2.all.v7.3.symbols_MouseOrtholog[c2_WhichAreRelevant$V2==1,]
str(c2.all.v7.3.symbols_MouseOrtholog_BrainRelevant)

c2.all.v7.3.symbols_RatOrtholog_BrainRelevant<-c2.all.v7.3.symbols_RatOrtholog[c2_WhichAreRelevant$V2==1,]
str(c2.all.v7.3.symbols_RatOrtholog_BrainRelevant)

write.table(c2.all.v7.3.symbols_BrainRelevant, "c2.all.v7.3.symbols_BrainRelevant.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

write.table(c2.all.v7.3.symbols_MouseOrtholog_BrainRelevant, "c2.all.v7.3.symbols_MouseOrtholog_BrainRelevant.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

write.table(c2.all.v7.3.symbols_RatOrtholog_BrainRelevant, "c2.all.v7.3.symbols_RatOrtholog_BrainRelevant.gmt.txt", na="", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

#To make a brain function GMT, for each species (mouse, rat, human) I will need to combine together StressGeneSets, InternalizingBehaviorGeneSets, PsychiatricGeneSets, co-expression patterns from CombinedHCspecific_mouse and CombinedHCspecific_human, and any brain relevant datasets from C2.

#I tried combining them in Excel first


#Add as a product to go along with the "With our powers combined commentary"?

#Let's see if it worked:

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("BrainFunction_CustomGeneSets_V1_Rat_GMT.gmt.txt")
str(GMT_ForRats)
# List of 260
# $ Gandal_2018_AlcoholAbuseDisorder_Downregulated_Cortex         : chr [1:2082] "Cyp51" "Slc25a5" "Pik3c2a" "Matr3" ...
# $ Gandal_2018_AlcoholAbuseDisorder_Upregulated_Cortex           : chr [1:2082] "Upf1" "Gas7" "Camk1g" "Med24" ...
# $ Gandal_2018_AutismSpectrumDisorder_Downregulated_Cortex       : chr [1:2082] "Mad1l1" "Ica1" "Ndufab1" "Pdk2" ...

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")
write.csv(temp1, "TestRun_NAcc_EE_MainEffectsModel_FGSEAResults_BrainFunctionGMT.csv")
#Works beautifully. :)

#Since these two gene sets are dinky, I also tried combining them:
setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325/GMTs_Rat")
GMT_ForRats<-gmtPathways("BrainCellTypesAndFunction_CustomGeneSets_V1_Rat_GMT.gmt.txt")
str(GMT_ForRats)
# List of 549
# $ Astrocyte_All_Cahoy_JNeuro_2008                                           : chr [1:2082] "Gfap" "Aqp4" "Pla2g7" "Slc39a12" ...
# $ Astrocyte_All_Darmanis_PNAS_2015                                          : chr [1:2082] "Fgfr3" "Aqp4" "Gja1" "Agt" ...
# $ Astrocyte_All_Doyle_Cell_2008                                             : chr [1:2082] "Adora2b" "Cyp4f4" "Fzd2" "Nr2e1" ...

temp1<-fgsea(GMT_ForRats, NAcc_EE_Betas_forGSEARanked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/NAcc_Remove2outliers/11_FGSEA/UpdatedGMT_20210325")
write.csv(temp1, "TestRun_NAcc_EE_MainEffectsModel_FGSEAResults_BrainCellsAndFunctionGMT.csv")



######################
