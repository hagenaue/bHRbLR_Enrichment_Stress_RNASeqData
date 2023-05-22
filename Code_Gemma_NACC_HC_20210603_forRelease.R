#Code for extracting gene sets from Gemma and GeneWeaver from the Hippocampus and Nucleus Accumbens to create a better .GMT file
#Megan Hagenauer, June 2021

library(plyr)
library(httr)
library(RCurl)
library(jsonlite)

library(devtools)
install_github("PavlidisLab/gemmaAPI.R")
library(gemmaAPI)

?endpointFunctions
?highLevelFunctions

setGemmaUser(username="hagenaue", password=[DeletedForPublicRelease]) 

httr::set_config(config(ssl_verifypeer = 0L))
# This code prevents this curl error from occurring with some of Gemma's API functions: 
# Error in curl::curl_fetch_memory(url, handle = handle) : 
# SSL certificate problem: unable to get local issuer certificate

NACC_Annotation<-annotationInfo(annotation='nucleus accumbens', request='datasets')

length(names(NACC_Annotation))
#[1] 103

NACC_Annotation_AsDF<-data.frame(t(sapply(NACC_Annotation,c)))

dim(NACC_Annotation_AsDF)
#[1] 103  36

head(NACC_Annotation_AsDF)

NACC_Annotation_AsDF$accession[sapply(NACC_Annotation_AsDF$accession, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$accession))

length(unlist(NACC_Annotation_AsDF$taxon))
NACC_Annotation_AsDF$taxon[sapply(NACC_Annotation_AsDF$taxon, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$taxon))

length(unlist(NACC_Annotation_AsDF$batchConfound))
NACC_Annotation_AsDF$batchConfound[sapply(NACC_Annotation_AsDF$batchConfound, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$batchConfound))

length(unlist(NACC_Annotation_AsDF$$batchEffect))
NACC_Annotation_AsDF$batchEffect[sapply(NACC_Annotation_AsDF$batchEffect, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$batchEffect))

length(unlist(NACC_Annotation_AsDF$bioMaterialCount))
NACC_Annotation_AsDF$bioMaterialCount[sapply(NACC_Annotation_AsDF$bioMaterialCount, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$bioMaterialCount))

length(unlist(NACC_Annotation_AsDF$processedExpressionVectorCount))
NACC_Annotation_AsDF$processedExpressionVectorCount[sapply(NACC_Annotation_AsDF$processedExpressionVectorCount, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$processedExpressionVectorCount))

length(unlist(NACC_Annotation_AsDF$isPublic))
NACC_Annotation_AsDF$isPublic[sapply(NACC_Annotation_AsDF$isPublic, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$isPublic))

length(unlist(NACC_Annotation_AsDF$name))
NACC_Annotation_AsDF$name[sapply(NACC_Annotation_AsDF$name, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$name))

length(unlist(NACC_Annotation_AsDF$technologyType))
NACC_Annotation_AsDF$technologyType[sapply(NACC_Annotation_AsDF$technologyType, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$technologyType))

length(unlist(NACC_Annotation_AsDF$troubled))
NACC_Annotation_AsDF$troubled[sapply(NACC_Annotation_AsDF$troubled, is.null)] <- NA
length(unlist(NACC_Annotation_AsDF$troubled))

NACC_Annotation_AsDF_Simpler<-cbind.data.frame(unlist(NACC_Annotation_AsDF$shortName), unlist(NACC_Annotation_AsDF$accession), unlist(NACC_Annotation_AsDF$name), unlist(NACC_Annotation_AsDF$isPublic), unlist(NACC_Annotation_AsDF$technologyType),  unlist(NACC_Annotation_AsDF$taxon), unlist(NACC_Annotation_AsDF$batchConfound), unlist(NACC_Annotation_AsDF$batchEffect), unlist(NACC_Annotation_AsDF$bioMaterialCount), unlist(NACC_Annotation_AsDF$processedExpressionVectorCount), unlist(NACC_Annotation_AsDF$troubled), stringsAsFactors=FALSE)

str(NACC_Annotation_AsDF_Simpler)

#let's rename those columns (although I'm sure there is a better way to do this...)
colnames(NACC_Annotation_AsDF_Simpler)<-c("shortName", "accession", "name", "isPublic", "technologyType", "taxon", "batchConfound", "batchEffect", "bioMaterialCount", "processedExpressionVectorCount","troubled")

#Maybe we should trim it down to what we are most likely to use before moving further:
#double checking basic exclusion criteria:

table(NACC_Annotation_AsDF_Simpler$isPublic)
# TRUE 
# 103 

table(NACC_Annotation_AsDF_Simpler$troubled)
# FALSE 
# 103 

table(NACC_Annotation_AsDF_Simpler$taxon)
# human mouse rat 
# 8    67    28

#How about more complicated criteria?
table(NACC_Annotation_AsDF_Simpler$batchConfound)
#It looks like this may need to be navigated on a case-by-case basis - most of the confounds are factor-specific, and especially show up in relationship to organism part, tissue, sex, and time point.  (although some are clearly more relevant to our questions)

table(NACC_Annotation_AsDF_Simpler$technologyType)
# GENELIST ONECOLOR TWOCOLOR 
# 43       59        1
#Gene list means RNA-Seq, so about half of these are RNA-Seq studies.

setwd("~/Documents/SideProjects/BrainGMT/Gemma")
write.csv(NACC_Annotation_AsDF_Simpler, "NACC_Annotation_AsDF_Simpler.csv")


#####################################

#I had trouble using the api to download differential expression results on mass. Here is the code that should work:
DEResults = 
  datasetInfo('GSE107999',
              request='degs', # we want this endpoint to return data. see documentation
              differential='94808', #from DEComparisons[[1]]$id
              return = TRUE, # TRUE by default, all functions have this. if false there'll be no return
              file = NULL # NULL by default, all functions have this. If specificed, output will be saved.
  )
#Error in loadNamespace(name) : there is no package called ‘chromote’


install_github("rstudio/chromote")
# ERROR: dependencies ‘websocket’, ‘fastmap’ are not available for package ‘chromote’
# * removing ‘/Users/mhh/Library/R/3.4/library/chromote’
# Installation failed: Command failed (1)

#I think I just need to update. Sigh. 

#I tried installing websocket and it failed, I think because of the SSL option that I changed to make the gemma code work. I'll need to come back to this.
#Configuration failed because openssl was not found. 

#############

#A temporary work around:
#I downloaded the results for the datasets of interest (Nucleus accumbens + variables/conditions resembling our own) and unzipped them, then deleted the zip files 


library(plyr)

setwd("~/Documents/SideProjects/BrainGMT/Gemma/NucleusAccumbens")

NACC_DatasetsOfInterest<-list.files()

for(i in c(16:length(NACC_DatasetsOfInterest))){
  
  setwd("~/Documents/SideProjects/BrainGMT/Gemma/NucleusAccumbens")
  
  print(i)
  print(NACC_DatasetsOfInterest[i])
  
  setwd(paste("./", NACC_DatasetsOfInterest[i], sep=""))
  
  TempResultsFiles<-list.files()
  
  TempAnalysisResults<-read.delim("analysis.results.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")
  
  TempResultsToJoin<-list(TempAnalysisResults)
  str(TempResultsToJoin)
  
  for(j in c(2:length(TempResultsFiles))){
    print(k)
    TempResultsToJoin[[j]]<-assign(paste("TempResultSet", j, sep=""), read.delim(TempResultsFiles[j], sep="\t", stringsAsFactors = FALSE, comment.char = "#"))
  }
  
  TempResultsJoined<-join_all(TempResultsToJoin, by="Element_Name")
  
  write.csv(TempResultsJoined, "TempResultsJoined.csv")
  
  rm(TempResultsFiles, TempAnalysisResults, TempResultsToJoin)
  
  TempIndicesForQValueColumns<-grep("QValue", colnames(TempResultsJoined),)

  for(k in c(TempIndicesForQValueColumns)){
    print(k)
    print(colnames(TempResultsJoined)[k])
    if(sum(TempResultsJoined[,k]<0.10 & TempResultsJoined[,k+1]<0.0001, na.rm=TRUE)>0){
      write.csv(TempResultsJoined[TempResultsJoined[,k]<0.10 & TempResultsJoined[,k+1]<0.0001,], paste(colnames(TempResultsJoined)[k], "_10.csv", sep=""))
    }else{}
    
  }
  
  rm(TempResultsJoined, TempIndicesForQValueColumns)
  
}

# [1] 15
# [1] "693_GSE8870_diffExpAnalysis_29489"
# Error in read.table(file = file, header = header, sep = sep, quote = quote,  : 
#                       more columns than column names

#for some reason this file isn't reading in properly, so I skipped it.
#it was the first dataset that I downloaded, and I think I might have messed around with the files so some of the comments weren't hashtagged. I added some hashtags and it seems to work now.

i<-15

#Other steps I should add to the automated code to be one step closer to actual gene sets.
#1. Remove results without gene symbols
#2. Remove results with multiple mappings
#3. Deal with "date genes"
#4. For results with multiple factor levels, pull out results that are nominally significant for each level.
#4. Count duplicates, and for any genes with duplicates showing the same direction of effect, choose the most extreme value.
#5. Decide whether there are sufficient results to split results into down and upregulated
#6. Split results when relevant.
#7. Output the final gene set sizes.

#Questions for Paul:
#1: Are they taking into account the multi-level nature of the data from multiple organism parts? (says ANOVA)
#2: Why do the contrast p-values for factors with only 2 levels not match the ANOVA output?
#3: How do we interpret the contrast output for models with interaction terms?  (main effect vs. effect within reference group?)



###########

setGemmaUser(username="hagenaue", password=[DeletedForPublicRelease]) 

httr::set_config(config(ssl_verifypeer = 0L))
# This code prevents this curl error from occurring with some of Gemma's API functions: 
# Error in curl::curl_fetch_memory(url, handle = handle) : 
# SSL certificate problem: unable to get local issuer certificate

HC_Annotation<-annotationInfo(annotation='hippocampus', request='datasets')

length(names(HC_Annotation))
#[1] 648

HC_Annotation_AsDF<-data.frame(t(sapply(HC_Annotation,c)))

dim(HC_Annotation_AsDF)
#[1] 648  36

head(HC_Annotation_AsDF)

HC_Annotation_AsDF$accession[sapply(HC_Annotation_AsDF$accession, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$accession))

length(unlist(HC_Annotation_AsDF$taxon))
HC_Annotation_AsDF$taxon[sapply(HC_Annotation_AsDF$taxon, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$taxon))

length(unlist(HC_Annotation_AsDF$batchConfound))
HC_Annotation_AsDF$batchConfound[sapply(HC_Annotation_AsDF$batchConfound, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$batchConfound))

length(unlist(HC_Annotation_AsDF$$batchEffect))
HC_Annotation_AsDF$batchEffect[sapply(HC_Annotation_AsDF$batchEffect, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$batchEffect))

length(unlist(HC_Annotation_AsDF$bioMaterialCount))
HC_Annotation_AsDF$bioMaterialCount[sapply(HC_Annotation_AsDF$bioMaterialCount, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$bioMaterialCount))

length(unlist(HC_Annotation_AsDF$processedExpressionVectorCount))
HC_Annotation_AsDF$processedExpressionVectorCount[sapply(HC_Annotation_AsDF$processedExpressionVectorCount, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$processedExpressionVectorCount))

length(unlist(HC_Annotation_AsDF$isPublic))
HC_Annotation_AsDF$isPublic[sapply(HC_Annotation_AsDF$isPublic, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$isPublic))

length(unlist(HC_Annotation_AsDF$name))
HC_Annotation_AsDF$name[sapply(HC_Annotation_AsDF$name, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$name))

length(unlist(HC_Annotation_AsDF$technologyType))
HC_Annotation_AsDF$technologyType[sapply(HC_Annotation_AsDF$technologyType, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$technologyType))

length(unlist(HC_Annotation_AsDF$troubled))
HC_Annotation_AsDF$troubled[sapply(HC_Annotation_AsDF$troubled, is.null)] <- NA
length(unlist(HC_Annotation_AsDF$troubled))

HC_Annotation_AsDF_Simpler<-cbind.data.frame(unlist(HC_Annotation_AsDF$shortName), unlist(HC_Annotation_AsDF$accession), unlist(HC_Annotation_AsDF$name), unlist(HC_Annotation_AsDF$isPublic), unlist(HC_Annotation_AsDF$technologyType),  unlist(HC_Annotation_AsDF$taxon), unlist(HC_Annotation_AsDF$batchConfound), unlist(HC_Annotation_AsDF$batchEffect), unlist(HC_Annotation_AsDF$bioMaterialCount), unlist(HC_Annotation_AsDF$processedExpressionVectorCount), unlist(HC_Annotation_AsDF$troubled), stringsAsFactors=FALSE)

str(HC_Annotation_AsDF_Simpler)
#'data.frame':	648 obs. of  11 variables:
#'
#let's rename those columns (although I'm sure there is a better way to do this...)
colnames(HC_Annotation_AsDF_Simpler)<-c("shortName", "accession", "name", "isPublic", "technologyType", "taxon", "batchConfound", "batchEffect", "bioMaterialCount", "processedExpressionVectorCount","troubled")

#Maybe we should trim it down to what we are most likely to use before moving further:
#double checking basic exclusion criteria:

table(HC_Annotation_AsDF_Simpler$isPublic)
# TRUE 
#648 

table(HC_Annotation_AsDF_Simpler$troubled)
# FALSE 
#648  

table(HC_Annotation_AsDF_Simpler$taxon)
# human mouse   rat 
# 44   460   143 

#How about more complicated criteria?
table(HC_Annotation_AsDF_Simpler$batchConfound)
#It looks like this may need to be navigated on a case-by-case basis - most of the confounds are factor-specific, and especially show up in relationship to organism part, tissue, sex, and time point.  (although some are clearly more relevant to our questions)

table(HC_Annotation_AsDF_Simpler$technologyType)
# DUALMODE GENELIST ONECOLOR TWOCOLOR 
# 13      254      371       10 
#Gene list means RNA-Seq, so about half of these are RNA-Seq studies.

setwd("~/Documents/SideProjects/BrainGMT/Gemma")
write.csv(HC_Annotation_AsDF_Simpler, "HC_Annotation_AsDF_Simpler.csv")
#Outputted 06-15-2021



#Adapting the code from the Nucleus Accumbens analysis to make even more steps automated...


setwd("~/Documents/SideProjects/BrainGMT/Gemma")

DateGeneDecoder<-read.csv("DateGeneDecoder.txt", header=TRUE, stringsAsFactors = FALSE)

DateGeneDecoder_Human<-DateGeneDecoder
DateGeneDecoder_Human$OldSymbol<-toupper(DateGeneDecoder$OldSymbol)
DateGeneDecoder_Human$NewSymbol<-toupper(DateGeneDecoder$NewSymbol)


write(paste(DateGeneDecoder[,1],"='", DateGeneDecoder[,2], "',", sep=""), "DateGeneDecoder_Mouse_String.txt")
write(paste(DateGeneDecoder_Human[,1],"='", DateGeneDecoder_Human[,2], "',", sep=""), "DateGeneDecoder_Human_String.txt")
#This is a situation where I could have used the "collapse" function from the fstrings package to remove the unwanted "", but I couldn't install it because my entire system is out of date. I really need to fix that. But for, now, I'm just going to quickly take care of things by outputting to .txt and then come back and take care of the bigger problem later.

#ooh - here's another way to do it:
toString(paste(DateGeneDecoder[,1],"='", DateGeneDecoder[,2], "'", sep=""))
toString(paste(DateGeneDecoder_Human[,1],"='", DateGeneDecoder_Human[,2], "'", sep=""))


library(plyr)

setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus")

HC_DatasetsOfInterest<-list.files()
#86

HC_GeneSets_Temp<-data.frame(Dataset_Name=character(), Factor_Contrast=character(), Direction=character(), Taxon=character(), GeneSetSize=numeric(), GeneSet=character())

#Loops over all folders (publications) of down-loaded results:

#These folders are problematic because they completely lack gene symbols:
c(28,65,69)

for(i in c(1:length(HC_DatasetsOfInterest))){
  
  setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus")
  
  print(i)
  print(HC_DatasetsOfInterest[i])
  
  #Pulling out relevant meta-data for the dataset analysis:
  Dataset_ID<-strsplit(HC_DatasetsOfInterest[i], "_")[[1]]
  names(Dataset_ID)<-c("Gemma_ID", "GEO_Accession", "Output", "Result_ID")
  
  #This code bugs out when Gemma has added a .1 etc to denote more than one entry for a dataset, so I had to add a string split.
  
  Dataset_Annotation<-HC_Annotation_AsDF_Simpler[which(HC_Annotation_AsDF_Simpler$accession==strsplit(Dataset_ID[2][[1]], ".", fixed=TRUE)[[1]][1]),]
  
  Dataset_Name<-paste(Dataset_Annotation$accession, "_", gsub( " ", "_", Dataset_Annotation$name, fixed=TRUE), sep="")
  
  setwd(paste("./", HC_DatasetsOfInterest[i], sep=""))
  
  TempResultsFiles<-list.files()
  
  #Pulls out the ANOVA results:
  TempAnalysisResults<-read.delim("analysis.results.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")
  
  TempResultsToJoin<-list(TempAnalysisResults)
  str(TempResultsToJoin)
  
  #Joins the ANOVA results to the results from the individual contrasts:
  for(j in c(2:length(TempResultsFiles))){
    print(k)
    TempResultsToJoin[[j]]<-assign(paste("TempResultSet", j, sep=""), read.delim(TempResultsFiles[j], sep="\t", stringsAsFactors = FALSE, comment.char = "#"))
  }
  
  TempResultsJoined<-join_all(TempResultsToJoin, by="Element_Name")
  
  write.csv(TempResultsJoined, "TempResultsJoined.csv")
  
  rm(TempResultsFiles, TempAnalysisResults, TempResultsToJoin)
  
  TempIndicesForQValueColumns<-grep("QValue", colnames(TempResultsJoined))
  TempIndicesForFoldChangeColumns<-grep("FoldChange", colnames(TempResultsJoined))
  
  #Removing rows lacking gene symbol annotation, with multi-mapped annotation, or "date" genes that need fixing:
  
  if(sum(is.na(TempResultsJoined$Gene_Symbol)==FALSE & (TempResultsJoined$Gene_Symbol=="")==FALSE & (TempResultsJoined$Gene_Symbol==" ")==FALSE)>0){
  
  TempResultsJoined_NoNA<-TempResultsJoined[is.na(TempResultsJoined$Gene_Symbol)==FALSE & (TempResultsJoined$Gene_Symbol=="")==FALSE & (TempResultsJoined$Gene_Symbol==" ")==FALSE & (TempResultsJoined$Gene_Symbol=="null")==FALSE & (TempResultsJoined$Gene_Symbol=="N/A")==FALSE,]
  
  rm(TempResultsJoined)
  
  TempResultsJoined_NoNA_NoMultimapped<-TempResultsJoined_NoNA[grep("[|]", TempResultsJoined_NoNA$Gene_Symbol, invert=TRUE),]
  
  rm(TempResultsJoined_NoNA)
  
  #To replace date genes, I need to first determine whether the gene symbols reflect humans or mice/rats, because the Excel perversion of the date genes is all lower case (i.e., the same in both species)
  #On the other hand, since this is all Gemma output, I'm pretty sure there won't be excel versions of the gene names. So maybe just skip that step and flag if there is an issue?
  
  sum(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol%in%DateGeneDecoder[,1])
  
  TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol<-recode(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, Sept1='Septin1',
                               Sept2='Septin2',
                               Sept3='Septin3',
                               Sept4='Septin4',
                               Sept5='Septin5',
                               Sept6='Septin6',
                               Sept7='Septin7',
                               Sept8='Septin8',
                               Sept9='Septin9',
                               Sept10='Septin10',
                               Sept11='Septin11',
                               Sept12='Septin12',
                               Sept13='Septin13',
                               Sept14='Septin14',
                               Sept7p2='Septin7p2',
                               March1='Marchf1',
                               March2='Marchf2',
                               March3='Marchf3',
                               March4='Marchf4',
                               March5='Marchf5',
                               March6='Marchf6',
                               March7='Marchf7',
                               March8='Marchf8',
                               March9='Marchf9',
                               March10='Marchf10',
                               March11='Marchf11',
                               Dec1='Delec1',
                               Nov='Ccn3')
  
  TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol<-recode(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, SEPT1='SEPTIN1',
                               SEPT2='SEPTIN2',
                               SEPT3='SEPTIN3',
                               SEPT4='SEPTIN4',
                               SEPT5='SEPTIN5',
                               SEPT6='SEPTIN6',
                               SEPT7='SEPTIN7',
                               SEPT8='SEPTIN8',
                               SEPT9='SEPTIN9',
                               SEPT10='SEPTIN10',
                               SEPT11='SEPTIN11',
                               SEPT12='SEPTIN12',
                               SEPT13='SEPTIN13',
                               SEPT14='SEPTIN14',
                               SEPT7P2='SEPTIN7P2',
                               MARCH1='MARCHF1',
                               MARCH2='MARCHF2',
                               MARCH3='MARCHF3',
                               MARCH4='MARCHF4',
                               MARCH5='MARCHF5',
                               MARCH6='MARCHF6',
                               MARCH7='MARCHF7',
                               MARCH8='MARCHF8',
                               MARCH9='MARCHF9',
                               MARCH10='MARCHF10',
                               MARCH11='MARCHF11',
                               DEC1='DELEC1',
                               NOV='CCN3')

  sum(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol%in%DateGeneDecoder[,1])
 
  write.csv(TempResultsJoined_NoNA_NoMultimapped, "TempResultsJoined_GeneSymbolsClean.csv")
  
  #Overwrite the other version of the df to reuse previously-written code:
  TempResultsJoined<-TempResultsJoined_NoNA_NoMultimapped

  rm(TempResultsJoined_NoNA_NoMultimapped)
  
  #Loops over each variable and pulls out the rows of genes with q<0.10 (FDR<0.10) and p<0.0001 (to help cut down on the # of genes in results files with thousands of "sig" results):
  
  for(k in c(TempIndicesForQValueColumns)){
    print(k)
    print(colnames(TempResultsJoined)[k])
    if(sum(TempResultsJoined[,k]<0.10 & TempResultsJoined[,k+1]<0.0001, na.rm=TRUE)>0){
      TempResultsJoined_ResultK<-TempResultsJoined[TempResultsJoined[,k]<0.10 & TempResultsJoined[,k+1]<0.0001,]
      write.csv(TempResultsJoined_ResultK, paste(colnames(TempResultsJoined)[k], "_10.csv", sep=""))
      
      #For the genes that show a significant relationship with a variable (q<0.10 & p<0.0001), this code loops over each level of contrast and pulls out the results that are p<0.05 for that contrast. 
      #Note that only a handful of these combinations are actually meaningful - I couldn't figure out a good way to automatically match variable + relevant contrasts.
      
      for(m in c(TempIndicesForFoldChangeColumns)){
        print(m)
      
        TempContrast<-paste(gsub("QValue_", "", colnames(TempResultsJoined)[k]), gsub("FoldChange_", "", colnames(TempResultsJoined)[m]), sep="_")
        print(TempContrast)
        
        #Pulling out down-regulated genes for each contrast:
        if(sum(TempResultsJoined_ResultK[,m]<0 & TempResultsJoined_ResultK[,m+2]<0.05, na.rm=TRUE)>0){
        TempResultsJoined_ResultK_FoldChangeM_Down<-TempResultsJoined_ResultK[TempResultsJoined_ResultK[,m]<0 & TempResultsJoined_ResultK[,m+2]<0.05,]
        write.csv(TempResultsJoined_ResultK_FoldChangeM_Down, paste(colnames(TempResultsJoined)[k], "_10_w_", colnames(TempResultsJoined)[m], "_Down.csv", sep=""))
        
        if(nrow(TempResultsJoined_ResultK_FoldChangeM_Down)>999){
          TempResultsJoined_ResultK_FoldChangeM_Down_Reordered<-TempResultsJoined_ResultK_FoldChangeM_Down[order(TempResultsJoined_ResultK_FoldChangeM_Down[,m]),]
          TempResultsJoined_ResultK_FoldChangeM_Down<-TempResultsJoined_ResultK_FoldChangeM_Down_Reordered[c(1:999),]
          rm(TempResultsJoined_ResultK_FoldChangeM_Down_Reordered)
          }else{}
        
        TempResultsJoined_ResultK_FoldChangeM_Down_Symbols<-unique(TempResultsJoined_ResultK_FoldChangeM_Down$Gene_Symbol)
        
        write.csv(TempResultsJoined_ResultK_FoldChangeM_Down_Symbols, paste(colnames(TempResultsJoined)[k], "_10_w_", colnames(TempResultsJoined)[m], "_Down_UniqueSymbols.csv", sep=""))
        
        TempResultsJoined_ResultK_FoldChangeM_DownDF<-data.frame(Dataset_Name=Dataset_Name, Factor_Contrast=TempContrast, Direction="Down", Taxon=Dataset_Annotation$taxon, GeneSetSize=length(TempResultsJoined_ResultK_FoldChangeM_Down_Symbols), GeneSet=toString(TempResultsJoined_ResultK_FoldChangeM_Down_Symbols))
        
        HC_GeneSets_Temp<-rbind.data.frame(HC_GeneSets_Temp, TempResultsJoined_ResultK_FoldChangeM_DownDF)
        
        rm(TempResultsJoined_ResultK_FoldChangeM_Down, TempResultsJoined_ResultK_FoldChangeM_Down_Symbols, TempResultsJoined_ResultK_FoldChangeM_DownDF)
        
        }else{}
        
        #Pulling out up-regulated genes for each contrast:
        if(sum(TempResultsJoined_ResultK[,m]>0 & TempResultsJoined_ResultK[,m+2]<0.05, na.rm=TRUE)>0){
        TempResultsJoined_ResultK_FoldChangeM_Up<-TempResultsJoined_ResultK[TempResultsJoined_ResultK[,m]>0,]
        write.csv(TempResultsJoined_ResultK_FoldChangeM_Up, paste(colnames(TempResultsJoined)[k], "_10_w_", colnames(TempResultsJoined)[m], "_Up.csv", sep=""))
        
        if(nrow(TempResultsJoined_ResultK_FoldChangeM_Up)>999){
          TempResultsJoined_ResultK_FoldChangeM_Up_Reordered<-TempResultsJoined_ResultK_FoldChangeM_Up[order(TempResultsJoined_ResultK_FoldChangeM_Up[,m], decreasing=TRUE),]
          TempResultsJoined_ResultK_FoldChangeM_Up<-TempResultsJoined_ResultK_FoldChangeM_Up_Reordered[c(1:999),]}else{}
        
        TempResultsJoined_ResultK_FoldChangeM_Up_Symbols<-unique(TempResultsJoined_ResultK_FoldChangeM_Up$Gene_Symbol)
        
        write.csv(TempResultsJoined_ResultK_FoldChangeM_Up_Symbols, paste(colnames(TempResultsJoined)[k], "_10_w_", colnames(TempResultsJoined)[m], "_Up_UniqueSymbols.csv", sep=""))
        
      
        TempResultsJoined_ResultK_FoldChangeM_UpDF<-data.frame(Dataset_Name=Dataset_Name, Factor_Contrast=TempContrast, Direction="Up", Taxon=Dataset_Annotation$taxon, GeneSetSize=length(TempResultsJoined_ResultK_FoldChangeM_Up_Symbols), GeneSet=toString(TempResultsJoined_ResultK_FoldChangeM_Up_Symbols))
        
        HC_GeneSets_Temp<-rbind.data.frame(HC_GeneSets_Temp, TempResultsJoined_ResultK_FoldChangeM_UpDF)
        
        rm(TempResultsJoined_ResultK_FoldChangeM_Up, TempResultsJoined_ResultK_FoldChangeM_Up_Symbols, TempResultsJoined_ResultK_FoldChangeM_UpDF)
        
        }else{}
      }
    }else{}
    
  }
  
  rm(TempResultsJoined, TempIndicesForQValueColumns)
  }else{}
}

dim(HC_GeneSets_Temp)
#[1] 570   6
head(HC_GeneSets_Temp)
#It worked bwa ha ha ha

setwd("~/Documents/SideProjects/BrainGMT/Gemma")

write.table(HC_GeneSets_Temp, "HC_GeneSets_Temp.txt")

#Now we need to cut out all of the less informative gene sets (and/or more likely to be confounded comparisons):

grep("organism_part", HC_GeneSets_Temp$Factor_Contrast)
grep("sex", HC_GeneSets_Temp$Factor_Contrast)
grep("developmental_stage", HC_GeneSets_Temp$Factor_Contrast)

HC_GeneSets_Temp_Informative<-HC_GeneSets_Temp[-c(grep("organism_part", HC_GeneSets_Temp$Factor_Contrast), grep("sex", HC_GeneSets_Temp$Factor_Contrast), grep("developmental_stage", HC_GeneSets_Temp$Factor_Contrast)),]

dim(HC_GeneSets_Temp_Informative)
#[1] 469   6

sum(HC_GeneSets_Temp_Informative$GeneSetSize<1)
#[1] 0

sum(HC_GeneSets_Temp_Informative$GeneSetSize<10)
#[1] 128

write.table(HC_GeneSets_Temp_Informative, "HC_GeneSets_Temp_Informative.txt")

#Hmm... I just found another bug: Any time I have a folder with more than one AnalysisResults (e.g., GSE113796) it only parses by the q-values in the first file.
#Doesn't seem like this happened though.


###hmmmm.... I suppose I could compile the gene sets into one document and delete out the ones that we don't want too.
###I think I'll add that one in next time around - for now, I'd like to take a look and see how things went.


#28
#"13263_GSE86392_diffExpAnalysis_105657"
#all gene symbols are NAs (??!!!) - should probably note this on Gemma. Looks like the annotation just didn't output fully.
#Same with 65
#[1] "17271_GSE140598_diffExpAnalysis_167383"
#And 69:
#"17355_GSE140595_diffExpAnalysis_167639"



#Other steps left to do:
#1. If the results for the individual contrasts are much different from the overall ANOVA, consider making gene sets for both.
#3. Deal with "date genes"
#4. Count duplicates, and for any genes with duplicates showing the same direction of effect, choose the most extreme value. - I'll save dealing with this for once we are figuring out GeneWeaver.
#7. Output the final gene set sizes.


#Questions for Paul:
#1: Are they taking into account the multi-level nature of the data from multiple organism parts? (says ANOVA)
#2: Why do the contrast p-values for factors with only 2 levels not match the ANOVA output?
#3: How do we interpret the contrast output for models with interaction terms?  (main effect vs. effect within reference group?)

#########################################################

#Re-doing the same code for the NACC to make everything standard:
library(plyr)
library(dplyr)

setwd("~/Documents/SideProjects/BrainGMT/Gemma/NucleusAccumbens")

NACC_DatasetsOfInterest<-list.files()
#16

NACC_GeneSets_Temp<-data.frame(Dataset_Name=character(), Factor_Contrast=character(), Direction=character(), Taxon=character(), GeneSetSize=numeric(), GeneSet=character())

#Loops over all folders (publications) of down-loaded results:

#These folders are problematic because they completely lack gene symbols:


for(i in c(2:length(NACC_DatasetsOfInterest))){
  
  setwd("~/Documents/SideProjects/BrainGMT/Gemma/NucleusAccumbens")
  
  print(i)
  print(NACC_DatasetsOfInterest[i])
  
  #Pulling out relevant meta-data for the dataset analysis:
  Dataset_ID<-strsplit(NACC_DatasetsOfInterest[i], "_")[[1]]
  names(Dataset_ID)<-c("Gemma_ID", "GEO_Accession", "Output", "Result_ID")
  
  #This code bugs out when Gemma has added a .1 etc to denote more than one entry for a dataset, so I had to add a string split.
  
  Dataset_Annotation<-NACC_Annotation_AsDF_Simpler[which(NACC_Annotation_AsDF_Simpler$accession==strsplit(Dataset_ID[2][[1]], ".", fixed=TRUE)[[1]][1]),]
  
  Dataset_Name<-paste(Dataset_Annotation$accession, "_", gsub( " ", "_", Dataset_Annotation$name, fixed=TRUE), sep="")
  
  setwd(paste("./", NACC_DatasetsOfInterest[i], sep=""))
  
  TempResultsFiles<-list.files()
  
  #Pulls out the ANOVA results:
  TempAnalysisResults<-read.delim("analysis.results.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")
  
  TempResultsToJoin<-list(TempAnalysisResults)
  str(TempResultsToJoin)
  
  #Joins the ANOVA results to the results from the individual contrasts:
  for(j in c(2:length(TempResultsFiles))){
    print(k)
    TempResultsToJoin[[j]]<-assign(paste("TempResultSet", j, sep=""), read.delim(TempResultsFiles[j], sep="\t", stringsAsFactors = FALSE, comment.char = "#"))
  }
  
  TempResultsJoined<-join_all(TempResultsToJoin, by="Element_Name")
  
  write.csv(TempResultsJoined, "TempResultsJoined.csv")
  
  rm(TempResultsFiles, TempAnalysisResults, TempResultsToJoin)
  
  TempIndicesForQValueColumns<-grep("QValue", colnames(TempResultsJoined))
  TempIndicesForFoldChangeColumns<-grep("FoldChange", colnames(TempResultsJoined))
  
  #Removing rows lacking gene symbol annotation, with multi-mapped annotation, or "date" genes that need fixing:
  
  if(sum(is.na(TempResultsJoined$Gene_Symbol)==FALSE & (TempResultsJoined$Gene_Symbol=="")==FALSE & (TempResultsJoined$Gene_Symbol==" ")==FALSE)>0){
    
    TempResultsJoined_NoNA<-TempResultsJoined[is.na(TempResultsJoined$Gene_Symbol)==FALSE & (TempResultsJoined$Gene_Symbol=="")==FALSE & (TempResultsJoined$Gene_Symbol==" ")==FALSE & (TempResultsJoined$Gene_Symbol=="null")==FALSE & (TempResultsJoined$Gene_Symbol=="N/A")==FALSE,]
    
    rm(TempResultsJoined)
    
    TempResultsJoined_NoNA_NoMultimapped<-TempResultsJoined_NoNA[grep("[|]", TempResultsJoined_NoNA$Gene_Symbol, invert=TRUE),]
    
    rm(TempResultsJoined_NoNA)
    
    #To replace date genes, I need to first determine whether the gene symbols reflect humans or mice/rats, because the Excel perversion of the date genes is all lower case (i.e., the same in both species)
    #On the other hand, since this is all Gemma output, I'm pretty sure there won't be excel versions of the gene names. So maybe just skip that step and flag if there is an issue?
    
    sum(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol%in%DateGeneDecoder[,1])
    
    TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol<-recode(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, Sept1='Septin1',
                                                             Sept2='Septin2',
                                                             Sept3='Septin3',
                                                             Sept4='Septin4',
                                                             Sept5='Septin5',
                                                             Sept6='Septin6',
                                                             Sept7='Septin7',
                                                             Sept8='Septin8',
                                                             Sept9='Septin9',
                                                             Sept10='Septin10',
                                                             Sept11='Septin11',
                                                             Sept12='Septin12',
                                                             Sept13='Septin13',
                                                             Sept14='Septin14',
                                                             Sept7p2='Septin7p2',
                                                             March1='Marchf1',
                                                             March2='Marchf2',
                                                             March3='Marchf3',
                                                             March4='Marchf4',
                                                             March5='Marchf5',
                                                             March6='Marchf6',
                                                             March7='Marchf7',
                                                             March8='Marchf8',
                                                             March9='Marchf9',
                                                             March10='Marchf10',
                                                             March11='Marchf11',
                                                             Dec1='Delec1',
                                                             Nov='Ccn3')
    
    TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol<-recode(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, SEPT1='SEPTIN1',
                                                             SEPT2='SEPTIN2',
                                                             SEPT3='SEPTIN3',
                                                             SEPT4='SEPTIN4',
                                                             SEPT5='SEPTIN5',
                                                             SEPT6='SEPTIN6',
                                                             SEPT7='SEPTIN7',
                                                             SEPT8='SEPTIN8',
                                                             SEPT9='SEPTIN9',
                                                             SEPT10='SEPTIN10',
                                                             SEPT11='SEPTIN11',
                                                             SEPT12='SEPTIN12',
                                                             SEPT13='SEPTIN13',
                                                             SEPT14='SEPTIN14',
                                                             SEPT7P2='SEPTIN7P2',
                                                             MARCH1='MARCHF1',
                                                             MARCH2='MARCHF2',
                                                             MARCH3='MARCHF3',
                                                             MARCH4='MARCHF4',
                                                             MARCH5='MARCHF5',
                                                             MARCH6='MARCHF6',
                                                             MARCH7='MARCHF7',
                                                             MARCH8='MARCHF8',
                                                             MARCH9='MARCHF9',
                                                             MARCH10='MARCHF10',
                                                             MARCH11='MARCHF11',
                                                             DEC1='DELEC1',
                                                             NOV='CCN3')
    
    sum(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol%in%DateGeneDecoder[,1])
    
    write.csv(TempResultsJoined_NoNA_NoMultimapped, "TempResultsJoined_GeneSymbolsClean.csv")
    
    #Overwrite the other version of the df to reuse previously-written code:
    TempResultsJoined<-TempResultsJoined_NoNA_NoMultimapped
    
    rm(TempResultsJoined_NoNA_NoMultimapped)
    
    #Loops over each variable and pulls out the rows of genes with q<0.10 (FDR<0.10) and p<0.0001 (to help cut down on the # of genes in results files with thousands of "sig" results):
    
    for(k in c(TempIndicesForQValueColumns)){
      print(k)
      print(colnames(TempResultsJoined)[k])
      if(sum(TempResultsJoined[,k]<0.10 & TempResultsJoined[,k+1]<0.0001, na.rm=TRUE)>0){
        TempResultsJoined_ResultK<-TempResultsJoined[TempResultsJoined[,k]<0.10 & TempResultsJoined[,k+1]<0.0001,]
        write.csv(TempResultsJoined_ResultK, paste(colnames(TempResultsJoined)[k], "_10.csv", sep=""))
        
        #For the genes that show a significant relationship with a variable (q<0.10 & p<0.0001), this code loops over each level of contrast and pulls out the results that are p<0.05 for that contrast. 
        #Note that only a handful of these combinations are actually meaningful - I couldn't figure out a good way to automatically match variable + relevant contrasts.
        
        for(m in c(TempIndicesForFoldChangeColumns)){
          print(m)
          
          TempContrast<-paste(gsub("QValue_", "", colnames(TempResultsJoined)[k]), gsub("FoldChange_", "", colnames(TempResultsJoined)[m]), sep="_")
          print(TempContrast)
          
          #Pulling out down-regulated genes for each contrast:
          if(sum(TempResultsJoined_ResultK[,m]<0 & TempResultsJoined_ResultK[,m+2]<0.05, na.rm=TRUE)>0){
            TempResultsJoined_ResultK_FoldChangeM_Down<-TempResultsJoined_ResultK[TempResultsJoined_ResultK[,m]<0 & TempResultsJoined_ResultK[,m+2]<0.05,]
            write.csv(TempResultsJoined_ResultK_FoldChangeM_Down, paste(colnames(TempResultsJoined)[k], "_10_w_", colnames(TempResultsJoined)[m], "_Down.csv", sep=""))
            
            if(nrow(TempResultsJoined_ResultK_FoldChangeM_Down)>999){
              TempResultsJoined_ResultK_FoldChangeM_Down_Reordered<-TempResultsJoined_ResultK_FoldChangeM_Down[order(TempResultsJoined_ResultK_FoldChangeM_Down[,m]),]
              TempResultsJoined_ResultK_FoldChangeM_Down<-TempResultsJoined_ResultK_FoldChangeM_Down_Reordered[c(1:999),]
              rm(TempResultsJoined_ResultK_FoldChangeM_Down_Reordered)
            }else{}
            
            TempResultsJoined_ResultK_FoldChangeM_Down_Symbols<-unique(TempResultsJoined_ResultK_FoldChangeM_Down$Gene_Symbol)
            
            write.csv(TempResultsJoined_ResultK_FoldChangeM_Down_Symbols, paste(colnames(TempResultsJoined)[k], "_10_w_", colnames(TempResultsJoined)[m], "_Down_UniqueSymbols.csv", sep=""))
            
            TempResultsJoined_ResultK_FoldChangeM_DownDF<-data.frame(Dataset_Name=Dataset_Name, Factor_Contrast=TempContrast, Direction="Down", Taxon=Dataset_Annotation$taxon, GeneSetSize=length(TempResultsJoined_ResultK_FoldChangeM_Down_Symbols), GeneSet=toString(TempResultsJoined_ResultK_FoldChangeM_Down_Symbols))
            
            NACC_GeneSets_Temp<-rbind.data.frame(NACC_GeneSets_Temp, TempResultsJoined_ResultK_FoldChangeM_DownDF)
            
            rm(TempResultsJoined_ResultK_FoldChangeM_Down, TempResultsJoined_ResultK_FoldChangeM_Down_Symbols, TempResultsJoined_ResultK_FoldChangeM_DownDF)
            
          }else{}
          
          #Pulling out up-regulated genes for each contrast:
          if(sum(TempResultsJoined_ResultK[,m]>0 & TempResultsJoined_ResultK[,m+2]<0.05, na.rm=TRUE)>0){
            TempResultsJoined_ResultK_FoldChangeM_Up<-TempResultsJoined_ResultK[TempResultsJoined_ResultK[,m]>0,]
            write.csv(TempResultsJoined_ResultK_FoldChangeM_Up, paste(colnames(TempResultsJoined)[k], "_10_w_", colnames(TempResultsJoined)[m], "_Up.csv", sep=""))
            
            if(nrow(TempResultsJoined_ResultK_FoldChangeM_Up)>999){
              TempResultsJoined_ResultK_FoldChangeM_Up_Reordered<-TempResultsJoined_ResultK_FoldChangeM_Up[order(TempResultsJoined_ResultK_FoldChangeM_Up[,m], decreasing=TRUE),]
              TempResultsJoined_ResultK_FoldChangeM_Up<-TempResultsJoined_ResultK_FoldChangeM_Up_Reordered[c(1:999),]}else{}
            
            TempResultsJoined_ResultK_FoldChangeM_Up_Symbols<-unique(TempResultsJoined_ResultK_FoldChangeM_Up$Gene_Symbol)
            
            write.csv(TempResultsJoined_ResultK_FoldChangeM_Up_Symbols, paste(colnames(TempResultsJoined)[k], "_10_w_", colnames(TempResultsJoined)[m], "_Up_UniqueSymbols.csv", sep=""))
            
            
            TempResultsJoined_ResultK_FoldChangeM_UpDF<-data.frame(Dataset_Name=Dataset_Name, Factor_Contrast=TempContrast, Direction="Up", Taxon=Dataset_Annotation$taxon, GeneSetSize=length(TempResultsJoined_ResultK_FoldChangeM_Up_Symbols), GeneSet=toString(TempResultsJoined_ResultK_FoldChangeM_Up_Symbols))
            
            NACC_GeneSets_Temp<-rbind.data.frame(NACC_GeneSets_Temp, TempResultsJoined_ResultK_FoldChangeM_UpDF)
            
            rm(TempResultsJoined_ResultK_FoldChangeM_Up, TempResultsJoined_ResultK_FoldChangeM_Up_Symbols, TempResultsJoined_ResultK_FoldChangeM_UpDF)
            
          }else{}
        }
      }else{}
      
    }
    
    rm(TempResultsJoined, TempIndicesForQValueColumns)
  }else{}
}

dim(NACC_GeneSets_Temp)
#[1] 196   6
head(NACC_GeneSets_Temp)
#It worked bwa ha ha ha

setwd("~/Documents/SideProjects/BrainGMT/Gemma")

write.table(NACC_GeneSets_Temp, "NACC_GeneSets_Temp.txt")

#Now we need to cut out all of the less informative gene sets (and/or more likely to be confounded comparisons):

grep("organism_part", NACC_GeneSets_Temp$Factor_Contrast)
grep("sex", NACC_GeneSets_Temp$Factor_Contrast)
grep("developmental_stage", NACC_GeneSets_Temp$Factor_Contrast)

NACC_GeneSets_Temp_Informative<-NACC_GeneSets_Temp[-c(grep("organism_part", NACC_GeneSets_Temp$Factor_Contrast), grep("sex", NACC_GeneSets_Temp$Factor_Contrast), grep("developmental_stage", NACC_GeneSets_Temp$Factor_Contrast)),]

dim(NACC_GeneSets_Temp_Informative)
#[1] 139   6

sum(NACC_GeneSets_Temp_Informative$GeneSetSize<1)
#[1] 0

sum(NACC_GeneSets_Temp_Informative$GeneSetSize<10)
#[1] 36

write.table(NACC_GeneSets_Temp_Informative, "NACC_GeneSets_Temp_Informative.txt")




#########################################################

#Trying to get GeneWeaver API to work:




#my API key is in this code - for public release I've replaced it with [MyAPIKey]

https://geneweaver.org/api/get/publication/byid/[MyAPIKey]/26/
  
  https://geneweaver.org/api/get/geneset/bygenesetid/8/
  #Not found
  #adding geneset->genesets doesn't help
  
  https://geneweaver.org/api/get/geneset/bygenesetid/719/
  
  
  #More documentation:
  https://geneweaver.org/api/swagger.json

https://geneweaver.org/api/genesets/get/
  
  https://geneweaver.org/api/genesets/bygenesetid/8/
  #Nope
  
  https://geneweaver.org/api/genesets/get/bygenesetid/8/
  
  https://geneweaver.org/api/genesets/8
#This brings me to the results management page (??) - changing out the number for other numbers does the same...


#Trying things out now that the api is supposedly fixed...



#This one seems to work (it pulls up metadata for the gene set,  not the gene set itself):  

  #Get Gene Set by Gene Set ID: This call returns all information about a specified gene set given that gene set ID.
  #/api/get/geneset/byid/<GeneSetID>/
  #Sample Call: https://geneweaver.org/api/get/geneset/byid/220592/
  #my attempt - using Gene set id for EtOH Preference Nucleus accumbens:
  https://geneweaver.org/api/get/geneset/byid/719/
  #seems to work in chrome - but doesn't include pubmed id :(
  
  #It says the associated publication is pub_id: 26, so let's try feeding that into:
    
  #Get Publication by Publication ID: Returns all the publication data for given publication ID.
  
  #/api/get/publication/byid/<apikey>/<PublicationID>/
    #Sample Call: https://geneweaver.org/api/get/publication/byid/Fw7J4GeAXE8CMVvLTKyrtBDk/26/
    
    https://geneweaver.org/api/get/publication/byid/[MyAPIKey]/26/
    
    #Includes full publication citation information including:
    "pub_year": "2007"
    "pub_pubmed": "17451403"
    
  #Crazy thought: since these pub_ids seem internally generated, can I just loop up through them until they stop?
    
    https://geneweaver.org/api/get/publication/byid/[MyAPIKey]/1/ #empty
    https://geneweaver.org/api/get/publication/byid/[MyAPIKey]/2/ #for gene set without publication?
      https://geneweaver.org/api/get/publication/byid/[MyAPIKey]/3/
      #It seems like everything above 8470 (up until around 100,000) is scrap or GO Ontology
      
      #Can I adapt the API to get the geneset id from the publications? (are all publications paired with a geneset?)
      
      

      
  #Get Geneset by Geneset ID: This call returns all the information about a given geneset given its geneset ID
  #/api/get/geneset/bygenesetid/<GeneSetID>/
  Sample Call: https://geneweaver.org/api/get/geneset/bygenesetid/8/
  #This sample call doesn't work 
  
  #Get Genes by Gene Set ID: This call returns all genes belonging to a given gene set.
  /api/get/genes/bygenesetid/<GeneSetID>/
  /api/get/genes/bygenesetid/<GeneSetID>/
  https://geneweaver.org/api/get/genes/bygenesetid/220592/
  
  https://geneweaver.org/api/get/genes/bygenesetid/719/
  #The gene symbols seem to be in "ode_ref_id":
  #It looks like the genes displayed as part of the geneset are "ode_pref": true
  
  https://geneweaver.org/api/get/genes/bygenesetid/166586/
  #Same pattern
  
  #Get Projects by User: Returns all the projects that are owned by a given user.
  #/api/get/project/byuser/<apikey>/
    
    https://geneweaver.org/api/get/project/byuser/[MyAPIKey]/
    #Works... yess
    
    #So let's try pulling up the gene sets for those projects:
   
     #Nucleus Accumbens project:
    https://geneweaver.org/api/get/geneset/byprojectid/[MyAPIKey]/32017/
    #This brings me to the results management page (??) - changing out the number for other numbers does the same...
    
    https://geneweaver.org/api/get/geneset/fromproject/[MyAPIKey]/32017/
      
      
      ########################
    
    #Publication information from ISI:
    
    setwd("~/Documents/SideProjects/BrainGMT/ISI")
    
    BrainMicroarrayPubs<-read.csv("BrainMicroarray_savedrecs_CompiledwCitationReport_forR.csv", header=TRUE, stringsAsFactors = FALSE)  
    str
    
    #Adding some filters:
    dim(BrainMicroarrayPubs[BrainMicroarrayPubs$Publication.Year<2012 & BrainMicroarrayPubs$Publication.Year>2000,])
    #[1] 3435   50
    
    table(BrainMicroarrayPubs$Publication.Year, BrainMicroarrayPubs$Associated.Data...T.F.)
    #        0   1
    # 1996   1   0
    # 1998   2   0
    # 1999   8   0
    # 2000  27   0
    # 2001  77   0
    # 2002 120   8
    # 2003 162   9
    # 2004 249  18
    # 2005 300  44
    # 2006 302  48
    # 2007 295  58
    # 2008 340  79
    # 2009 337  89
    # 2010 333  82
    # 2011 378  79
    # 2012 353 124
    # 2013 307 113
    # 2014 344 104
    # 2015 365  75
    # 2016 384  25
    # 2017 378  22
    # 2018 326  18
    # 2019 308  10
    # 2020 272  18
    # 2021 125   1
    
    #Unlike Pubmed, it doesn't seem like most publications have associated data after 2011. I suspect they are being massively undercounted.
    #I wonder why the numbers drop after 2014
    #Join with Pubmed output to get better quality data?
    
    ##############
    
    setwd("~/Documents/SideProjects/BrainGMT/Pubmed")
    
    Pubmed_PubsWData<-read.csv("BrainMicroarrayWData_PubsCompiled.csv", header=TRUE, stringsAsFactors = FALSE)
    str(Pubmed_PubsWData)
    
    head(BrainMicroarrayPubs$Pubmed.Id)
    
    sum(BrainMicroarrayPubs$Pubmed.Id%in%Pubmed_PubsWData$PMID)
    #[1] 2116
    
    BrainMicroarrayPubs$Has.Associated.Data.InPubmed<-BrainMicroarrayPubs$Pubmed.Id%in%Pubmed_PubsWData$PMID
    
    table(BrainMicroarrayPubs$Associated.Data...T.F., BrainMicroarrayPubs$Has.Associated.Data.InPubmed)
    #   FALSE TRUE
    # 0  4431 1690
    # 1   598  426
    
    #Interesting- each dataset has a slightly different set of pubs that are indicated to have associated data.
    
    sum(BrainMicroarrayPubs$Associated.Data...T.F.==TRUE|BrainMicroarrayPubs$Has.Associated.Data.InPubmed==1)
    #[1] 2714
    
    dim(BrainMicroarrayPubs)
    #[1] 7145   51
    
    #I wonder how much the discrepancy is explained by Publications lacking PubmedIDs
    
    BrainMicroarrayPubs_wPMID<-BrainMicroarrayPubs[is.na(BrainMicroarrayPubs$Pubmed.Id)==FALSE,]
    
    dim(BrainMicroarrayPubs_wPMID)
    #[1] 6921   51
    
    table(BrainMicroarrayPubs_wPMID$Associated.Data...T.F., BrainMicroarrayPubs_wPMID$Has.Associated.Data.InPubmed)
    #   FALSE TRUE
    # 0  4209 1690
    # 1   596  426
    
    #That mostly just eliminated publications from the "no associated data" category.
    
    #I wonder if it would be an even higher percentage if we only looked at particular years (and maybe with some filtering criteria that might generally indicate publication quality)
    
    table((BrainMicroarrayPubs_wPMID$Associated.Data...T.F.| BrainMicroarrayPubs_wPMID$Has.Associated.Data.InPubmed), BrainMicroarrayPubs_wPMID$Publication.Year)
    
    #       1996 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021
    # FALSE    1    2    8   26   69  114  157  245  286  283  280  323  320  317  358  127   70   71  194  111   70  105  276  253  117
    # TRUE     0    0    0    1    6   11   10   19   48   56   60   82   91   84   85  333  340  369  232  281  311  228   33   24    7
    
  #I wonder how many of these also have PMC IDs (i.e., NIH funding)
  #I'll have to crossreference it with the Pubmed database to gather that information
    
    
  
    ###############
    
    #What's in GeneWeaver:
    
    setwd("~/Documents/SideProjects/BrainGMT/GeneWeaver")
    GeneWeaver_GeneSets<-read.delim("publications2symbols.dat", sep="\t", header=FALSE, stringsAsFactors = FALSE)
    str(GeneWeaver_GeneSets)
    head(GeneWeaver_GeneSets)
    
    # 'data.frame':	693 obs. of  4 variables:
    #   $ V1: chr  "Differential effects in hippocampus and amygdala of high and low anxiety selected lines" "Hippocampus vs. Amygdala main effect Contextual fear conditioning Selected Lines -Differential Expression" "Ethanol preference related in Amygdala of congenic rat strains" "EtOH Preference related in Hippocampus of Congenic Rat Strains" ...
    # $ V2: num  17309658 17309658 17451403 17451403 17451403 ...
    # $ V3: int  10090 10090 10116 10116 10116 10116 10116 10090 10116 10116 ...
    # $ V4: chr  "Igfals,Polr1b,Noct,Col9a3,Ptprq,Ybx2,Tex261,Nrg3,Ackr1,Klrc2,Pak3,Galnt7,Zw10,Sdhaf3,Myo19,Tmem33,Akirin1,G6pc3"| __truncated__ "Apod,Cd19,Hspa9-ps1,Rab7,Igfals,Polr1b,Noct,Col9a3,Ptprq,Ybx2,Tex261,Nrg3,Ackr1,Klrc2,Pak3,Galnt7,Zw10,Sdhaf3,M"| __truncated__ "" "" ...
    
    #Can I write these out into a file type that is readable in Excel?

    write.table(GeneWeaver_GeneSets, "publications2symbols.txt", sep="\t")
    
    #Sanity Check:
    GeneWeaver_GeneSets[GeneWeaver_GeneSets$V1=="Differentially expressed genes upregulated in Amygdala after single-prolonged stress treatment",]
    #PSMB8,RHBDF1,Rtl8a,LSM3P5,LOC100361572
    #On GeneWeaver website: Ace, Htr2c, Kl, Nfx1, Otx2
    #Website matches pubication.
    
    GeneWeaver_GeneSets[GeneWeaver_GeneSets$V1=="Downregulation of gene expression in the lateral hypothalamus of mu opioid Knockout (KO) mice following administration of chronic morphine",]
    #Frs2,Gpr33,Efcab1,Cibar1
    #On GeneWeaver website: Apln, Aqp4, RGS16, RGS4
    #Website matches pubication.
    
    GeneWeaver_GeneSets[GeneWeaver_GeneSets$V1=="Upregulated Nicotine-dependent genes in human neuroblastoma (SH-SY5Y) cells",]
    #Gm44226,Gm52104,PDYN
    #On GeneWeaver website: RFXAP, TFPI2, ZFR
    
    
    #One of those columns of numbers is PubmedIDs, the other is probably taxon (maybe?)
    dim(table(GeneWeaver_GeneSets$V2))
    #[1] 298
    #Looks like there are 298 pubmed ids associated with the 693 gene sets.
    
    sum(GeneWeaver_GeneSets$V2%in%BrainMicroarrayPubs$Pubmed.Id)
    #[1] 250
    #Nice...
    
    table(GeneWeaver_GeneSets$V3)
    # 4932  6239  7227  7955  9606 10090 10116 
    # 4    10     8     4   296   254   117
    #These must be taxon
    
    GeneWeaver_GeneSets$V4[1]
    # [1] "Igfals,Polr1b,Noct,Col9a3,Ptprq,Ybx2,Tex261,Nrg3,Ackr1,Klrc2,Pak3,Galnt7,Zw10,Sdhaf3,Myo19,Tmem33,Akirin1,G6pc3,Tmem140,Lce1b,Dipk2a,Mtfp1,1700027A15Rik,Tspo2,Sdr9c7,Zfp995,Cyp4f16,Adgra3,Mex3b,Cenpu,Noxo1,4932437C15Rik,Atp13a2,4930597L12Rik,4833422M21Rik,Dcp1a,Rab11fip4os2,Efcab11,Tnfrsf23,AI480526,Arrdc3,Xxylt1,Nsun3,Ash1l,Plch2,Jph4,Rccd1,Dhtkd1,Hspa9-ps1,Rab7,Apod,Cd19"
    
    GeneWeaver_GeneSets$V4[2]
    #[1] "Apod,Cd19,Hspa9-ps1,Rab7,Igfals,Polr1b,Noct,Col9a3,Ptprq,Ybx2,Tex261,Nrg3,Ackr1,Klrc2,Pak3,Galnt7,Zw10,Sdhaf3,Myo19,Tmem33,Akirin1,G6pc3,Tmem140,Lce1b,Dipk2a,Mtfp1,1700027A15Rik,Tspo2,Sdr9c7,Zfp995,Cyp4f16,Adgra3,Mex3b,Cenpu,Noxo1,4932437C15Rik,Atp13a2,4930597L12Rik,4833422M21Rik,Dcp1a,Rab11fip4os2,Efcab11,Tnfrsf23,AI480526,Arrdc3,Xxylt1,Nsun3,Ash1l,Plch2,Jph4,Rccd1,Dhtkd1"
    
    #Looks great. :)
    
   
    table(BrainMicroarrayPubs_wPMID$Pubmed.Id%in%GeneWeaver_GeneSets$V2, (BrainMicroarrayPubs_wPMID$Associated.Data...T.F.| BrainMicroarrayPubs_wPMID$Has.Associated.Data.InPubmed))
    
    #       FALSE TRUE
    # FALSE  4164 2682
    # TRUE     45   30
    
    #Huh- why does that only add up to 75?
    
    table(BrainMicroarrayPubs_wPMID$Pubmed.Id%in%GeneWeaver_GeneSets$V2)
    
    # FALSE  TRUE 
    # 6846    75 
    
    #Oh wait - there must be only a few Pubmed IDs mapping to lots of gene sets in GeneWeaver.
    
    length(GeneWeaver_GeneSets$V2[GeneWeaver_GeneSets$V2%in%BrainMicroarrayPubs_wPMID$Pubmed.Id])
    #[1] 245
    
    dim(table(GeneWeaver_GeneSets$V2[GeneWeaver_GeneSets$V2%in%BrainMicroarrayPubs_wPMID$Pubmed.Id]))
    #[1] 75
    #There we go
    
    #Double-checking against the larger Pubmed database:
    
    sum(GeneWeaver_GeneSets$V2%in%Pubmed_PubsWData$PMID)
    #[1] 52
    
    dim(table(GeneWeaver_GeneSets$V2))
    #[1] 298
    
    #So roughly 3/4 of the GeneWeaver gene sets are either not explicitly brain related or don't have associated data in Pubmed. 
    #Since there are so few, I should probably just glance over the list of gene sets and see which look like they're not brain - I need to do that anyway for my .gmt.
    #Eyeballing the list, there are a few related to ulcerative colitis or GWAS candidates, but mostly brain transcriptional profiling. There are also some that are empty (no genes associated with the gene set)
    #Ah - I forgot that these publications can also be RNA-Seq - although the main question I was trying to answer was regarding how many of the old microarray studies are covered by GeneWeaver.
    #Interesting - most of the gene sets in GeneWeaver aren't directional.
    
    #When I come back:
    #Prune the GeneWeaver .gmt file down to brain transcriptional profiling studies
    #Try plotting some of the Publications w/ data vs. without for citation counts 1-2 years post publication and 5-6 years post publication
    #Finish setting up files for Evelyn to make boxplots (joining QC-ed subject info with drug info, outputting DeltaDeltaCq for GabaGlu and DA5HT)
    
    
    
    