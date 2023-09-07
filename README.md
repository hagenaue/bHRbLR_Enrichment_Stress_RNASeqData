# Code Repository: bHRbLR_Enrichment_Stress_RNASeqData

This repository includes the code for the analysis of RNA-Seq data from the nucleus accumbens (NACC) and hippocampus (HC) within the publication "Adolescent social and environmental enrichment produces long-term increases in social resilience in a genetic rodent model of stress vulnerability: Impact on behavior, circulating hormones, and brain gene expression" by *Angela M. O'Connor, *Megan H. Hagenauer, Liam Cannon Thew Forrester, Pamela Patterson, Keiko Arakawa, Huzefa Khalil, Evelyn R. Richardson, Elaine K. Hebda-Bauer, Farizah I. Rob, Yusra Sannah, Stanley J. Watson, Jr., Huda Akil. (DOI: TBA)

This repository was created by Megan H. Hagenauer (ORCID: 0000-0002-3715-9475). This work was completed at the University of Michigan between 10/2019-09/2023.

## Input:

The input subject metadata for the datasets used in this code can be found in the respository here: 
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/Angela_Sample_MetaData_Updated20200910.txt

The input nucleus accumbens (NACC) RNA-Seq data (count matrix) used in this code can be found in the respository here: 
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/NAcc_remove2samples_gene_featureCounts_counts.txt

The input hippocampal (HC) RNA-Seq data (count matrix) used in this code can be found in the respository here: 
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/Hippocampus_remove6samples_gene_featureCounts_counts.txt

The input gene annotation for the RNA-Seq data can be found in the repository here:
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/Annotation_Fan.txt

The raw RNA-Seq data has also been uploaded to SRA/GEO (DOI: TBA) with detailed metadata.


## Analysis Code:

The code for the primary analyses used in the paper is located in a single file per brain region:

Nucleus Accumbens (NACC) Analysis Code:
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/Angela_NACC_20210129_ForCodeRelease.R

Hippocampal (HC) Analysis Code:
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/Angela_HC_20210129_ForCodeRelease.R

Code for exploratory analyses directly examining the relationship between gene expression and behavior or plasma hormone levels (instead of treatment group). These results were not included in the paper:
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/Angela_QuickAnalysis_HC_BehavHormones_20210129.R
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/Angela_QuickAnalysis_NACC_HormonesBehav_20210129.R

The code for the construction of the custom gene set database (Brain.GMT) is located here for the purposes of manuscript review - it is going to be submitted as its own publication:
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/Code_ImprovingGMTforSTRHC_forCodeRelease.R
https://github.com/hagenaue/bHRbLR_Enrichment_Stress_RNASeqData/blob/bbcb5c436b09afda2ac462e8d6fec4b9b047fd73/Code_Gemma_NACC_HC_20210603_forRelease.R

## Notes/Warnings:

There are several naming conditions in the code file that differ from what was used in the final paper (and GEO/SRA data release):

Adolescent Enrichment: EE=Social and Environmental ("Enhanced) Enrichment - sometimes just referred to as environmental enrichment in the code SE=Social Enrichment - referred to as EC or "cage enrichment" in the code NIL=Standard Housing

Bred Line: bLR=bred Low Responder rat line - sometimes just called LR in the code bHR=bred High Responder rat line - sometimes just called HR in the code


All analyses should be ignored that include ultrasonic vocalization data (we eventually realized this data was recorded incorrectly - 10x differences in gain) or time in the center of the open field (we realized this data cannot be interpreted as usual in terms of exploration/anxiety because we used a non-standard protocol in which the rats were placed in the center at start of the test causing anxious animals to sometimes freeze in the center instead of showing traditional thigmotaxis).
