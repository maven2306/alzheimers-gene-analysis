# MOST SIGNIFICANT GENES IN LATE-ONSET ALZHEIMER’S DISEASE
<div align="center">
 <h2>This project was conducted collaboratively by Matteo Venturini and Ana Esteve Garcia.</h2>
</div> 
#
Alzheimer's disease (AD) is marked by β-amyloid plaques and tau protein neurofibrillary tangles. Genetic factors play a significant role in AD, with polymorphisms in the TOMM40, APOE, PVRL2, and APOC1 genes linked to increased risk. TOMM40 is crucial for mitochondrial function and is downregulated in AD, with a poly T repeat associated with the disease's onset. APOE, particularly the ε4 allele, is the strongest genetic risk factor, influencing neurofibrillary tangles and amyloid-peptide aggregates. PVRL2, involved in adherence junctions, may interact with TOMM40 to impact Aβ metabolism, increasing AD risk. APOC1, responsible for cholesterol metabolism and neuronal processes, has an insertion allele linked to AD. This report presents evidence of significant expression differences in these genes between AD patients and healthy controls, suggesting they are key targets for drug development to treat AD.

## GWAS ANALYSIS
We first conducted a GWAS analysis. Genome-Wide-Association-Studies are helpful tools to analyse which
mutations in which SNPs are correlated to a given phenotype or disease. Out of 18950 genes analysed, we found 44
genes whose frequency was significantly different in AD patients compared to the healthy control.

<div align="center">
 <img width="407" alt="image" src="https://github.com/maven2306/genes_Alzheimer-s/assets/169473359/498e02d3-2d33-40f8-9e3a-7e64aaf63475">
</div> 

<div align="center">
 Figure1. Manhattanplot. The X-axis shows genes found on each chromosome, and the Y-axis shows the log of the
 P-value. All the genes above the significance threshold line are correlated with AD.
</div> 

# 
Out of the 44, we selected 4 genes that had the lowest p-value and hence the ones that were more correlated to AD. Our genes of interest were all on Chromosome 19.
<div align="center">
 <img width="536" alt="image" src="https://github.com/maven2306/genes_Alzheimer-s/assets/169473359/3efc415a-71db-47c4-b719-519ba382fc83">
</div> 

<div align="center">
 Figure 2. This table shows the 44 most significant genes of the GWAS analysis, along with their p-value and on which chromosome they're found. The 4 most significant ones are highlighted.
</div> 

##  TRANSCRIPTOME ANALYSIS
We then ran a transcriptome analysis to investigate the gene expression level of each of these 4 genes in each of the following brain regions: Entorhinal Cortex, Hippocampus, Middle Temporal Gyrus, Posterior Cingulate, and Superior Frontal Gyrus. We compared the gene expression of these genes between patients and controls, and we made a violin plot that shows in which region each of our selected genes is significantly differentially expressed.

<div align="center">
 <img width="701" alt="image" src="https://github.com/maven2306/genes_Alzheimer-s/assets/169473359/963c5c2f-d9f6-4612-a8ba-77ccbc0b27fc">
</div> 
<div align="center">
 Figure 3. Four violin plots show in which region each gene’s expression pattern is statistically different between patients and controls. TOMM40 shows a statistically different pattern in SFG, APOC1 in PC, PVRL2 in HIP, APOE in MTG.
</div> 

 ## PPI-NETWORKS
Afterwards, we looked for protein-protein interactions between the genes found in the GWAS analysis. This gave us a
 broader perspective and understanding of how these proteins work, and how these interactions and associations
 between each other could be one of the targets for future research in Alzheimer's Disease. To do this we used the
 datasets given by Sügis et al, 2019. We found 128 genes that interact with our 44 significant AD genes found in the
 GWAS. One significant AD gene we focused on was APOE, and to determine its predicted protein interaction
 partners, we used the database STRING. In the figure below, the known interactions from curated databases and the
 experimentally determined ones are shown. They are all proteins that interact with APOE affecting Alzheimer’s.

