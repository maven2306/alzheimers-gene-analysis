<div align="center">
 <h1>MOST SIGNIFICANT GENES IN LATE-ONSET ALZHEIMER’S DISEASE</h1>
</div> 

<div align="center">
 <h2>This project was conducted collaboratively by Matteo Venturini and Ana Esteve Garcia.</h2>
</div> 


### Table of Contents

1. [GWAS ANALYSIS](#GWAS-ANALYSIS)
2. [TRANSCRIPTOME ANALYSIS](#TRANSCRIPTOME-ANALYSIS)
3. [PPI-NETWORKS](#PPI-NETWORKS)

---

Alzheimer's disease (AD) is marked by β-amyloid plaques and tau protein neurofibrillary tangles. Genetic factors play a significant role in AD, with polymorphisms in the TOMM40, APOE, PVRL2, and APOC1 genes linked to increased risk. TOMM40 is crucial for mitochondrial function and is downregulated in AD, with a poly T repeat associated with the disease's onset. APOE, particularly the ε4 allele, is the strongest genetic risk factor, influencing neurofibrillary tangles and amyloid-peptide aggregates. PVRL2, involved in adherence junctions, may interact with TOMM40 to impact Aβ metabolism, increasing AD risk. APOC1, responsible for cholesterol metabolism and neuronal processes, has an insertion allele linked to AD. This report presents evidence of significant expression differences in these genes between AD patients and healthy controls, suggesting they are key targets for drug development to treat AD.

## GWAS ANALYSIS
We first conducted a GWAS analysis. Genome-Wide-Association-Studies are helpful tools to analyse which
mutations in which SNPs are correlated to a given phenotype or disease. Out of 18950 genes analysed, we found 44
genes whose frequency was significantly different in AD patients compared to the healthy control.

<div align="center">
 <img width="407" alt="image" src="https://github.com/maven2306/genes_Alzheimer-s/assets/169473359/498e02d3-2d33-40f8-9e3a-7e64aaf63475">
</div> 

<div align="center">
 Figure 1. Manhattanplot. The X-axis shows genes found on each chromosome, and the Y-axis shows the log of the
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

## AHB EXPRESSION ANALYSIS
APOE is the most important risk factor for developing late-onset Alzheimer’s disease. That’s why we decided to take
 a closer look at which cortical regions it is more expressed. We decided to select the 10 regions that show the
 highest expression levels of this gene. We then decided to plot them in a cortical map for better visualisation. The light green areas are the regions where
 APOE is more expressed. We used the data from only the left hemisphere because the data from the right
 hemisphere was not complete. We can see how APOE is mostly expressed in the entorhinal cortical region. 
 
<div align="center">
 <img width="709" alt="image" src="https://github.com/maven2306/genes_Alzheimer-s/assets/169473359/94ee51aa-923b-44ab-abd6-6b276d2ec38c">
</div> 
<div align="center">
  Figure 4. We can see a cortical map of the left hemisphere and the right hemisphere. Given that we used data from
  only the left hemisphere, the right hemisphere doesn’t show any differential colours. Instead, we can see how in the
  left hemisphere APOE is mostly expressed in the entorhinal cortical region
</div> 

##  MRI ANALYSIS
 Some very interesting insights were obtained using data from MRI scans of patients with AD and controls. Here we
 can see the most shrunk cortical regions in patients affected with AD. The effect size is the t-value obtained from
 t-tests ran to test for a difference in the brain volume between controls and patients. The most shrunk cortical region
 is the parahippocampal region in the left hemisphere, which is also confirmed by scientific literature (Thompson et
 al., 2003), followed by the isthumscingulate in the left hemisphere and the parsobitalis also in the left hemisphere.
 We also performed a statistical evaluation, most specifically a regression between APOE expression in the left
 hemisphere of the brain with the pathology score of AD in the left and right hemispheres of the brain. The result of
 the test was 4.4457e-25, which is extremely significant. This suggests a very strong and meaningful relationship
 between APOE expression and the pathology scores of AD.
 
<div align="center">
  <img width="820" alt="image" src="https://github.com/maven2306/genes_Alzheimer-s/assets/169473359/e2c1e174-db52-4dc8-ab36-47e8567cb849">
</div> 
<div align="center">
 Figure 5. In this map of the cortical regions, we can see highlighted the most shrunk ones. The darker is the colour, the most shrunk the region is. 
</div> 

# 
 Furthermore, to show our resulting association between the APOE gene expression and its effect on brain volume, a
 scatter plot with the best-fit line and with a 95% confidence interval was made. This, presented in Figure 6, shows
 how a decrease in APOE expression leads to a bigger effect size on the volume of the brain.
<div align="center">
 <img width="533" alt="image" src="https://github.com/maven2306/genes_Alzheimer-s/assets/169473359/8c1be862-7eb4-48b3-9d32-4885cdde9237">
</div> 
<div align="center">
 Figure 6. Correlation between the normal APOE expression
  of cortical regions and the effect size difference in the brain
  volume of Alzheimer’s patients, with a confidence interval of
  95%.
</div> 

## MACHINE LEARNING
What if we had ascript that could determine whether someone will develop Late-Onset AD based on their MRI data?
 Well, that’s what we did. For this, we used two different models. In the first one, we split our groups of patients and
 controls into train and test sets. Then the classifier was assessed and its effectiveness was evaluated. The trained
 classifier showed a promising performance with only getting incorrect 1.67% of the cases in the training set. However,
 it also had 9 false positives and 2 false negatives. We could conclude that the performance of this model was very
 accurate and it would work to identify if a new participant is a potential patient. But further optimization of the
 classifier could be done to minimise the false positives and the overall performance. We then made a second
 prediction model, using a neural network. This model showed fewer false negatives, only around 2, while more false
 positives, around 30. By changing the setting parameters we got an 87.56% validation accuracy which shows mostly
 no signs of overfitting. Therefore, even though the number of wrongly predicted participants with Alzheimer's was
 quite high, false negatives were exceptionally low, making this model a safe and significantly accurate model to use.
 However, further research should be made to minimise the number of false positives and make it even more
 accurate.
<div align="center">
 <img width="989" alt="image" src="https://github.com/maven2306/genes_Alzheimer-s/assets/169473359/efc8a588-62ab-4f8c-8c3c-4fe26c11a0c3">
</div>  
<div align="center">
  Figure 7. The following parameters were used: options = {'Verbose', true, 'Plots',' training progress,
  'shuffle', 'every-epoch', 'InitialLearnRate', 0.001, 'L2Regularization', 0.1, 'Momentum', 0.99, 'ValidationPatience', 50}.
</div>  

## DISCUSSION
The purpose of this report was to explore deeper which roles play several of the most significantly associated genes
 in Alzheimer’s disease. Thanks to the publicly available datasets, we ran analyses that allowed us to obtain valuable
 insights into the expression patterns underlying AD.
 Our GWAS analysis identified 44 genes that showed statistically significant differential expression between AD
 patients and healthy controls. Among these genes, we selected the four most significant ones based on their low
 p-values. These genes include TOMM40, APOE, PVRL2, and APOC1, all of which have been previously implicated in
 ADasfound in the literature.
 To further explore the expression patterns of these four genes, we conducted transcriptome analyses in different
 brain regions, specifically the Entorhinal Cortex, Hippocampus, Middle Temporal Gyrus, Posterior Cingulate, and
 Superior Frontal Gyrus. The analysis showed significant differential expression patterns for each gene in specific brain
 regions. TOMM40 exhibited a statistically different expression pattern in the Superior Frontal Gyrus (SFG), APOC1 in
 the Posterior Cingulate (PC), PVRL2 in the Hippocampus (HIP), and APOE in the Middle Temporal Gyrus (MTG).
 In addition, we explored the regional expression profile of APOE in the brain using data from the Allan Human Brain
 Institute. Our analysis highlighted the top ten cortical regions with the highest APOE expression, showing that APOE
 is most expressed in the hippocampus. This observation aligns with the well-established involvement of the
 hippocampus in AD pathology and proves once more the significance of APOE in AD. These findings, supported also
 by the transcriptome analysis mentioned above, suggest that the genes involved in AD pathogenesis are not equally
 dysregulated across the brain.
 Lastly, MRI analysis provided insights into the structural differences (specifically the brain volume) between AD
 patients and healthy controls. The analysis revealed that the most shrunk cortical region in patients with AD is the
 parahippocampal region. These findings complement the genetic and expression studies, providing a deeper
 understanding of the pathological changes associated with AD.
 Based on the significant differential expression patterns of TOMM40, APOE, PVRL2, and APOC1 in AD patients and the
 collective findings from the various analyses conducted in this study, we propose the exploration of targeted drug
 development. These genes represent promising targets for therapeutic interventions to mitigate the progression or
 prevent AD onset. By collaborating with pharmaceutical companies, resources can be allocated to investigate further
 the potential of developing drugs that modulate the expression of these genes in the relevant brain regions.
 In conclusion, this data report has identified the most significant genes associated with late-onset Alzheimer's
 disease and characterized their expression patterns in specific brain regions. The findings contribute to our
 understanding of the genetic and molecular factors underlying AD and provide insights for potential therapeutic
 interventions. Further research and drug development efforts in targeting TOMM40, APOE, PVRL2, and APOC1 could
 pave the way for improved treatments and management strategies for AD patients.

