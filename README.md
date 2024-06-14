<div align="center">
 <h1>MOST SIGNIFICANT GENES IN LATE-ONSET ALZHEIMER’S DISEASE</h1>
</div> 

<div align="center">
 <h2>This project was conducted collaboratively by Matteo Venturini and Ana Esteve Garcia.</h2>
</div> 


### Table of Contents

1. [GWAS ANALYSIS](#GWAS-ANALYSIS)
2. [TRANSCRIPTOME ANALYSIS](#TRANSCRIPTOME-ANALYSIS)
3. [AHB EXPRESSION ANALYSIS](#AHB-EXPRESSION-ANALYSIS)
4. [MRI ANALYSIS](#MRI-ANALYSIS)
5. [MACHINE LEARNING](#MACHINE-LEARNING)
6. [DISCUSSION](#DISCUSSION)
7. [MANAGEMENT PLAN](#MANAGEMENT-PLAN)
8. [REFERENCES](#REFERENCES)

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


##  MANAGEMENT PLAN
We used several datasets. The first one comes from studies on AD by Jansen et al. (2019), obtained from the
 database FUMA GWAS (https://fuma.ctglab.nl/). We performed a GWAS analysis on the “AD_magma.genes.txt” file.
 Then we performed a transcriptome analysis on the file “GSE5281_series_matrix.txt.gz”. This dataset (which is based
 on the studies of Liang WS et al., and Readhead B et al.) was obtained from this website:
 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281,
 and
 we found the metadata file
 “GPL570-55999_metadata.txt” on this website: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570. To
 perform our PPI analysis, we used the file “intact_int.txt’” from a study by Sügis et al.
 (https://www.nature.com/articles/s41597-019-0152-0).
 On the same website, we found the file
 “'ENGS_SYMBOL_list_intact.txt” which is the metadata file, and also the dataset “alz_intact_int.txt”. The
 https://string-db.org/ website was used to compare our PPIs analysis with their database. We then used a
 transcriptomic dataset from the Allan Human Brain Institute (http://human.brain-map.org/) to check the expression
 levels of more than 20.000 genes across 34 cortical regions in both hemispheres. The dataset used is named
 “AHBA_transcriptomics_114atlas.mat”. We also used MRI data taken from the OASIS and ADNI databases
 (https://www.oasis-brains.org/,
 http://adni.loni.usc.edu/).
 We used the data from this file
 “SURASIS_grouped_regionProperties_FS_aparc2.mat” and these files for metadata: “SURASIS_mr_sessions.csv” and
 “SURASIS_clinical_data.csv”. We decided to run a GWAS analysis, an AD gene expression analysis, a protein-protein
 interaction analysis, an AHB expression analysis, and an MRI analysis. All the analyses were run using MATLAB, Version: 24.1.0.2578822. The Curve Fitting Toolbox (version 24.1), the Deep Learning Toolbox
 (version 24.1), and the Statistics and Machine Learning Toolbox (version 24.1) were used. 

 ## REFERENCES 
 - Knopman, D. S., Amieva, H., Petersen, R. C., Chételat, G., Holtzman, D. M., Hyman, B.
 T., Nixon, R. A., & Jones, D. T. (2021). Alzheimer disease. Nature reviews. Disease
 primers, 7(1), 33. https://doi.org/10.1038/s41572-021-00269-y
 - Lee, E. G., Chen, S., Leong, L., Tulloch, J., & Yu, C. E. (2021). TOMM40 RNA
 Transcription in Alzheimer's Disease Brain and Its Implication in Mitochondrial
 Dysfunction. Genes, 12(6), 871. https://doi.org/10.3390/genes12060871
 - Liang, X., Liu, C., Liu, K., Cong, L., Wang, Y., Liu, R., Fa, W., Tian, N., Cheng, Y., Wang, N.,
 Hou, T., Du, Y., & Qiu, C. (2022). Association and interaction of TOMM40 and PVRL2
 with plasma amyloid-β and Alzheimer's disease among Chinese older adults: a
 population-based study. Neurobiology of aging, 113, 143–151.
 https://doi.org/10.1016/j.neurobiolaging.2021.12.013
- Long, J. M., & Holtzman, D. M. (2019). Alzheimer Disease: An Update on
 Pathobiology and Treatment Strategies. Cell, 179(2), 312–339.
 https://doi.org/10.1016/j.cell.2019.09.001
 - Serrano-Pozo, A., Das, S., & Hyman, B. T. (2021). APOE and Alzheimer's disease:
 advances in genetics, pathophysiology, and therapeutic approaches. The Lancet.
 Neurology, 20(1), 68–80. https://doi.org/10.1016/S1474-4422(20)30412-9
 - Sügis, E., Dauvillier, J., Leontjeva, A. et al. HENA, heterogeneous network-based data
 set for Alzheimer’s disease. Sci Data 6, 151 (2019).
 https://doi.org/10.1038/s41597-019-0152-0
 - Thompson, P. M., Hayashi, K. M., de Zubicaray, G., Janke, A. L., Rose, S. E., Semple, J.,
 Herman, D., Hong, M. S., Dittmer, S. S., Doddrell, D. M., & Toga, A. W. (2003).
 Dynamics of gray matter loss in Alzheimer's disease. The Journal of neuroscience :
 the official journal of the Society for Neuroscience, 23(3), 994–1005.
 https://doi.org/10.1523/JNEUROSCI.23-03-00994.2003
- Zhou, Q., Zhao, F., Lv, Z. P., Zheng, C. G., Zheng, W. D., Sun, L., Wang, N. N., Pang, S., de
 Andrade, F. M., Fu, M., He, X. H., Hui, J., Jiang, W., Yang, C. Y., Shi, X. H., Zhu, X. Q.,
 Pang, G. F., Yang, Y. G., Xie, H. Q., Zhang, W. D., … Yang, Z. (2014). Association
 between APOC1 polymorphism and Alzheimer's disease: a case-control study and
 meta-analysis. PloS one, 9(1), e87017.
 https://doi.org/10.1371/journal.pone.0087017
