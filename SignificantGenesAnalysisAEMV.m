%% Matteo Venturini (2688785) & Ana Esteve Garcia (2725481)

%% GWAS ANALYSIS
%week 2
% Load in the GWAS GENE data
addpath("DATA")
addpath("scripts")
tgtable=readtable("DATA/AD_magma.genes.txt");

% Analyze table so get rid of the columns we don't want

columns_to_remove={'START','STOP','NSNPS','NPARAM', 'N'};
tgtable_new=removevars(tgtable, columns_to_remove);

%make a GWAS plot
figure, hold on
pos=[1:size(tgtable_new.P)];
plot(pos,-log10(tgtable_new.P(:)),'.')
odd= find(mod(tgtable.CHR,2)==1);
plot(pos(odd),-log10(tgtable_new.P(odd)),'.r')
for ii=1:22
tmp=find(tgtable_new.CHR==ii);
XT(ii)=round((tmp(end)-tmp(1))./2) + tmp(1);
XTlabels{ii}=(['CH',num2str(ii)]);
end                 
xticks(double(XT));
xticklabels(XTlabels);
y = -log10(10.^-8);
line([1,pos(end)],[y,y],'LineStyle','--')

% P-value distribution in histogram
figure,
histogram(tgtable_new.P)
xlabel('p-value')
ylabel('frequency')
title('P-value distribution')

%location where APOE is found in the text file and its ZSTAT score
APOE_gene = find(contains(tgtable_new.SYMBOL, 'APOE'));
z_score = tgtable_new.ZSTAT(APOE_gene);

%highest Zscore with correspondant gene SYMBOL 
max_z_score = max(tgtable_new.ZSTAT);
max_index = find(tgtable_new.ZSTAT == max_z_score);
highest_symbols = tgtable_new.SYMBOL(max_index);

%make code that finds the genes that reach the significance thresholds

significance_AD_GENES = find(tgtable_new.P<10^-8);% find the index of the genes with a p-value below the threshold
num_sign_AD_GENES = length(find(tgtable_new.P<10^-8));% number of significant genes


tgtable_significant= tgtable_new(significance_AD_GENES,:); % creates a table with only the significant genes
sortedTable = sortrows(tgtable_significant, 'P');%sort them out in order of P values from most significant to the least

significant_genes_ad = tgtable_significant.SYMBOL

selectedSignificantGenes_index = find(tgtable_new.P<10^-50); % Select the four genes with the lowest p-value. The number of genes selected was arbitrary. 

SignificantGenesTable = tgtable_new(selectedSignificantGenes_index,:); % make a table of the four selected genes 

SignificantGenesNames = SignificantGenesTable.SYMBOL %shows four significant genes

%% week 3 - there is a function
% Examine whether the level of expression of the list genes we found in the
% GWAS are different between AD patients and controls in regions commonly
% affected in AD

%load meta data shown in a table
MetaDataTable = readtable('DATA/SamplesMetaData.txt')

BrainRegionNames = unique(MetaDataTable.regionDescriptions, 'rows');%6 brain regions we are going to focus on

% run a for loop to find the expression values in the 6 different brain regions of the 4 selected genes. 

% newTable = table();
% 
% for i = 1:length(SignificantGenesNames);
% 
%     for j = 1:length(BrainRegionNames);
% 
%         [p, t] = geneExpressionDifference(SignificantGenesNames{i},BrainRegionNames{j});% this calls the geneExpressionDifference function to give the p value and t test value for each brain region in the 4 significant ad genes
% 
%         Gene = [SignificantGenesNames{i}];%name of gene
%         Region = [BrainRegionNames{j}];%name of region
% 
%         T = table({Gene},{Region},p,t);
%         T.Properties.VariableNames = ["Gene","Brain Region", "P value", "T value"];%sets the variable names to each column
% 
%         newTable = [newTable;T];
%         newTable.Properties.VariableNames = ["Gene","Brain Region", "P value", "T value"];% new table made with the 4 variable names
%     end
% 
% end

% newTable %tabe showing the p value and t values of each of the 6 brain regions in the 4 most significant AD genes

%EC, HIP, MTG, PC, SFG, VCX

% Get all the info data we will use 
TRAWGeneExpdata=readtable('ExpressionData.txt');
TMETAGeneExpdata=readtable('SamplesMetaData.txt');
Tprobes=readtable('ProbeNames.txt');


%make a table with probe_id and gene_symbol
tidx1 = find(contains(Tprobes.Properties.VariableNames,'probid'));
tidx2 = find(contains(Tprobes.Properties.VariableNames,'gene_symbol'));
tbl_probeNames = Tprobes(:,[tidx1:tidx2]);
TRAWGeneExpdata_renamed = renamevars(TRAWGeneExpdata, 'probe_id', 'probid');%rename the column probe_id to probid
tbl_expression_probeNames = join(tbl_probeNames,TRAWGeneExpdata_renamed,"Keys","probid");% table showing the probid names with its correspondent gene_symbol and samples 
tbl_probeNames_clean = rmmissing(tbl_expression_probeNames);% clean table removing the ones that have a missing value, so it will not give an error

%all the other coding is in the function including the violin plots and
%histograms for each brain region in the four most significant ad genes

%% Week 4

% What we could do is examine the level of involvement of SNP and
% corresponding gene variants of genes that show the largest levels of differentiating gene expression

RegionOfInterest='HIP'; %our gene of interest
meta_REGidx = (strcmp(MetaDataTable.regionDescriptions, RegionOfInterest)); % Here we find the index of HIP in the metadata table
tbl_meta_RegionOfInterest = TMETAGeneExpdata(meta_REGidx, :); % here we display the table with the meta data (Tissue Sample, Group ID, subject ID) in the region HIP

%select the samples out of the total set of probes and samples
tmpmatch_tissueSamples_indx = find(ismember(tbl_expression_probeNames.Properties.VariableNames, tbl_meta_RegionOfInterest.tissueSampleDescriptions));
tmpSampleVals = tbl_expression_probeNames(:, tmpmatch_tissueSamples_indx);

% define which samples are from the controls and which are from the
% patients by doing their indexes
meta_CONindx = strcmp(tbl_meta_RegionOfInterest.groupID, 'control');
meta_PATindx = strcmp(tbl_meta_RegionOfInterest.groupID, 'affected');

%select the values of these samples and perform a t-test between the two
tmpSampleVals_val=tmpSampleVals{:,:}
ConVal=tmpSampleVals_val(:, meta_CONindx)
PatVal=tmpSampleVals_val(:, meta_PATindx)
[a,b,c,d]=ttest2(ConVal',PatVal')
[aa,bb] = sort(b) %sort the p-values in ascending order
tmptop10_probes = bb(1:10) %select the 10 probes with the lowest p-values
GENElist_top10 = tbl_expression_probeNames.gene_symbol(tmptop10_probes)'%list of the 10 genes that match these top10_probes of tbl_probeNames

% find the genes in the GWAS file and get their Z score

clear tmpzstats
for ii = 1:10
    a = contains(tgtable_new.SYMBOL, GENElist_top10{ii})
    tmpzstats(ii) = mean(tgtable_new.ZSTAT(a))
end

%average the values to get a mean 'GWAS ZSCORE'; we take the absolute
mean_gwas_zscore=nanmean(abs(tmpzstats))

%WEEK 4.2- PROTEOMICS-PPI interactions

%read files as tables
TPPI_intact = readtable('intact_int.txt')
TENGS_SYMBOL = readtable('ENGS_SYMBOL_list_intact.txt')

%we are going to look which protein protein interactions interact with APOE
SIGNIFICANT_AD_GENES = TENGS_SYMBOL.gene_name(significance_AD_GENES) % 44 genes found earlier, found from the TENGS_SYMBOL.gene_name
GENE_of_INTEREST = 'APOE'%gene of interest APOE
clear PPIs_of_interest; tel=0;
id=find(strcmp(TENGS_SYMBOL.gene_name, GENE_of_INTEREST));% Finds where the TENGS_SYMBOL table where the gene names match the APOE gene
id2=find(strcmp(TPPI_intact.ensg1, TENGS_SYMBOL.ensg{id})); %Finds in the TPPI_intact table where the ENSG1 matches the ENSG of the APOE gene
d = ismember(TENGS_SYMBOL.ensg, TPPI_intact.ensg2(id2)); %New array made where true corresponds to genes in TENGS_SYMBOL that match the ENSG2 in the TPPI_intact table where ENSG1 matches the APOE gene
PPIs_of_interest_A = TENGS_SYMBOL.gene_name(d) %list of genes that interact with APOE in PPI

%repeat other way around, PPI from B to A
id2=find(strcmp(TPPI_intact.ensg2, TENGS_SYMBOL.ensg{id})); %Finds in the TPPI_intact table where the ENSG2 matches the ENSG of the APOE gene
d = ismember(TENGS_SYMBOL.ensg, TPPI_intact.ensg1(id2)); %New array made where true corresponds to genes in TENGS_SYMBOL that match the ENSG1 in the TPPI_intact table where ENSG1 matches the APOE gene
PPIs_of_interest_B = TENGS_SYMBOL.gene_name(d) %list of genes that interact with APOE in PPI

PPIs_of_interest = [PPIs_of_interest_A;PPIs_of_interest_B]%total PPIs, A+B

% PPI performed in all the other perviously found significant genes by
% creating a loop

clear PPIs_of_interest; %clear any previous value given to PPIs_of_interest
tel=0; % this keeps track of the number of PPIs of interest found
for jj = 1:length(significant_genes_ad)%loop into each significant AD gene
    PROT_of_INTEREST = significant_genes_ad(jj); %SIGNIFICANT_AD_GENES  now are assigned as PROT_of_INTEREST
    id=find(strcmp(TENGS_SYMBOL.gene_name, PROT_of_INTEREST));%search the index in the TENGS_SYMBOL.gene_name that matches the PROT_of_INTEREST
    if nnz(id)==0; continue; end %checks whether a match is found, if there is non, the loop continues to the next itineration
    
    id1=find(ismember(TPPI_intact.ensg1, TENGS_SYMBOL.ensg(id)));%looks for the indices in TPPI_intact.ensg1 that match the ENSG TENGS_SYMBOL 
    d1=find(ismember(TENGS_SYMBOL.ensg, TPPI_intact.ensg2(id1)));%looks for the indices in TENGS_SYMBOL.ensg that match the ENSG2 in TPPI_intact(esng1 column)
    
    id2=find(ismember(TPPI_intact.ensg2, TENGS_SYMBOL.ensg(id)));%looks for the indices in TPPI_intact.ensg2 that match the ENSG TENGS_SYMBOL 
    d2=find(ismember(TENGS_SYMBOL.ensg, TPPI_intact.ensg1(id2)));%looks for the indices in TENGS_SYMBOL.ensg that match the ENSG2 in TPPI_intact(esng2 column)
    
    d=[d1;d2]; %combines all indices
    d=unique(d);
    
    for ii=1:size(d,1)
        tel=tel+1;
        PPIs_of_interest(tel) = TENGS_SYMBOL.gene_name(d(ii));
    end
end

PPIs_of_interest = unique(PPIs_of_interest);%list of genes that interact(PPI) with the genes in SIGNIFICANT_AD_genes

%put the genes into a table
TPPI = table(PPIs_of_interest(:));
TPPI.Properties.VariableNames = {'gene_name'}

%% 4.3 PPI-NETWORKS

TPPI_AD_intact = readtable ('alz_intact_int.txt');%load data in a table
ADproteins = unique([TPPI_AD_intact.ensg1, TPPI_AD_intact.ensg2])% Concatanate ensg1 and ensg2 in one list called AD proteins

unique_ensg1 = unique(TPPI_AD_intact.ensg1);%unique will remove any repeated elements
unique_ensg2 = unique(TPPI_AD_intact.ensg2);
unique_ensg_all = unique([unique_ensg1; unique_ensg2]);%list of unique IDs  involved in the AD PPIs

PPI_AD_matrix = zeros(size(unique_ensg_all, 1));%starting new matrix

%for each itineration it gets the values in the esng1 and ensg2 columns
%using the index in ii and assigns this values to t1 and t2
for ii = 1:size(TPPI_AD_intact, 1)
    t1 = TPPI_AD_intact.ensg1{ii};
    t2 = TPPI_AD_intact.ensg2{ii};
    PPI_AD_matrix(contains(unique_ensg_all,t1),contains(unique_ensg_all,t2))=TPPI_AD_intact.score(ii);%makes matrix
end
 
%making matrix bidirectional
PPI_AD_matrix_bidirectional = PPI_AD_matrix + PPI_AD_matrix.'

%make the bidirectional matrix, binary
%PPI_AD_matrix_bidirectional_binary = PPI_AD_matrix_bidirectional>0
PPI_AD_matrix_bidirectional_binary = logical(PPI_AD_matrix_bidirectional);


% plot graph- to see if the Proteins are interacting
G = graph(double(PPI_AD_matrix_bidirectional_binary))
G.Nodes.Names = unique_ensg_all;
figure;
plot(G, 'NodeLabel', G.Nodes.Names, 'LineWidth', G.Edges.Weight, 'Layout','force')

%Find the most interactive protein
PPI_network = sum(PPI_AD_matrix_bidirectional_binary);
[max_degree, max_index] = max(PPI_network);%max degree and index where the most interactive protein is located
most_interactive_protein = G.Nodes.Names(max_index);
%make histogram with all the PII values
allPPIvalues = nonzeros(PPI_AD_matrix);
figure, histogram(allPPIvalues);

%graph with the strengths of the PPIs
PPI_AD_matrix_bidirectional_weighted = PPI_AD_matrix_bidirectional*10;
G_weighted = graph(PPI_AD_matrix_bidirectional_weighted);
G_weighted.Nodes.Names = unique_ensg_all;
figure;
plot(G_weighted, 'NodeLabel', G_weighted.Nodes.Names, 'LineWidth',G_weighted.Edges.Weight, 'Layout', 'force');

%Are there any subclusters of proteins within our network
PPI_AD_matrix_bidirectional = PPI_AD_matrix | PPI_AD_matrix.';
[a, b] = modularity_und(PPI_AD_matrix_bidirectional);%a and b represent the module assignments and the modularity value
numModules = numel(unique(a));% count of distinct modules in the network
nodes_grouped_together = splitapply(@(x) {x}, G.Nodes.Names, a);%each cell contains the names of the nodes belonging to the same module

%We found there are 15 clusters of different proteins. Here we are going to
%list them

% Count the number of proteins in each module group
moduleCounts = histcounts(a, 1:numModules+1);

% Sort the 15 groups in descending order based on how many protein counts they have
[sortedCounts, sortedIndices] = sort(moduleCounts, 'descend');

% Select the 15 modules with the highest protein counts
ModulesGroups = sortedIndices(1:15);

% Create a cell array to store the groups of proteins
proteinGroups = cell(15, 1);

% Iterate over the top modules and store the proteins in each group
for i = 1:15
    moduleIndex = ModulesGroups(i);
    proteins = nodes_grouped_together{moduleIndex};
    proteinGroups{i} = proteins;
end

% Display the list of protein groups
for i = 1:15
    fprintf('Group %d: %s\n', i, strjoin(proteinGroups{i}, ', '));
end





%% Week 5 
%5.1.1-5.1.4
AHBA = load('AHBA_transcriptomics_114atlas.mat')%load AHBA data

% We look for APOE in AHBA
APOE_AHBA = find(ismember(AHBA.gene_symbol, "APOE"));

% We look for the left hemisphere regions index 
lh_region = find(contains(AHBA.regionDescription, 'lh'));

% We look for the names of the left hemisphere regions
lh_region_genes = AHBA.regionDescription(lh_region);

% We look for the expression levels of the brain regions in left hemisphere
AHBA_lefthemi = AHBA.gene_expression_region_FS_z(lh_region, :, :);


% Now we calculate the average of the 3rd dimension (we average out the
% values of the 6 patients, so we go from a 3d matrix to a 2d matrix)

AHBA_lefthemi_group_average = nanmean(AHBA_lefthemi, 3);

% Q5.1.5) 
% Now extract the transcriptome profile of APOE across the 57 regions:

APOE_expression_DK114_lh = AHBA_lefthemi_group_average(:, APOE_AHBA);


% Q 5.1.6
% Create a histogram with 20 bins
histogram(APOE_expression_DK114_lh, 20)

% Add x- and y-labels
xlabel('Data values')
ylabel('Frequency')

% Q 5.1.7
greenColormap = [zeros(1001,1),linspace(0,1,1001)',zeros(1001,1)];%to make the colormap of the brain green(the lighter the green the more is expressed)
colormap=greenColormap;
regionNames = AHBA.regionDescription; %list of names of the brain regions

% Now we create a cortical map
corticalMap_LR = [APOE_expression_DK114_lh; zeros(57,1)]'; 

plotSurfaceBoth(regionNames, corticalMap_LR, colormap);

% Q 5.1.9 

[sorted_expression, sorted_indices] = sort(APOE_expression_DK114_lh, 'descend');%sort the expressions in descending order
top_10_regions = regionNames(sorted_indices(1:10))%top 10 most expressed regions from the APOE_expression_DK114_lh

% Q 5.1.10

% Now we make a table with the 10 most expressed regions and their
% expression levels
top_10_expressions = sorted_expression(1:10);

rank = [1:10];%rank them 1 to 10

top_10_regions_table = table(rank', top_10_regions, top_10_expressions, 'VariableNames', {'RANK','REGION_NAME','EXPRESSION_LEVEL'})%put them into a table


% Now we make a table with the significant AD genes and the brain region of
% the left hemisphere where they are most expressed

topExpressionTable = table();

% This is to find the names of the significant ad genes 
AD_genes = tgtable_new.SYMBOL(significance_AD_GENES);

% We run a while loop with the code used before 
kk = 1;
while kk < length(AD_genes);

    kk;   

    current_gene = AD_genes{kk};% current gene is each of the 44 AD genes

    current_gene_str = string(current_gene);%Cconvert the gene to a string

    
    APOE_AHBA = find(ismember(AHBA.gene_symbol, current_gene_str)); % Find the index of current_gene_str in AHBA.gene_symbol
    
    if ~isempty(APOE_AHBA)
        GENE_expression_DK114_lh = AHBA_lefthemi_group_average(:, APOE_AHBA); %Get the gene expression values for current_gene_str
        [sorted_expression, sorted_indices] = sort(GENE_expression_DK114_lh, 'descend'); %sort them in descending order
        top_1_region = regionNames(sorted_indices(1));% top region with highest expression

        expression_table = table(cellstr(current_gene_str),top_1_region); %Make an expression table with the gene and the most expressed region
        expression_table.Properties.VariableNames = ["Gene","Most Expressed Region"];%name variable names

        topExpressionTable = [topExpressionTable;expression_table]; %add the expression table to topExpressionTable
        topExpressionTable.Properties.VariableNames = ["Gene","Most Expressed Region"]; %add varaible names
    else
        disp(['Gene not found: ' current_gene_str])%display 'gene not found' if the gene is not found in the AHBA.gene_symbol
    end

    kk = kk+1;%move to next gene index

end

topExpressionTable


% 5.2

% It's the violin plot- added in week 3 


% 5.3 
%load data
IPS_DATA = readtable('GSE102956_COUNTS_IPS.txt');
MIC_DATA = readtable('GSE102956_COUNTS_MIC.txt');
NEU_DATA = readtable('GSE102956_COUNTS_NEU.txt');
AST_DATA = readtable('GSE102956_COUNTS_AST.txt');


%Q5.3.4 Select the columns of APOE3 and APOE4
TYPEofData='MIC'
filename=['GSE102956_COUNTS_',TYPEofData,'.txt'];
iPSC_DATA=readtable(filename);

%read in the variable names separately
fid=fopen(filename);
varNames = strsplit(fgetl(fid), '\t');
fclose(fid);

%add a unique identifier to each varName, otherwise you cant use it as columnnames 
for ii=1:length(varNames)
    varNames_updated{ii} = [varNames{ii},'_',num2str(ii)];
end
%alternative:%You can do this by means of a for loop yourself, or google and find this matlab inbuild
varNames_updated=matlab.lang.makeUniqueStrings(varNames); 
varNames=varNames_updated; %updated the varNames to the varNames_updated

%use the varNames_updated to update the varNames of table
iPSC_DATA.Properties.VariableNames=["GENE", varNames];        
%To continue we can then:
%Q5.3.5 - %select which columns are APOE3 and which are APOE4 neurons
APOE3neurons=strfind(iPSC_DATA.Properties.VariableNames,[TYPEofData,'_E3']);
APOE4neurons=strfind(iPSC_DATA.Properties.VariableNames,[TYPEofData,'_E4']);
%select which columns are APOE3 and which are APOE4 neurons - a for loop solution
tmp=[];
for ii=1:size(APOE3neurons,2)
  if(APOE3neurons{ii} == 1)
      tmp(ii)=1;
  else
      tmp(ii)=0;
  end
end
APOE3columnsInTable=find(tmp==1);
tmp=[];
for ii=1:size(APOE4neurons,2)
  if(APOE4neurons{ii} == 1)
      tmp(ii)=1;
  else
      tmp(ii)=0;
  end
end
APOE4columnsInTable=find(tmp==1);

%Alternatively to reduce our code to just 2 lines
APOE3columnsInTable=find(~cellfun(@isempty, APOE3neurons));
APOE4columnsInTable=find(~cellfun(@isempty, APOE4neurons));

%select the data we want
%Q5.3.6 - select the data matching APOE3 
Data3=iPSC_DATA{:,APOE3columnsInTable};
Data4=iPSC_DATA{:,APOE4columnsInTable};

%remove all the empty cells, and set them to NaN so they don't affect the
Data3(Data3 == 0) = NaN;
Data4(Data4 == 0) = NaN;
%ttest
[H, P, CI, STATS] = ttest2(Data3', Data4');

%this is done to see how many tests we can perform
NumberOfTests = size(Data3, 1);%Data 3 and Data4 have same number of rows so doesn't matter which one you use
alpha = 0.05 / NumberOfTests;%alpha level 0.05

iPSCdiffExpGenes = iPSC_DATA.GENE(P<alpha);
iPSCdiffExpGenes_SIGN = sign(STATS.tstat(P<alpha));

%% Week 6 

MRIdata=load('SURASIS_grouped_regionProperties_FS_aparc2.mat');

%load in the meta data
tbl_meta_data_mri=readtable('SURASIS_mr_sessions.csv');
tbl_meta_data_clinical=readtable('SURASIS_clinical_data.csv');

subjectsMRI = MRIdata.subjects'

% Here we find the index of the patients that are in MRIdata but also in tbl_meta_data_mri 
findSubjects_metaMRI = find(ismember(tbl_meta_data_mri.MRID, subjectsMRI))

% Here we make a new table of tbl_meta_data_mri that contains only the
% subjects also in MRIdata
tbl_meta_data_mri_new = tbl_meta_data_mri(findSubjects_metaMRI, :);

%gets only the subject column
subjectsClinical = tbl_meta_data_mri_new.Subject

%index of all clinical subjects which are present (multiple) times in the
%tbl_meta_data_clinical table
findSubjects_metaClinical = find(ismember(tbl_meta_data_clinical.Subject, subjectsClinical))

%show all columns
tbl_meta_data_clinical_new = tbl_meta_data_clinical(findSubjects_metaClinical, :)

%find for all subjects in MRIdata their PATIENT CONTROL STATUS and add this
%to the meta_data_mri_sessions table
for s = 1: numel(MRIdata.subjects);
    tmpsub_mri_session = MRIdata.subjects{s};

    % step 1: find the index of tmpsub_mri_session in tbl_meta_data_mri_filered. 
    
    tmp1=find(strcmp(tbl_meta_data_mri_new.MRID, tmpsub_mri_session));
    tmpsubcode = tbl_meta_data_mri_new.Subject{tmp1};

    % step 3: find the clinical status of this subject in tbl_meta_data_clinical

    tmp2=find(strcmp(tbl_meta_data_clinical.Subject, tmpsubcode));


    % step 4: extract the status
    tmpstatus=tbl_meta_data_clinical.dx1(tmp2);
    
    % Each patient has multiple timepoint measurments, so we need to
    % differentiate into each one of them
    if strcmp(tmpstatus, 'Cognitively normal')
    tbl_meta_data_mri_new.status{tmp1} ='CON';
    elseif strcmp(tmpstatus, 'AD Dementia')
    tbl_meta_data_mri_new.status{tmp1} = 'AD';
    else
    tbl_meta_data_mri_new.status{tmp1} = 'Other';
    end
end


% find the region middle temporal gyrus
regions = find(contains(MRIdata.regionDescriptions,'middletemporal'));
%gets the total volume of the brain region
dataProp = find(contains(MRIdata.propertyDescriptions,'GrayVol'));
%one region for the left hemisphere and one for the right
ROIvol = squeeze(MRIdata.regionProperties(regions,dataProp,:));

ROIvol_mean = mean(ROIvol);

%remove outliers, which are 3sd from the mean

outliers = ROIvol_mean < nanmean(ROIvol_mean)-3*nanstd(ROIvol_mean) | ROIvol_mean > nanmean(ROIvol_mean)+3*nanstd(ROIvol_mean);

ROIvol_mean(outliers) = NaN; %we exclude these values from further analysis



% determine who are Controls and who are the Patients
CONindx = find(strcmp(tbl_meta_data_mri_new.status,'CON'));
PATindx = find(strcmp(tbl_meta_data_mri_new.status,'AD'));
% Plot the data - use your Violin plot for this
plotdata.Controls = ROIvol_mean(CONindx);
plotdata.Patients = ROIvol_mean(PATindx);
figure; violinplot(plotdata);
ylabel('Volume (mm3)');
title('Regional gray matter volume');

[h,p,~,stats]=ttest2(ROIvol_mean(CONindx),ROIvol_mean(PATindx));

disp (p)

%difference in volume?

ROIvol_everyregion=squeeze(MRIdata.regionProperties(:,dataProp,:)); % Extract region volumes from MRIdata for each region
outliers_everyregion =  ROIvol_everyregion < nanmean(ROIvol_everyregion)-3*nanstd(ROIvol_everyregion) | ROIvol_everyregion > nanmean(ROIvol_everyregion)+3*nanstd(ROIvol_everyregion); % Identify outliers based on 3 standard deviations from the mean

CONindx_everyregion = find(strcmp(tbl_meta_data_mri_new.status,'CON')); % Find indices of control subjects
PATindx_everyregion = find(strcmp(tbl_meta_data_mri_new.status,'AD')); % Find indices of AD patients

CONindx_everyregion(find(ismember(CONindx, outliers))) = [];%remove outliers 
PATindx_everyregion(find(ismember(PATindx, outliers))) = [];%remove outliers
[a,b,c,TSTATS]=ttest2(ROIvol_everyregion(:,CONindx_everyregion)',ROIvol_everyregion(:,PATindx_everyregion)'); % Perform a two-sample t-test between control and AD groups using region volumes
AD_pathology_score_DK68 = TSTATS.tstat;%tTEST as AD pathology scores for each region

%% top 10 regions

% Sort regions based on effect size
[sorted_effect, sorted_idx] = sort(AD_pathology_score_DK68, 'descend');
top_10_regions = MRIdata.regionDescriptions(sorted_idx(1:10));

% Create table
top_10_regions_table = table(top_10_regions, sorted_effect(1:10)', 'VariableNames', {'Brain Region', 'Effect Size'});

% Display the table
disp(top_10_regions_table);

% Create a colormap for the plot
colormap('cool');

% Plot the effects
bar(AD_pathology_score_DK68);
xlabel('Brain Region');
ylabel('Effect Size');
title('Distribution of Effects across Cortical Regions');

%colormao of the cortical map of the brain shown in pink where is most
%shown the difference
cm = [linspace(1, 0, 1001); zeros(1, 1001); linspace(1, 0.6, 1001)]'
cm = flipud(cm);

plotSurfaceBoth_DK68(MRIdata.regionDescriptions, AD_pathology_score_DK68, cm);


% 68 regions we have in the cortical atlas we have to get them from the 114
% Cammoun DK atlas, which are subregions of the 68 cortical regions
for ii=1:size(MRIdata.regionDescriptions)/2 
DK68region = MRIdata.regionDescriptions{ii}% Get back the DK68 region at index ii from MRIdata
AHBAregions = find(contains(AHBA.regionDescription, DK68region)) % Find the indices of AHBA regions that contain DK68region
APOE_expression_DK68_lh(ii) = nanmean(APOE_expression_DK114_lh(AHBAregions))% Calculate the mean APOE expression for the AHBA regions related to DK68region
end

lhregions_idx = find(contains(MRIdata.regionDescriptions, 'lh-'));%index of left hemisphere regions
rhregions_idx = find(contains(MRIdata.regionDescriptions, 'rh-'));%index of right hemisphere regions

AD_pathology_score_DK68_lhrh = nanmean([AD_pathology_score_DK68(lhregions_idx); AD_pathology_score_DK68(rhregions_idx)]);%index of DK68 left region and right ones as well

% we are plotting the expression of APOE in x axis, in the cortical regions
% of the left hemispheres while on the y axis we plot the t test scores of
% the difference of the brain volumes between patients and controls

APOE_expression_DK68_lh=unique(APOE_expression_DK68_lh)
figure,
plot(APOE_expression_DK68_lh', AD_pathology_score_DK68_lhrh, '.')


%% 6.2 DATA VISUALISATION
f = figure
hold on
x = APOE_expression_DK68_lh';
y = AD_pathology_score_DK68_lhrh';
s = scatter(x,y)
%let's place a line through the points
l = lsline;
f.Color = 'cyan';

%time for colours
color1 = [255 0 0] ./ 255;
s.MarkerFaceColor = color1;
s.MarkerEdgeColor = color1;
s.MarkerEdgeAlpha = .3;
s.MarkerFaceAlpha = .3;

l.Color = 'black';
l.LineWidth = 1;

%labelling axis
xlabel('APOE expression pattern');
ylabel('effect size on volume');
%Thickness of axis
a=gca;
a.XAxis.LineWidth = 1;
a.YAxis.LineWidth = 1;
%font size of axis labels
a.FontSize = 12;

% ADD CONFIDENT INTERVALS
[~,I] = sort (x);
xsorted = x(I);
ysorted = y(I);

fitresult = fit(xsorted, ysorted, 'poly1')
p95 = predint(fitresult, xsorted, 0.95, 'observation', 'off');
%p90 = predint(fitresult, xsorted, 0.90, 'observation', 'off');

plot(xsorted,p95,'r--')
%plot(ysorted,p90,'r--')

%patch([xsorted; flipud(xsorted)], [p90(:,1); flipud(p90(:,2))], 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.05);
patch([xsorted; flipud(xsorted)], [p95(:,1); flipud(p95(:,2))], 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.1);
% replot the points to bring them on top
plot(xsorted,ysorted,'.b')

%% 6.3 statistics
% regression between APOE expression in lh and the t test scores of the lh and rh
% x axis is apoe because is the independent variable
s = regstats(AD_pathology_score_DK68_lhrh,APOE_expression_DK68_lh)
%we get two p values but which one is best? y on predictions of x, ttest on
%predictions of apoe
pval = s.tstat.pval(1);

disp(pval)

%6.4 CONNECTOMICS: CONNECTIVITY IN AD BRAIN
%Examining whether there is a difference in connectivity in APOE-4
%non-carriers(APOE 3) and APOE 4 carriers in the hippocampus, precuneus and
%whether there is a difference in aging effect

% define the number of APOE3 and APOE4 carriers.
nAPOE3_carriers=30;
nAPOE4_carriers=25;
% load in the META data on AGE
tmp1 = readtable('ages_list.txt','ReadVariableNames',false);
tmp1.Properties.VariableNames{1}='Age';
% load in the META data on GENDER
tmp2 = readtable('genders_list.txt','ReadVariableNames',false);
tmp2.Properties.VariableNames{1}='Gender';
% make a binary variable for GENDER, FEMALE=1, MALE=0
tmp3=tmp2;
gen=double(strcmp(tmp2.Gender,'Female'));
tmp3.Gender=gen;
tmp3.Properties.VariableNames{1}='GenderBin';
%put all demographics in one table
tblDemoGraphics = [tmp1,tmp2, tmp3];
clear tmp*
%split the Demographics table into 2 tables 1 for APOE3 and 1 for APOE4. We can also make 1 table with patient status in it as an alternative
tblDemoGraphics_APOE3 = tblDemoGraphics(1:nAPOE3_carriers,:);
tblDemoGraphics_APOE4 = tblDemoGraphics(nAPOE3_carriers+1:end,:);
%load in the META data on the RegionNames
regionNames = readtable('APOE-4_9_region_names_full_file.txt','ReadVariableNames',false);
regionNames_cleaned=[];
for ii=1:size(regionNames,1)
 tmp=regionNames{ii,:};
 totaltmp=[];
  for jj=1:size(tmp,2)
   if isempty(tmp{jj});continue;end
    totaltmp=[totaltmp,tmp{jj},'_'];
  end
%we can remove all extra ' just because we want it to look pretty
%totaltmp(strfind(totaltmp,''''))=[];
%remove the last
 totaltmp(end)=[];
%our regionNames
 regionNames_cleaned{ii}=totaltmp;
end
%put regionNames in a table and clean all intermediate steps
regionNames=regionNames_cleaned';
tbl_regionNames=cell2table(regionNames);
clear regionNames_cleaned *tmp*
%load in the connectivity matrices for the APOE3 and APOE4 group
for ii=1:nAPOE3_carriers
tmpconnectivity_matrix=load(['APOE-3_',num2str(ii),'_connectivity_matrix_file.txt']);
APOE3_matrices(ii,:,:)=tmpconnectivity_matrix;
end
for ii=1:nAPOE4_carriers
tmpconnectivity_matrix=load(['APOE-4_', num2str(ii), '_connectivity_matrix_file.txt']); %hint almost the same as in the APOE3 carrier loop
APOE4_matrices(ii,:,:)=tmpconnectivity_matrix;
end

%6.4.3 Define our region of interest
RegionOfInterest='hippocampus'
idx_roi = (strfind(tbl_regionNames.regionNames, RegionOfInterest))%indeces where in the tbl_regionNames the hippocampus is located, first left then right hemisphere
idx_roi =find(~cellfun(@isempty, idx_roi))

%6.4.5 Get the APOE4 and APOE5 connectivity values of our RegionOfInterest,
%the hippocampus in this case
%get the connectivity of all APOE3_carriers of the left hippocampus
CONN_ROI_APOE4 =nansum(squeeze(nansum(squeeze(APOE4_matrices(:,idx_roi,:)),2)),2);
CONN_ROI_APOE3 =nansum(squeeze(nansum(squeeze(APOE3_matrices(:,idx_roi,:)),2)),2);

%effect of gender
s=regstats(CONN_ROI_APOE4, tblDemoGraphics_APOE4.GenderBin);
CONN_ROI_APOE4 = s.r;
s=regstats(CONN_ROI_APOE3, tblDemoGraphics_APOE3.GenderBin);
CONN_ROI_APOE3 = s.r;

%6.4.7 Plot the effect of connectivity atrophy of the Hippocampus across
%age betwwen the two groups 

figure,plot(tblDemoGraphics_APOE3.Age,CONN_ROI_APOE3,'.'); lsline
hold on
plot(tblDemoGraphics_APOE4.Age,CONN_ROI_APOE4,'.r'); lsline

%6.4.3 Define our region of interest
RegionOfInterest='precuneous'
idx_roi = (strfind(tbl_regionNames.regionNames, RegionOfInterest))%indeces where in the tbl_regionNames the hippocampus is located, first left then right hemisphere
idx_roi =find(~cellfun(@isempty, idx_roi))

%6.4.5 Get the APOE4 and APOE5 connectivity values of our RegionOfInterest,
%the hippocampus in this case
%get the connectivity of all APOE3_carriers of the left hippocampus
CONN_ROI_APOE4 =nansum(squeeze(nansum(squeeze(APOE4_matrices(:,idx_roi,:)),2)),2);
CONN_ROI_APOE3 =nansum(squeeze(nansum(squeeze(APOE3_matrices(:,idx_roi,:)),2)),2);

%effect of gender
s=regstats(CONN_ROI_APOE4, tblDemoGraphics_APOE4.GenderBin);
CONN_ROI_APOE4 = s.r;
s=regstats(CONN_ROI_APOE3, tblDemoGraphics_APOE3.GenderBin);
CONN_ROI_APOE3 = s.r;

%6.4.7 Plot the effect of connectivity atrophy of the Hippocampus across
%age betwwen the two groups 

figure,plot(tblDemoGraphics_APOE3.Age,CONN_ROI_APOE3,'.'); lsline
hold on
plot(tblDemoGraphics_APOE4.Age,CONN_ROI_APOE4,'.r'); lsline

%%  WEEK 7- MACHINE LEARNING

%7.1-UNSUPERVISED LEARNING

MRIdata = load('SURASIS_grouped_regionProperties_FS_aparc2.mat')

%Here we select two regions of interest that we suspect to describe a good
%measurement for whether a participant is a CONTROL or a PATIENT
Region1 = 'ctx-lh-medialorbitofrontal'
Region2 = 'ctx-lh-middletemporal'
idxRegion1=find(strcmp(MRIdata.regionDescriptions, Region1))
idxRegion2=find(strcmp(MRIdata.regionDescriptions, Region2))
% REGIONAL VOLUME OF THIS REGIONS of all the patients
Reg1val=squeeze(MRIdata.regionProperties(idxRegion1, dataProp, :));
Reg2val=squeeze(MRIdata.regionProperties(idxRegion2, dataProp, :));
%normalizing data
Reg1val= (Reg1val - mean(Reg1val)) ./ std(Reg1val);
Reg2val= (Reg2val - mean(Reg2val)) ./ std(Reg2val);


%outliers
outliers_reg1 =  Reg1val < nanmean(Reg1val)-3*nanstd(Reg1val) | Reg1val > nanmean(Reg1val)+3*nanstd(Reg1val);
outliers_reg2 =  Reg2val < nanmean(Reg2val)-3*nanstd(Reg2val) | Reg2val > nanmean(Reg2val)+3*nanstd(Reg2val);

%plot
figure, plot(Reg1val, Reg2val, '.')
X=[Reg1val, Reg2val]
[indxcluster]=kmeans(X,2);%plot it in 2 groups
figure, plot(Reg1val, Reg2val, '.') %plot all together
hold on
plot(Reg1val((indxcluster==2)),Reg2val((indxcluster==2)),'.r') %ctx-lh-medialorbitalfrontal
plot(Reg1val((indxcluster==1)),Reg2val((indxcluster==1)),'.b') %ctx-lh-middletemporal 

%% 7.2 supervised learning
%We are going to train the performance of the model in the training set and
%its test performance in the test set

%MRIdata = load('SURASIS_grouped_regionProperties_FS_aparc2.mat');


%seprate them from patients and controls to train and test-set


% Calculate set sizes
testsetsize_CON = round((length(CONindx)/100)*40);%40 test controls
trainsetsize_CON = round((length(CONindx)/100)*60); %60 train controls

testsetsize_PAT = round((length(PATindx)/100)*40); %40 test patients
trainsetsize_PAT = round((length(PATindx)/100)*60); %60 train patients

% Check bounds for CONindx array
if testsetsize_CON + trainsetsize_CON > length(CONindx)
    error('Set sizes exceed the number of elements in CONindx array.');
end

% Check bounds for PATindx array
if testsetsize_PAT + trainsetsize_PAT > length(PATindx)
    error('Set sizes exceed the number of elements in PATindx array.');
end

% Randomly permute indices
tmp = randperm(length(CONindx));
trainsetCON = CONindx(tmp(1:trainsetsize_CON));%trainset control index
testsetCON = CONindx(tmp(trainsetsize_CON+1:testsetsize_CON+trainsetsize_CON));%testset control index

tmp = randperm(length(PATindx));
trainsetPAT = PATindx(tmp(1:trainsetsize_PAT));%trainset patient index
testsetPAT = PATindx(tmp(trainsetsize_PAT+1:testsetsize_PAT+trainsetsize_PAT));%testset patient index

%define input data
ROIvol = squeeze(MRIdata.regionProperties(:,dataProp,:));
ROIvol = rescale(ROIvol);

%Prepare output so if they are PATIENTS or CONTROLS and get the input of them

% make an output vector Y telling who are patients and who are controls
Ytrain = zeros(length([trainsetPAT; trainsetCON]),1);
Ytrain(1:length(trainsetPAT))=1;
% make it a boolean (the fitclinear function wants it that way)
Ytrain=(Ytrain==1);
% set the input data X
Xtrain=[ROIvol(:,trainsetPAT), ROIvol(:,trainsetCON)];
Xtrain=Xtrain';
% do the same for the test set, call this Ytest and Xtest
Ytest = zeros(length([testsetPAT; testsetCON]),1);
Ytest(1:length(testsetPAT))=1;
Ytest=(Ytest==1);
Xtest=[ROIvol(:,testsetPAT), ROIvol(:,testsetCON)];
Xtest=Xtest';

%classify data
%MDL represents the trained SVM model
%fitinfo gives information about the fitting process of the SVM model
[MDL, fitinfo] = fitclinear(Xtrain, Ytrain, 'Regularization', 'ridge')

%you get the predicted labels for the training data, 1=correct, 0=mistake
predictedLabelsTest = predict(MDL, Xtrain)
falsePositives = nnz(predictedLabelsTest(Ytest == 0) == 1);
falseNegatives = nnz(predictedLabelsTest(Ytest == 1) == 0);
%proportion of mistakes made - 2ways to do it
inssampleerror = 1-nnz((predict(MDL,Xtrain) == Ytrain))./length(Xtrain)

inssampleerror = loss(MDL,Xtrain,Ytrain);


%2nd prediction model 
% divide data into development/test/validation sets
trainRatio = 0.60;
valRatio = 0.20;
testRatio = 0.20

% Create train/test/val set.
% Note we set the default seed to make sure we get the same train/test/val subjects every run.
fprintf('random seed set to default.\n');
rng('default')

% divide the subjects of the dataset into training/validation/testsets
nsubjects = height(tbl_meta_data_mri_new);
subjects_study.train = false(nsubjects, 1);
subjects_study.val = false(nsubjects, 1);
subjects_study.test = false(nsubjects, 1);
[xtrain, xval, xtest] = dividerand(1:nsubjects, trainRatio, valRatio, testRatio);
subjects_study.train(xtrain) = true;
subjects_study.val(xval) = true;
subjects_study.test(xtest) = true;

dataProp = find(contains(MRIdata.propertyDescriptions,'GrayVol'));
ROIvol = MRIdata.regionProperties(:,dataProp,:);
% ROIvol = Array ( regions X properties X subjects )
train_data = ROIvol;
% We normalize the data
train_data = (train_data - mean(train_data, 3)) ./std(train_data, [], 3);
n_regions = size(train_data, 1);
n_properties = size(train_data, 2);

% define network architecture: we are defining how our network looks like, with our own parameters,
% all brain region has input neurones, and we define 10 in a hidden layer
% with hidden nodes . 1st node 68 , 2nd 68...
solvername = 'sgdm';
layers = [
imageInputLayer([n_regions n_properties 1], 'name','inputlayer')
fullyConnectedLayer(10, 'name', 'first layer')
reluLayer('name', 'relu1');%evaluation function
fullyConnectedLayer(2, 'name', 'second layer');
softmaxLayer('name', 'softmax');
classificationLayer('name', 'classification')];
% Set the options for the training
options = {'Verbose', true, ...
'Plots','training-progress', ...
'shuffle', 'every-epoch'};

options_study = trainingOptions(solvername,options{:},'ValidationData', ...
{permute(train_data(:,:,subjects_study.val), [1 2 4 3]), ...
categorical(strcmp(tbl_meta_data_mri_new.status(subjects_study.val), 'AD'))});
[trained_net, training_results] = trainNetwork( ...
permute(train_data(:,:,subjects_study.train), [1 2 4 3]), ...
categorical(strcmp(tbl_meta_data_mri_new.status(subjects_study.train), 'AD')),layers, options_study);

%xaxis how many times we psh the data through the network
% y is the accuracy of the model
%blue is the accuracy of the training data-starts gessing and starts to
%increase because model is learning until is plateau, the model stops
%learning

%dotted line is the performance on the validation data set
%the red is the loss, dropping in the validation data set but not as much as in the training

%evaluation on out trained model on the TEST dataset
XTest = permute(train_data(:,:,subjects_study.test), [1 2 4 3]);
YTest = categorical(strcmp(tbl_meta_data_mri_new.status(subjects_study.test), 'AD'));
[YPredicted,probs] = classify(trained_net, XTest);
% testError = 1 - accuracy
testError = 1 - mean(YPredicted == YTest)%0.0769
accuracy = 1 - testError %0.9231
% we think is quite accurate this result, which means that the train model
% actually works

confusionmat(YTest, YPredicted)

% use the function make_fair which balances the number of controls and patients during the training phase.
% Ensure same number of patients and controls in the val and train set (NOT necessary in test set).
fair_val = make_fair(strcmp(tbl_meta_data_mri_new.status,'AD'), subjects_study.val);
fair_train = make_fair(strcmp(tbl_meta_data_mri_new.status,'AD'), subjects_study.train);
% Make 3 lists of CON and CASES (i.e. Patients) out of this
CON_CASE_list_fairval = strcmp(tbl_meta_data_mri_new.status(fair_val), 'AD');
CON_CASE_list_fairtrain = strcmp(tbl_meta_data_mri_new.status(fair_train), 'AD');
CON_CASE_list_fairtest = strcmp(tbl_meta_data_mri_new.status(subjects_study.test),'AD');

% You can try out some parameter settings by including or excludingthem:
% (help trainingOptions)
options = {'Verbose', true, ...
'Plots','training-progress', ...
'shuffle', 'every-epoch', ...
'InitialLearnRate', 0.001, ...
'L2Regularization', 0.1, ...
'Momentum', 0.99, ...
'ValidationPatience', 50};
options_study = trainingOptions(solvername,...
options{:}, ...
'ValidationData', ...
{permute(train_data(:,:,fair_val), [1 2 4 3]), ...
categorical(CON_CASE_list_fairval)});
[trained_net, training_results] = trainNetwork( ...
permute(train_data(:,:,fair_train), [1 2 4 3]), ...
categorical(CON_CASE_list_fairtrain), ...
layers, options_study);
% Evaluate again!
XTest = permute(train_data(:,:,subjects_study.test), [1 2 4 3]);
YTest = categorical(CON_CASE_list_fairtest);
[YPredicted,probs] = classify(trained_net, XTest);
testError = 1 - mean(YPredicted == YTest);
accuracy = 1 - testError

%confusion matrix done to see the number of how many false labeled patients
%and controls are false negatives or false positives
confusionmat (YTest, YPredicted)







