function [p_value, t_value] = geneExpressionDifference(GeneName, BrainRegion)

% WEEK 3 % 

%EC, HIP, MTG, PC, SFG, VCX

% Step by step plan on how to get all the info
TRAWGeneExpdata=readtable('ExpressionData.txt');
TMETAGeneExpdata=readtable('SamplesMetaData.txt');
Tprobes=readtable('ProbeNames.txt');


%make a table with probe_id and gene_symbol
tidx1 = find(contains(Tprobes.Properties.VariableNames,'probid'));
tidx2 = find(contains(Tprobes.Properties.VariableNames,'gene_symbol'));
tbl_probeNames = Tprobes(:,[tidx1:tidx2]);
TRAWGeneExpdata_renamed = renamevars(TRAWGeneExpdata, 'probe_id', 'probid');
tbl_expression_probeNames = join(tbl_probeNames,TRAWGeneExpdata_renamed,"Keys","probid");
tbl_probeNames_clean = rmmissing(tbl_expression_probeNames);


% We select a REGION of interest
RegionOfInterest = BrainRegion;
meta_REGidx = (strcmp(TMETAGeneExpdata.regionDescriptions, RegionOfInterest));
tbl_meta_RegionOfInterest = TMETAGeneExpdata(meta_REGidx, :);

% We select a GENE of interest
GeneOfInterest = GeneName;
data_GENindx = strcmp(tbl_probeNames_clean.gene_symbol, GeneOfInterest); % here we look for APOE 

tmpmatch_tissueSamples_indx = find(ismember(tbl_probeNames_clean.Properties.VariableNames, tbl_meta_RegionOfInterest.tissueSampleDescriptions)); % here we look for the index of APOE in the tissue samples 
tmpSampleVals = tbl_probeNames_clean(data_GENindx, tmpmatch_tissueSamples_indx); % here we display the expression levels of APOE in each tissue sample 

%mean of each HIP region sample
%tmpSampleVals = mean(tmpSampleVals{:,:}) % here we make the mean of the 3 expression levels of APOE, given that there are 3 APOE genes

% if there's a matrix, MATLAB calculates the mean of the rows (column 1,
% mean, column 2, mean etc). But if there is only one row, MATLAB takes the
% mean of the columns and not the mean of a single row. So we need to tell
% it to take the mean of the rows, even if there's one gene. 

tmpSampleVals = mean(tmpSampleVals{:,:},1);


tmpT = table(tmpSampleVals', 'VariableNames', {'vals'});% here we flip the table tmpSampleVals, so we can add it to tbl_meta_RegionOfInterest)

tbl_data_GEN_RegionOfInterest = [tbl_meta_RegionOfInterest, tmpT]; % Now we add the APOE expression values to their matching sample AND their matching region AND their matching patient


meta_CONindx = (strcmp(tbl_data_GEN_RegionOfInterest.groupID, 'control')); % finds which patients are in the control group
meta_PATindx = (strcmp(tbl_data_GEN_RegionOfInterest.groupID, 'affected')); % find which patients are in the affected group

ConVal=tbl_data_GEN_RegionOfInterest.vals(meta_CONindx,:); % tells us the value of APOE in control patients 
PatVal=tbl_data_GEN_RegionOfInterest.vals(meta_PATindx,:); % tells us the value of APOE in affected patients 


% Now we perform a test statistic 
[h,p,ci,stats]=ttest2(ConVal, PatVal);


t_value = stats.tstat;


p_value = p; 

% Now we make a barplot
x = 1:2;
data=[mean(ConVal), mean(PatVal)];
errlow=[mean(ConVal)-std(ConVal) mean(PatVal)-std(PatVal) ];
errhigh=[mean(ConVal)-std(ConVal) mean(PatVal)-std(PatVal) ];
figure, bar(x, data)
hold on
er=errorbar(x,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
xticks(x);
xticklabels({'Control', 'Affected'});
title(sprintf('Gene expression of %s in %s (p = %.3f)', GeneName, BrainRegion, p));

% Now we make a violin plot

plotdata.Controls = ConVal;
plotdata.Patients = PatVal;
figure; violinplot(plotdata);
ylabel('Average expression level');
title(sprintf('%s gene expression in the %s', GeneOfInterest, RegionOfInterest));

end