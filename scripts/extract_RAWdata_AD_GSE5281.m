clear, clc, close all; %this comments clears all your variables, cleans up your terminal to make a fresh start and closes all figures.

tic

% read tables
%read in the table of the gene expression data 
tbl = readtable('GSE5281_series_matrix.txt','ReadVariableNames',false);

% read in the meta-data table about the name and information about the probes
tbl_probe_metadata = readtable('GPL570-55999_metadata.txt');

% extract gene expression table
%find the index that starts with ID_REF and take all subsequent elements of
%the Table
[~, J] = ismember('ID_REF', tbl.Var1);
tbl_expression = tbl(J+1:end, :);

% extract tissue names
tissueSampleDescriptions = table2cell(tbl(J, 2:end))';

% make a META DATA tissue-subject-region table
% select the line that starts with !Sample_title; you can do this by
% looking for the line that has the member 
[~, J] = ismember('!Sample_title', tbl.Var1);
labelNames = table2cell(tbl(J, 2:end))'; %labelNames are the Names of the brain regions examined in this experiment

%creates empty lists for regionDescriptions, subjectID and groupID of
%number of elements (numel) of labelNames
regionDescriptions = cell(numel(labelNames), 1);
subjectID = cell(numel(labelNames), 1);
groupID = cell(numel(labelNames), 1);
for i = 1:numel(labelNames)
    if strfind(labelNames{i}, ' ')
        tmp = strsplit(labelNames{i}, ' ');
        regionDescriptions{i,1} = tmp{1};
        groupID{i,1} = tmp{2};
        subjectID{i,1} = [tmp{2}, '_', tmp{3}];
    else strfind(labelNames{i}, '_')
        tmp = strsplit(labelNames{i}, '_');
        regionDescriptions{i,1} = tmp{1};
        groupID{i,1} = tmp{2};
        subjectID{i,1} = [tmp{2},'_',tmp{3}];
    end
end
%put all meta data in a meta table tbl_meta
tbl_meta = table(tissueSampleDescriptions, regionDescriptions, groupID, subjectID);
%write down the table as a txt file and save it as a mat file
%output 1
writetable(tbl_meta, 'SamplesMetaData.txt', 'Delimiter','\t');

%make a data table
tbl_expression.Properties.VariableNames = ['probe_id', tissueSampleDescriptions'];
%output 2
writetable(tbl_expression, 'ExpressionData.txt', 'Delimiter','\t');
%save('ExpressionData.mat', 'tbl_expression');
%save('SamplesMetaData.mat', 'tbl_meta');

% make probe-gene table
probe_id = tbl_probe_metadata.ID;
gene_entrez_id = tbl_probe_metadata.ENTREZ_GENE_ID;
gene_symbol = tbl_probe_metadata.GeneSymbol;

tmp = cellfun(@(X) strsplit(X,' /// '), gene_symbol,'UniformOutput',false);
N = cellfun(@(X) numel(X), tmp, 'UniformOutput', true);
gene_symbol_alt = cell(numel(gene_symbol), max(N));
for i = 1:numel(gene_symbol)
    tmp = gene_symbol{i};
    if strfind(tmp, ' /// ')
        tmp1 = strsplit(tmp,' /// ');
        for j = 1:numel(tmp1)
            gene_symbol_alt{i,j} = tmp1{j};
        end
    else
        gene_symbol_alt{i,1} = tmp;
    end
end
gene_symbol = gene_symbol_alt(:, 1);
gene_symbol_alt = gene_symbol_alt(:, 2:end);
tbl_probe_metadata = table(probe_id, gene_symbol, gene_symbol_alt);
tbl_probe_metadata.Properties.VariableNames{1}='probid';
%output 3
writetable(tbl_probe_metadata, 'ProbeNames.txt', 'Delimiter','\t');

toc

