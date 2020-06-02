
%% Summary: CorrPerGene
%
% Script to call Correlation Per Gene 
%          for
%                a) 16p11.2 Deletion Functional Connectivity (FC) and 
%                b) 22q11.2 Deletion Functional Connectivity (FC)
% 
%%
% Step1: load data: GeneExpression, FC beta maps, and CNVgeneSets
% Step2: compute CorrPerGene 
%        *(use pre-computed Label-Shuffle null FC maps to compute p-value)
%        -> save .xslx file with Corr, pval, pval-fdr, CentilesCorr
%        -> use "scriptR_plot_HistogramCorrPerGene_16p22q.R" to make
%        histogram plots for Correlations 
% Step3: compute stats: a) median Corr per CNVgeneSet 
%                       b) Extreme Rank (mean correlation for CNVgenes > 95th percentile for all genes)
%
%%
% Dependencies: This script uses functions defined in BasicFunc... .m files
% Data: all input files and output files are assumed in "../data" directory       
% 
%% Paths

data_dir = '../data';
plots_dir = '../plots';

%% Gene Expression data

% Get AHBA gene expression data (nROI x nGenes: 64 x 15633)
csv_filename = [data_dir '/df_abagen_roi_x_genes_ROInames_MIST_64.csv' ];
tab_data_mat = readtable(csv_filename,'ReadRowNames',true,'PreserveVariableNames',true);

GeneSymbAll = (tab_data_mat.Properties.VariableNames)';
ROInames = tab_data_mat.Properties.RowNames;
ROIxGeneExp_AHBA_all = table2array(tab_data_mat);

%% CNV gene lists

cell_array_cnv_names = { '16p112'; '22q112';};
[opt_CNVgeneSets] = BasicFunc_load_cnv_gene_lists(data_dir,cell_array_cnv_names);


%% Get Beta maps and RowMean (Nodal) profiles

% Beta maps and RowMean (nodal) profiles
xlsx_filename = [data_dir '/Supplementary_Table_S1.xlsx']; % %xlsx filname
cell_beta_map_names = {'1_16p11.2Del'; '3_22q11.2Del'; } ; %the beta map sheet names

[opt_FCmaps] = BasicFunc_load_FC_maps(xlsx_filename,cell_beta_map_names);
% Note: warning appears due to length of table description in row 1

%% Call Corr for FC RowMeans (Nodal) profiles 

fdr_threshold = 0.05;
flag_Null = 2;  % 1: random permutaiton, 2: Label shuffle based Null FC beta maps; 
cell_FC_beta_map_names = {'16p11del'; '22q11del'; };  %simplified names excluding 1_ or 3_

nMriMaps = length(cell_FC_beta_map_names);
cell_corr = cell(1,nMriMaps);
cell_pval = cell(1,nMriMaps);
cell_pval_fdr = cell(1,nMriMaps);
cell_corr_centiles = cell(1,nMriMaps);

% count number of genes with significant pval (/pval_fdr) <0.05
count_genes_signif_pval = zeros(1,nMriMaps);
count_genes_signif_pval_fdr = zeros(1,nMriMaps);

for i=1:nMriMaps
    
    % load pre-computed null maps or random-permutation indices
    if(flag_Null == 2)  
        % pre-computed nodal profile based on label-shuffle null maps
        load([ data_dir '/cell_null_nodal_profiles_' cell_FC_beta_map_names{i,1} '.mat'],'cell_null_nodal_profiles'); 
        cell_null_maps = cell_null_nodal_profiles;
    else
        % pre-computed 10k random-permutation indices
        load([ data_dir '/cell_random_perm_indices_10k.mat'],'cell_nROI_perms');
        cell_null_maps = cell_nROI_perms;    
    end
    
    % FC beta map nodal profile
    in_mri_profile = opt_FCmaps.cell_nodal_profiles{i,1};    
    
    % compute CorrPerGene and PvalperGene
    [arrayCorrPerGene,arrayPvalperGene] = fComputeCorr_inMRIprofile_v_GeneExp_perGene(ROIxGeneExp_AHBA_all,in_mri_profile,flag_Null,cell_null_maps);
    
    cell_corr{1,i} = arrayCorrPerGene(:);
    cell_pval{1,i} = arrayPvalperGene(:);   
    [~, crit_p, ~, cell_pval_fdr{1,i}]=fdr_bh(arrayPvalperGene(:),fdr_threshold); % fdr
    cell_corr_centiles{1,i} = BasicFunc_compute_centile_per_array(arrayCorrPerGene(:));

    
    count_genes_signif_pval(1,i) = length(find(arrayPvalperGene(:) < 0.05));
    count_genes_signif_pval_fdr(1,i) = length(find(cell_pval_fdr{1,i}(:) < 0.05));
end

% make tables for corr, pval, pval_fdr
table_corr = array2table(cell2mat(cell_corr),'VariableNames',cell_FC_beta_map_names,'RowNames',GeneSymbAll);
table_corr_pval = array2table(cell2mat(cell_pval),'VariableNames',cell_FC_beta_map_names,'RowNames',GeneSymbAll);
table_corr_pval_fdr = array2table(cell2mat(cell_pval_fdr),'VariableNames',cell_FC_beta_map_names,'RowNames',GeneSymbAll);
table_corr_centile = array2table(cell2mat(cell_corr_centiles),'VariableNames',cell_FC_beta_map_names,'RowNames',GeneSymbAll);

% save as excel sheets
xlsx_Filename = [data_dir '/tab_CorrPerGene_Pval_16p22q_MIST64.xlsx'];
sheet_name = 'Corr';
writetable(table_corr,xlsx_Filename,'sheet',sheet_name,'Range','A1','WriteRowNames',true)
sheet_name = 'Pval';
writetable(table_corr_pval,xlsx_Filename,'sheet',sheet_name,'Range','A1','WriteRowNames',true)
sheet_name = 'FDRPval';
writetable(table_corr_pval_fdr,xlsx_Filename,'sheet',sheet_name,'Range','A1','WriteRowNames',true)
sheet_name = 'CentilesCorr';
writetable(table_corr_centile,xlsx_Filename,'sheet',sheet_name,'Range','A1','WriteRowNames',true)


%% Stats: Median Correlation for CNVgeneSets

% compute Median Correlation for 16p11.2 and 22q11.2 CNV gene sets using
% CorrPerGene values for FC16p and FC22q 

nCNVgeneSets = length(opt_CNVgeneSets.cell_array_cnv_names);

array_FC_v_CNVgeneSet_names = {'FC16p_v_16pGene'; 'FC16p_v_22qGene'; 
                    'FC22q_v_16pGene'; 'FC22q_v_22qGene'; };

nIterNull = 10000;
flag_twosided_pval = false;  % pvalue: two tail vs one tail

cell_array_rand = cell(4,1);
array_stat = zeros(4,1);
array_stat_pval = zeros(4,1);

counter =1;
for ind_mri =1:nMriMaps
in_corr = table2array(table_corr(:,ind_mri));
in_stat = in_corr(:); 
    for loop_geneSet =1:nCNVgeneSets

    % find indices for CNVgeneSets from GeneSymAll
    [~, indx_set_CNVgenes, ~] = intersect(string(GeneSymbAll),opt_CNVgeneSets.cell_array_geneSets{loop_geneSet,1});

    [stat_cnv,pval_stat,array_stat_rand] = BasicFunc_compute_median_stat_and_pval(in_stat,indx_set_CNVgenes,nIterNull,flag_twosided_pval);

    cell_array_rand{counter,1} = array_stat_rand;
    array_stat(counter,1) = stat_cnv;
    array_stat_pval(counter,1) = pval_stat;
    
    counter = counter + 1;
    end
end


%% Stats: Extreme ranked CNV genes

% compute mean correlation for CNVgenes > 95th percentile
nIterNull = 10000;
in_prctile_threshold = 95;

cell_array_rand_stat = cell(nMriMaps,1);
array_stat_topk = zeros(nMriMaps,1);
array_stat_pval_topk = zeros(nMriMaps,1);

counter =1;
for ind_mri =1:nMriMaps
    in_corr = table2array(table_corr(:,ind_mri));
    in_stat = in_corr(:); 

    % find indices for CNVgeneSets from GeneSymAll
    [~, indx_set_CNVgenes, ~] = intersect(string(GeneSymbAll),opt_CNVgeneSets.cell_array_geneSets{ind_mri,1});

    % count CNV genes in top decile (>95th percentile of all genes)
    top_n_prctile_all_genes = prctile(in_stat(:),in_prctile_threshold );
    max_k = length(find(in_stat(indx_set_CNVgenes(:)) > top_n_prctile_all_genes ));
    
    [mean_stat,pval_stat,array_mean_stat_rand] = BasicFunc_compute_mean_stat_topk_and_pval(in_stat,max_k,indx_set_CNVgenes,nIterNull);
    
    cell_array_rand_stat{ind_mri,1} = array_mean_stat_rand;
    array_stat_topk(ind_mri,1) = mean_stat;
    array_stat_pval_topk(ind_mri,1) = pval_stat;
    
 end



