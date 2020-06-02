
%% Summary: PLSR geneExpression 
%
% Script to call Partial Least Square Regression (PLSR) 
%          for
%                a) 16p11.2 Deletion Functional Connectivity (FC) and 
%                b) 22q11.2 Deletion Functional Connectivity (FC)
% 
% 
%%
% Step1: load data: GeneExpression, FC beta maps, and CNVgeneSets
% Step2: Call PLSR for FCrowmeans (nodal profile)
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
%%
% Reference: we reffered to the PLSR analysis from Morgan et al. 2019
%  https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md
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

%% Call PLSR for FC RowMeans (Nodal) profiles and Regional Connectivity profiles

fdr_threshold = 0.05;
flag_Null = 2;  % 1: random permutaiton, 2: Label shuffle based Null FC beta maps; 
flag_zscore = false;  % false => de-mean only;  true => zscore (de-mean and scale)
ncomp =2;    % for FC rowmeans (nodal profile) ncomp varies from 1 to 4
             % keep results for ncomp=2

cell_FC_beta_map_names = {'16p11del'; '22q11del'; };  %simplified names excluding 1_ or 3_
array_FC_v_CNVgeneSet_names = {'FC16p_v_16pGene'; 'FC16p_v_22qGene'; 
                    'FC22q_v_16pGene'; 'FC22q_v_22qGene'; };


nMriMaps = length(cell_FC_beta_map_names);
nCNVgeneSets = length(opt_CNVgeneSets.cell_array_cnv_names);

% cell array to save the return structures
cell_res_nodal_profile_v_cnv_genes = cell(nMriMaps,nCNVgeneSets);
cell_res_regional_profile_v_cnv_genes = cell(nMriMaps,nCNVgeneSets);

% cell array to save %var and pval for nodal profiles
cell_pval_nodal_profile_v_cnv_genes = cell(nMriMaps,nCNVgeneSets);
cell_PCTVAR_nodal_profile_v_cnv_genes = cell(nMriMaps,nCNVgeneSets);

% cell array to save %var and pval for regional connectivity profiles
cell_pval_regional_profile_v_cnv_genes = cell(nMriMaps,nCNVgeneSets);
cell_PCTVAR_regional_profile_v_cnv_genes = cell(nMriMaps,nCNVgeneSets);

% save the PLSR FC rowmeans (nodal) profile results  (PCTVAR and pval) in .xlsx file
xlsx_Filename_nodal = [data_dir '/tab_PLSR_Ncomp2_NodalFC_PctVar_Pval_MIST64.xlsx'];

% save the PLSR regional connectivity results in .xlsx file
xlsx_Filename_regional = [data_dir '/tab_PLSR_Ncomp2_RegionalFC_PctVarPvalFDR_MIST64.xlsx'];

counter =1;
% loop over FC beta maps
for i=1:nMriMaps
    
    % load pre-computed null maps or random-permutation indices
    if(flag_Null == 2)  
        % pre-computed nodal profile based on label-shuffle null maps
        load([ data_dir '/cell_null_nodal_profiles_' cell_FC_beta_map_names{i,1} '.mat'],'cell_null_nodal_profiles'); 
        load([ data_dir '/cell_null_beta_maps_ROIxROI_' cell_FC_beta_map_names{i,1} '.mat'],'cell_null_beta_maps_ROIxROI'); 
        cell_null_beta_maps = cell_null_beta_maps_ROIxROI;
        cell_null_maps = cell_null_nodal_profiles;
    else
        % pre-computed 10k random-permutation indices
        load([ data_dir '/cell_random_perm_indices_10k.mat'],'cell_nROI_perms');
        cell_null_maps = cell_nROI_perms;  
        cell_null_beta_maps = cell_nROI_perms; 
    end
    
    % FC beta map nodal profile
    in_nodal_profile = opt_FCmaps.cell_nodal_profiles{i,1};
    
    % FC beta map (regional connectivity profiles)
    in_BetaMap = opt_FCmaps.cell_beta_maps_ROIxROI{i,1};

    % loop over CNVgeneSets
    for j=1:nCNVgeneSets

        % CNV gene list
        CNVgeneList = opt_CNVgeneSets.cell_array_geneSets{j,1};  
        
        % PLSR for FC rowmeans(nodal) profiles
        [res_PLSR_nodal] = fComputePLSR_inMRIprofile_v_CNVgeneExp(ROIxGeneExp_AHBA_all,GeneSymbAll,CNVgeneList,in_nodal_profile,flag_zscore,flag_Null,cell_null_maps);
        
        cell_res_nodal_profile_v_cnv_genes{i,j} = res_PLSR_nodal;        
        cell_pval_nodal_profile_v_cnv_genes{i,j} = res_PLSR_nodal.Pval_nPLScomp(ncomp);
        cell_PCTVAR_nodal_profile_v_cnv_genes{i,j} = res_PLSR_nodal.PCTVAR_nPLScomp(ncomp);

        % PLSR for FC regional connectivity profiles
        [res_PLSR_regional] = fComputePLSR_regionalConnectivity_vs_CNVgeneExp(ROIxGeneExp_AHBA_all,GeneSymbAll,CNVgeneList,in_BetaMap,flag_zscore,flag_Null,cell_null_beta_maps);

        cell_res_regional_profile_v_cnv_genes{i,j} = res_PLSR_regional;        
        cell_pval_regional_profile_v_cnv_genes{i,j} = res_PLSR_regional.cell_PLSR_pval;
        cell_PCTVAR_regional_profile_v_cnv_genes{i,j} = res_PLSR_regional.cell_PLSR_Rsq;
        
        % save PLSR for FC regional connectivity profiles in xlsx file
        % make tables for ROInames, PCTVAR, Pval, FDRpval
        PCTVAR = cell2mat(res_PLSR_regional.cell_PLSR_Rsq);
        Pval = cell2mat(res_PLSR_regional.cell_PLSR_pval);
        [~, crit_p, ~, FDRpval]=fdr_bh(Pval(:),fdr_threshold); % fdr
        
        table_plsr = table(ROInames,PCTVAR,Pval,FDRpval);

        sheet_name = array_FC_v_CNVgeneSet_names{counter,1};
        writetable(table_plsr,xlsx_Filename_regional,'sheet',sheet_name,'Range','A1','WriteRowNames',false);
                
        counter = counter + 1;
    end
    
end


% Save results for PLSR for FC rowmeans (nodal) profiles: PCTVAR and pval
tempRowNames = {'FC16p11del'; 'FC22q11del'; };
tempVariableNames = {'Genes16p11'; 'Genes22q11'; };

tab_pval_nodal = array2table(cell2mat(cell_pval_nodal_profile_v_cnv_genes),'VariableNames',tempVariableNames,'RowNames',tempRowNames);
sheet_name = 'PVALnodalPLSRncomp2';
writetable(tab_pval_nodal,xlsx_Filename_nodal,'sheet',sheet_name,'Range','A1','WriteRowNames',true);


tab_PCTVAR_nodal = array2table(cell2mat(cell_PCTVAR_nodal_profile_v_cnv_genes),'VariableNames',tempVariableNames,'RowNames',tempRowNames);
sheet_name = 'PCTVARnodalPLSRncomp2';
writetable(tab_PCTVAR_nodal,xlsx_Filename_nodal,'sheet',sheet_name,'Range','A1','WriteRowNames',true);        


