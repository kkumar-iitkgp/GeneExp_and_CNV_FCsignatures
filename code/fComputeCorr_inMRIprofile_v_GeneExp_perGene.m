function [arrayCorrPerGene,arrayPvalPerGene] = fComputeCorr_inMRIprofile_v_GeneExp_perGene(ROIxGeneExp_AHBA_all,in_mri_profile,flag_Null,cell_null_maps,ind_roi)


%% Summary:
% Function to compute Correlation between in_mri_profile and GeneExpression
% per gene
%
%  % Input: 1) ROIxGeneExp_AHBA_all: nROI x nGenes geneExpression matrix
%           2) in_mri_profile: nROI x 1, input MRI (fMRI) alteration
%                              profile
%           3) flag_Null: either "1" -> random permutation of
%                         in_mri_profile (pre-computed re-orderings to be
%                         used to speed up)
%                         or "2" -> Label shuffling based null maps
%           4) cell_null_maps: pre-computed (either the random permutation 
%                              indices if flag_Null == 1
%                              or the label shuffle based null map profiles
%                              if flag_Null ==2
%           5) ind_roi: (OPTIONAL) regional connectivity profile index (to
%                       be used for null Label Shuffle beta map 
%                       regional connectiivty profile), a value between 1
%                       and nROI; value of 0 implies a nodal profile
%
%   % Output: 1) arrayCorrPerGene: Correlation value for each Gene
%                               (correlation between spatial pattern of 
%                               expression per gene and in_mri_profile)
%             2) arrayPvalPerGene: empirical p-value computed for each gene
%                                 using either Random permutation of 
%                                 in_mri_profile or label shuffling based 
%                                 null beta map profiles
%

%% 

% if rowmean (nodal profile) is used; else use roi_ind 
if ~exist('ind_roi','var')
     % ind_roi parameter does not exist, so default it to something
      ind_roi = 0;
end

t_test_sides =2;
corr_type = 'Pearson';   % 'Spearman'

%% Compute Correlation per Gene ( in_mri_profile vs GeneExp per gene)

X = in_mri_profile(:);
Y = ROIxGeneExp_AHBA_all;
arrayCorrPerGene =  corr(X,Y,'Type',corr_type);

%% Null permtation (Label shuffling or in_mri_profile permutation) and p-value

nGenes = size(ROIxGeneExp_AHBA_all,2);
arrayPvalPerGene = zeros(nGenes,1);

% Label Shuffle based Null profiles
if( flag_Null ==2)
    
    nIterNull = length(cell_null_maps);
    cell_Corr_all_gene_scores_mri_shuffle = cell(nIterNull,1);    
    for i=1:nIterNull
        
        if(ind_roi ==0)
            X = cell_null_maps{i,1};  % label shuffed profile 
        else
            X = cell_null_maps{i,1}(:,ind_roi);  % label shuffed regional connectivity profile
        end
        
        cell_Corr_all_gene_scores_mri_shuffle{i,1} = corr(X,Y,'Type',corr_type);
    end   
    
else
% random permutation of in_mri_profile 
% (use pre-computed random permutation indices)    

    nIterNull = length(cell_null_maps);
    cell_Corr_all_gene_scores_mri_shuffle = cell(nIterNull,1);
    
    for i=1:nIterNull
        X = in_mri_profile(cell_null_maps{1,i});  % in_mri_profile order permutation
        cell_Corr_all_gene_scores_mri_shuffle{i,1} = corr(X,Y,'Type',corr_type);
    end   
    
end


%% Get p-value (per Gene)

% cell2mat null distribution 
mat_corr_null = cell2mat(cell_Corr_all_gene_scores_mri_shuffle);

% loop over each Gene
for loop_g=1:nGenes
    temp_null_array = mat_corr_null(:,loop_g);
    temp_null_array(end) = arrayCorrPerGene(1,loop_g);  % replace last value by self (to avoid 0 p-value)

    if( t_test_sides ==2)
        % two sided t-test (absolute values)
        arrayPvalPerGene(loop_g,1) = length( find( abs(temp_null_array) >= abs(arrayCorrPerGene(1,loop_g))))/nIterNull;
    else
         % one sided t-test
        t_test_higher = length( find(temp_null_array >= arrayCorrPerGene(1,loop_g) ))/nIterNull;
        t_test_lower = length( find(temp_null_array <= arrayCorrPerGene(1,loop_g)))/nIterNull;
        arrayPvalPerGene(loop_g,1) = min(t_test_higher,t_test_lower);
    end
end

% replace pvalue == 0 with the minimum empirical p-value
arrayPvalPerGene(arrayPvalPerGene == 0) = 1/nIterNull;

end