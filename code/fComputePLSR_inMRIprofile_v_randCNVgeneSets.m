function [res_PLSR_rand_cnv_genes] = fComputePLSR_inMRIprofile_v_randCNVgeneSets(ROIxGeneExp_AHBA_all,GeneSymbAll,CNVgeneList,in_mri_profile,flag_zscore,nIterNull)

%% Summary: PLSR Specificity for FC row-mean (nodal) profiles vs CNV's GeneExpression
%          re-run PLSR using  pseudo CNVgeneSets (randomly drawn gene sets)
%
% Script to compute PLSR between FC row-mean (nodal) profiles and 
%        CNV's Gene expression profiles 
% 
%  % Input: 1) ROIxGeneExp_AHBA_all: nROI x nGenes geneExpression matrix
%           2) GeneSymAll: list of all gene symbols (15633)
%           3) CNVgeneList: list of CNV region genes
%           4) in_mri_profile: nROI x 1, input MRI (fMRI) alteration
%                              profile
%           5) flag_zscore: true/false -> zscore X and Y matrices else the
%                                       algorithm performs de-meaning only
%           6) nIterNull: number of null iterations (each iteration will
%              draw same number of genes as in CNVregion and re-compute PLSR)
%
%
%   % Output: res_PLSR_rand_cnv_genes 
%             structure containing: 1) PLSR %variance explained for
%             CNVgeneSet for two PLS components
%             2) array_rand_pctvar : array containing PCTVAR (%variance) for
%                the randomly drawn pseudo CNVgeneSets
%             3) P-value: empirical p-value computed using PCTVAR (%variance) 
%                   as test stat by counting the number of times PCTVAR for
%                   the pseudo CNVgeneSets is higher than the original
%                   PCTVAR. 
%
%%
% Reference for PLSR: https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md
%
% 
%% CNV geneSet GeneExpression

[~, indx_CNVgeneSet, ~] = intersect(string(GeneSymbAll),CNVgeneList);

n_genes = length(GeneSymbAll);
ncnv_genes = length(indx_CNVgeneSet);

Y = in_mri_profile(:); % Response variable
X = ROIxGeneExp_AHBA_all(:,indx_CNVgeneSet); % Predictors (GeneExpression)

% z-score if flag_zscore is true else PLSR de-means
if( flag_zscore )
    X=zscore(X);
    Y=zscore(Y);
end

%% PLSR 2 dim
ncomp=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,ncomp);

% cumulative PCTVAR for ncomp
temp=cumsum(100*PCTVAR(2,1:ncomp));
res_PLSR_rand_cnv_genes.PCTVAR_ncomp2 = temp(ncomp);

%% PLSR specificity: run PLSR for nIterNull pseudo CNVgeneSets

array_rand_pctvar = zeros(nIterNull,1);

for i=1:nIterNull
    
    rng(i);
    rand_gene_ids=randperm(n_genes,ncnv_genes);
    Xrand = ROIxGeneExp_AHBA_all(:,rand_gene_ids(:)); % Predictors

    % z-score if true else PLSR de-means
    if( flag_zscore )
        Xrand=zscore(Xrand);
    end
    
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xrand,Y,ncomp);
    % cumulative PCTVAR for ncomp
    temp=cumsum(100*PCTVAR(2,1:ncomp));

    array_rand_pctvar(i,1) = temp(ncomp);  
    
end

% return the pseudo CNVgeneSets (random) PCTVAR array
res_PLSR_rand_cnv_genes.array_rand_pctvar = array_rand_pctvar;

%% Empirical p-value

res_PLSR_rand_cnv_genes.PLSR_specificity_pval = length(find(array_rand_pctvar(:) >= res_PLSR_rand_cnv_genes.PCTVAR_ncomp2))/nIterNull ;


end

