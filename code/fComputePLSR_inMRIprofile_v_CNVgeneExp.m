function [res_PLSR_cnv_genes] = fComputePLSR_inMRIprofile_v_CNVgeneExp(ROIxGeneExp_AHBA_all,GeneSymbAll,CNVgeneList,in_mri_profile,flag_zscore,flag_Null,cell_null_maps)

%% Summary: PLSR for FC row-mean (nodal) profiles vs CNV's GeneExpression
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
%           6) flag_Null: either "1" -> random permutation of
%                         in_mri_profile (pre-computed re-orderings to be
%                         used to speed up)
%                         or "2" -> Label shuffling based null maps
%           7) cell_null_maps: pre-computed (either the random permutation 
%                              indices if flag_Null == 1
%                              or the label shuffle based null map profiles
%                              if flag_Null ==2
%
%
%   % Output: res_PLSR_cnv_genes 
%             structure containing: 1) PLSR %variance explained for
%             increasing number of PLS components
%             2) P-value: empirical p-value computed using %variance 
%                   as test stat for CNVgeneSet using either Random permutation of 
%                   in_mri_profile or label shuffling based null beta map profiles
%             3) PLS weights (value per CNV gene): W1 and W2
%             4) PLS predictor scores XS, that is, the PLS components that
%                  are linear combinations of the variables in X.
%                  (XS1 and XS2)

%%
% Reference for PLSR: https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md
%
% 
%% CNV geneSet GeneExpression

nIterNull = length(cell_null_maps);

[~, indx_CNVgeneSet, ~] = intersect(string(GeneSymbAll),CNVgeneList);

Y = in_mri_profile(:); % Response variable
X = ROIxGeneExp_AHBA_all(:,indx_CNVgeneSet); % Predictors (GeneExpression)

% z-score if flag_zscore is true
if( flag_zscore )
    X=zscore(X);
    Y=zscore(Y);
end

%% Permutation testing + impact of number of PLSR components
% permutation testing to assess significance of PLS result as a function of
% the number of components (ncomp) included


max_ncomp = 4;

Rsq_nPLScomp = zeros(max_ncomp,1);
Pval_nPLScomp = zeros(max_ncomp,1);
for ncomp=1:max_ncomp
    
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,ncomp);
temp=cumsum(100*PCTVAR(2,1:ncomp));
Rsquared = temp(ncomp);

array_null_Rsq = zeros(nIterNull,1);
    for j=1:nIterNull
        
        if( flag_Null ==2)       
             Yp = cell_null_maps{j,1};  % label shuffed profile 
        else
            Yp=Y(cell_null_maps{1,j});  % random permutation of Y      
        end

        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp(:),ncomp);

        temp=cumsum(100*PCTVAR(2,1:ncomp));
        array_null_Rsq(j) = temp(ncomp);
    end

Rsq_nPLScomp(ncomp)=Rsquared ;
Pval_nPLScomp(ncomp)=length(find(array_null_Rsq>=Rsquared))/nIterNull ;
end

% PLSR results: %var + p-val for varying nPLScomp from 1 to 4
res_PLSR_cnv_genes.Pval_nPLScomp = Pval_nPLScomp;
res_PLSR_cnv_genes.PCTVAR_nPLScomp = Rsq_nPLScomp;

%% PLSR 2 dim
ncomp=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,ncomp);

res_PLSR_cnv_genes.PCTVAR_2dim = PCTVAR;

%% align PLS components
%align PLS components with desired direction for interpretability 
% REF: https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md

[R1,~]=corr([XS(:,1),XS(:,2)],in_mri_profile(:));

if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

%% Return PLS weights and PLS predictor scores (XS)

% PLS weights (one value per CNV Gene)
res_PLSR_cnv_genes.PLS_W1 = stats.W(:,1);
res_PLSR_cnv_genes.PLS_W2 = stats.W(:,2);

% PLS predictor scores XS
res_PLSR_cnv_genes.XS1 = XS(:,1);
res_PLSR_cnv_genes.XS2 = XS(:,2);


end

