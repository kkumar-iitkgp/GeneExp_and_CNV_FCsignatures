function [res_PLSR_cnv_genes] = fComputePLSR_regionalConnectivity_vs_CNVgeneExp(ROIxGeneExp_AHBA_all,GeneSymbAll,CNVgeneList,in_BetaMap,flag_zscore,flag_Null,cell_null_maps)


%% Summary: PLSR for regional connectivity profiles vs CNV's GeneExpression
% 
% Script to compute PLSR between regional connectivity profiles and 
%        CNV's Gene expression profiles 
% 
%  % Input: 1) ROIxGeneExp_AHBA_all: nROI x nGenes geneExpression matrix
%           2) GeneSymAll: list of all gene symbols (15633)
%           3) CNVgeneList: list of CNV region genes
%           4) in_BetaMap: nROI x nROI, FC beta map, call PLSR for each
%                          regional connectivity profiles 
%                          (each column of FC beta map)
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

ncomp=2;   % ncomp for PLSR
nIterNull = length(cell_null_maps);

[~, indx_CNVgeneSet, ~] = intersect(string(GeneSymbAll),CNVgeneList);


%% PLSR per regional connectivity profiles vs GeneExpression for CNV genes

nROI = size(in_BetaMap,1);

% Ensure diagonal values are 1 for FC Beta maps
in_BetaMap(1:(size(in_BetaMap,1)+1):end) = 1; % diagonal = 1

cell_PLSR_pval = cell(nROI,1);
cell_PLSR_Rsq = cell(nROI,1);

cell_XS1_XS2 = cell(nROI,2);
cell_PLS_Weights_W1_W2 = cell(nROI,2);

for loop_r=1:nROI
    
    % X = CNV genes Gene Expression
    X = ROIxGeneExp_AHBA_all(:,indx_CNVgeneSet); % Predictors (CNV genes GeneExpression)
   
    % Y = regional connectivity profile (loop_r column)
    Y = in_BetaMap(:,loop_r);      
    
    % z-score:
    if( flag_zscore )
        X=zscore(X);
        Y=zscore(Y);
    end
        
    %% PLSR weights and XS1 XS2
        
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,ncomp);

     %% align PLS components 
        %store regions' IDs and weights in descending order of weight for both components:
        [R1,p1]=corr([XS(:,1),XS(:,2)],Y(:));

        %align PLS components with desired direction for interpretability 
        if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
            stats.W(:,1)=-1*stats.W(:,1);
            XS(:,1)=-1*XS(:,1);
        end
        if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
            stats.W(:,2)=-1*stats.W(:,2);
            XS(:,2)=-1*XS(:,2);
        end

        % PLS weights (one value per CNV Gene)
        cell_PLS_Weights_W1_W2{loop_r,1} = stats.W(:,1);
        cell_PLS_Weights_W1_W2{loop_r,2} = stats.W(:,2);

        % PLS predictor scores XS
        cell_XS1_XS2{loop_r,1} = XS(:,1);
        cell_XS1_XS2{loop_r,2} = XS(:,2);


        %% Permutation testing  (% variance explained)
        % permutation testing to assess significance of PLS result
        % (ncomp=2)
        % permutation of Y (mri profile)
        
        % cumulative PCTVAR explained
        temp=cumsum(100*PCTVAR(2,1:ncomp));
        Rsquared = temp(ncomp);
        
        array_null_Rsq = zeros(nIterNull,1);
        
        for j=1:nIterNull

            if( flag_Null ==2)       
                 Yp = cell_null_maps{j,1}(:,loop_r);  % label shuffed profile 
            else
                Yp=Y(cell_null_maps{1,j});  % random permutation of Y      
            end

            [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp(:),ncomp);

            temp=cumsum(100*PCTVAR(2,1:ncomp));
            array_null_Rsq(j) = temp(ncomp);
        end
      
        cell_PLSR_pval{loop_r,1} = length(find(array_null_Rsq>=Rsquared))/nIterNull ;
        cell_PLSR_Rsq{loop_r,1} = Rsquared;


end

%res_PLSR_cnv_genes.cnv_genes = cnv_genes;
res_PLSR_cnv_genes.cell_PLSR_pval = cell_PLSR_pval;
res_PLSR_cnv_genes.cell_PLSR_Rsq = cell_PLSR_Rsq ;

res_PLSR_cnv_genes.cell_PLS_Weights_1_2 = cell_PLS_Weights_W1_W2; % PLS weights (one value per CNV Gene)
res_PLSR_cnv_genes.cell_XS1_XS2 =cell_XS1_XS2;  % PLS predictor scores XS


end

