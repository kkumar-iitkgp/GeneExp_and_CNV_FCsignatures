function [stat_cnv,pval_stat,array_stat_rand] = BasicFunc_compute_median_stat_and_pval(in_stat,in_ia,nIterNull,flag_twosided_pval)

%% Summary:
% Script to compute median stat for CNVgeneSet and pval 
% w.r.t Psuedo CNV sets (random set of nCNVgenes)

%%

nGenesTotal = length(in_stat);
nGenesCNV = length(in_ia);

array_stat_rand = cell(nIterNull,1);
stat_cnv = median(in_stat(in_ia(:)));

for i=1:nIterNull
    rng(i);
    rand_idx = randperm(nGenesTotal,nGenesCNV);    
    array_stat_rand{i,1} = median(in_stat(rand_idx(:)));    
end

% compute p-value
if(flag_twosided_pval)
    pval_stat = length( find( abs(cell2mat(array_stat_rand)) >= abs(stat_cnv)))/nIterNull;
else
     pval_stat1 = length( find( cell2mat(array_stat_rand) >= stat_cnv))/nIterNull;
     pval_stat2 = length( find( cell2mat(array_stat_rand) <= stat_cnv))/nIterNull;
     pval_stat = min(pval_stat1,pval_stat2);
end


end

