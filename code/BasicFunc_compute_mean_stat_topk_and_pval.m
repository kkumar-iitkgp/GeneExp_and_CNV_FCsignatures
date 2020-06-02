function [mean_stat,pval_stat,array_mean_stat_rand] = BasicFunc_compute_mean_stat_topk_and_pval(in_stat,max_k,indx_set_CNVgenes,nIterNull)

%% Summary:
% Script to compute the mean of stat of top k genes 
%         out of randomly sampled nCNVgenes

%%

% Stat for the top k values
k =max_k;

nTotalgenes = length(in_stat);
nCNVgenes = length(indx_set_CNVgenes);

% Sort the stat values for CNVgeneSets
temp_sorted_r_cnv =sort(in_stat(indx_set_CNVgenes(:)) ,'descend'); 

% compute mean stat for top k
mean_stat = mean(temp_sorted_r_cnv(1:k,1));

% Null: compute mean stat for top k after randomly sampling nCNVgenes
array_mean_stat_rand = cell(nIterNull,1);
    for i=1:nIterNull
        rng(i);
        rand_idx = randperm(nTotalgenes,nCNVgenes);

        temp1 =sort(in_stat(rand_idx(:)) ,'descend');
        array_mean_stat_rand{i,1} = mean(temp1(1:k)) ;    
    end

% compute p-value
pval_stat = length( find( cell2mat(array_mean_stat_rand) >= mean_stat ))/nIterNull;

end


