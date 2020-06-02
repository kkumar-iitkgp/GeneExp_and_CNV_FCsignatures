function [opt_CNVgeneSets] = BasicFunc_load_cnv_gene_lists(data_dir,cell_array_cnv_names)


%% Summary:
% Script to load CNV gene set lists
% 
% Input:  1) data_dir: path to csv files (CNV gene list)
%         2) cell_array_cnv_name : cell array list of cnv names
%
% Output: 1) opt_CNVgeneSets: structure with cell array for gene sets per
%                             CNV and cell_array_cnv_names
% 
%%

% data_dir = '../data';
% cell_array_cnv_names = { '16p112'; '22q112';};

csv_filename_suffix  = 'cnv_genes_cnv';  % see ../data for the CNV gene list files

n_cnvs = length(cell_array_cnv_names);
cell_array_geneSets = cell(n_cnvs,1);

for loop_cnv = 1:n_cnvs
    
    % read from csv file per CNV
    csv_filename = [ data_dir '/' csv_filename_suffix cell_array_cnv_names{loop_cnv,1} '.csv' ];
    tab_data_mat = readtable(csv_filename,'ReadRowNames',false,'ReadVariableNames',false);
   
    cell_array_geneSets{loop_cnv,1} = table2array(tab_data_mat(:,1)) ;
end

    opt_CNVgeneSets.cell_array_cnv_names = cell_array_cnv_names ;
    opt_CNVgeneSets.cell_array_geneSets = cell_array_geneSets;

end