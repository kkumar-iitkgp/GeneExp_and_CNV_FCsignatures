function [opt_FCmaps] = BasicFunc_load_FC_maps(xlsx_filename,cell_beta_map_names)

%% Summary:
% Script to load Functional Connectivity maps from Supplement (.xlsx) file
% output: Structure opt_FCmaps with cell arrays for 
%         1) Row Mean profiles per FC map
%         2) FC beta map (nROIxnROI) => (regional connectivity profiles)
%         3) Beta map names 
%
%%

% %xlsx filname
% xlsx_filename = [data_dir '/Supplementary_Table_S1.xlsx'];
% %the beta map sheet names
%cell_beta_map_names = {'1_16p11.2Del'; '3_22q11.2Del'; } ;

nMaps = length(cell_beta_map_names);
                    
cell_beta_maps_ROIxROI = cell(nMaps,1);   
cell_nodal_profiles = cell(nMaps,1);
                    
for i=1:nMaps
    in_sheetname = cell_beta_map_names{i,1} ;
    tab_data_mat = readtable(xlsx_filename,'FileType','spreadsheet',...
        'Sheet',in_sheetname,'PreserveVariableNames',true);
    
    % Exclude the ROInames and header info (1st col and 2 rows)
    % and convert table2array -> string to num --> cell2mat
    beta_map_ROIxROI = cell2mat(cellfun(@str2num,table2array(tab_data_mat(3:end,2:end)),'un',0));
    
    cell_beta_maps_ROIxROI{i,1} = beta_map_ROIxROI;
    
    % RowMeans (nodal profiles): set diagonal to 0
    beta_map_ROIxROI(1:(size(beta_map_ROIxROI,1)+1):end)= 0;
    nodal_profile = mean(beta_map_ROIxROI,2);
    cell_nodal_profiles{i,1} = nodal_profile(:);
    
end
    
opt_FCmaps.cell_beta_maps_ROIxROI = cell_beta_maps_ROIxROI;
opt_FCmaps.cell_nodal_profiles = cell_nodal_profiles;
opt_FCmaps.cell_beta_map_names = cell_beta_map_names;

end