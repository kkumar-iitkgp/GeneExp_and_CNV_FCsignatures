
%% Summary
% Script for reproducing results for association between 
%  AHBA Gene expression spatial patterns and FC signatures for 16p11.2
%  deletion and 22q11.2 deletion.

%%
% Kuldeep Kumar
% 31 May 2020
% 
% Please cite: Moreau, Clara, et al. "Neuropsychiatric mutations delineate functional brain connectivity dimensions contributing to autism and schizophrenia." BioRxiv (2019): 862615.
% 
%


%% Step 1. PLSR

cd code

script1_call_PLSR_nodal_and_regional;  %Elapsed time is 343.963346 seconds.


%% Step 2. Correlation per Gene

script2_call_CorrPerGene;   % Elapsed time is 109.778165 seconds.

%% Step 3. Plot CorrPerGene 

% Run R-script:  "scriptR_plot_HistogramCorrPerGene_16p22q.R" 


