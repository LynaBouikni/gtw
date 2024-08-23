% A demo file for comparing temporal alignment algorithms on your real data.
clear variables;
prSet(1);  % Preset for managing debug or verbosity information 


% Load the data

loaded_data = load('C:\Users\lynab\OneDrive\Bureau\UAE INTERNSHIP\prepared_gtw_data_chunk.mat');
Xs = loaded_data.Xs;  % This contains your EEG-fMRI pairs

% If you don't have a ground-truth alignment, set aliT to empty
aliT = [];

% Parameters for the Algorithms
inp = 'linear'; 
qp = 'matlab'; 
dp = 'c';

parGtw = st('nItMa', 50, 'th', 0);  % Parameters for Generalized Time Warping (GTW)
parDtw = st('nItMa', 50, 'th', 0, 'inp', inp);  % Parameters for Dynamic Time Warping (DTW)
parPimw = st('nItMa', 50, 'th', 0, 'lA', 1, 'lB', 1);  % Parameters for Pairwise IMW (PIMW)
parCca = st('d', 2, 'lams', 0);  % Parameters for Canonical Correlation Analysis (CCA)
parFtw = st('nItMa', 2, 'th', 0, 'lam', 0, 'nor', 'n', 'qp', qp, 'inp', inp);  % Parameters for Fast Time Warping (FTW)

%% Basis for Temporal Warping
% The basis is crucial for defining how the sequences can
% be warped in time

ns = cellDim(Xs, 2);  % Number of time points in each sequence
bas = baTems(length(Xs{1}), ns, 'pol', [3 .4], 'tan', [3 .6 1]);  % Use your basis function implementation

%% Initial Alignment using Unaligned Time Warping (UTW)
aliUtw = utw(Xs, bas, aliT);  % Perform initial alignment

%% Perform Pairwise Dynamic Time Warping (PDTW)
aliPdtw = pdtw(Xs, aliUtw, aliT, parDtw);  % Perform DTW alignment

%% Perform Pairwise Derivative Dynamic Time Warping (PDDTW)
aliPddtw = pddtw(Xs, aliUtw, aliT, parDtw);  % Perform Derivative DTW alignment

%% Perform Pairwise IMW (PIMW)
aliPimw = pimw(Xs, aliUtw, aliT, parPimw, parDtw);  % Perform Pairwise IMW alignment

%% Perform Generalized Time Warping (GTW)
aliGtw = gtw(Xs, bas, aliUtw, aliT, parGtw, parCca, parFtw);  % Perform GTW alignment

%% Show Results
shAliCmp(Xs, Xs, {aliPdtw, aliPddtw, aliPimw, aliGtw}, aliT, parCca, parDtw);  % Compare the alignments

%% Show Basis Functions
shGtwPs(bas{1}.P);  % Show the basis functions used in GTW
