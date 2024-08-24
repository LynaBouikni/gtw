% Clear workspace
clear variables;
prSet(1);  % Set for managing debug or verbosity information 

% Load the first pair of EEG-fMRI data
load_path = 'C:\Users\lynab\OneDrive\Bureau\UAE INTERNSHIP\segmented_pairs\pair_1_daughter_eeg_fmri.mat';  % Adjust the path accordingly
loaded_data = load(load_path);

% Extract the EEG and fMRI data from the loaded file
Xs = loaded_data.segment_pair;  % 'segment_pair' is a cell array with EEG and fMRI data

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

%% Correct Basis Function Selection
% Compute the number of time points in each sequence
ns = cellDim(Xs, 2);

% Choose a more flexible basis for real data (e.g., splines)
num_splines = 5;  % Number of spline basis functions
spline_order = 3;  % Order of the spline
bas = baTems(ns(1), ns, 'spl', [num_splines, spline_order]);

%% Initial Alignment using Unaligned Time Warping (UTW)
aliUtw = utw(Xs, bas, aliT);

%% Perform Pairwise Dynamic Time Warping (PDTW)
aliPdtw = pdtw(Xs, aliUtw, aliT, parDtw);

%% Perform Pairwise Derivative Dynamic Time Warping (PDDTW)
aliPddtw = pddtw(Xs, aliUtw, aliT, parDtw);

%% Perform Pairwise IMW (PIMW)
aliPimw = pimw(Xs, aliUtw, aliT, parPimw, parDtw);

%% Perform Generalized Time Warping (GTW)
aliGtw = gtw(Xs, bas, aliUtw, aliT, parGtw, parCca, parFtw);

%% Show Results
shAliCmp(Xs, Xs, {aliPdtw, aliPddtw, aliPimw, aliGtw}, aliT, parCca, parDtw);

%% Show Basis Functions
shGtwPs(bas{1}.P);  % Show the basis functions used in GTW
