% A demo file for comparing temporal alignment algorithms on your real data.

clear variables;
prSet(1);  % Preset for visualization and debugging (you can customize or remove this)
tag = 3;  % This might be used for a specific dataset tag, but we'll leave it as is for now
l = 300;  % Length of the sequences (time points)
m = 3;  % Number of sequences (modalities or subjects)

%% Parameters for the algorithms
inp = 'linear'; 
qp = 'matlab'; 
dp = 'c';

parGtw = st('nItMa', 50, 'th', 0);  % Parameters for Generalized Time Warping (GTW)
parDtw = st('nItMa', 50, 'th', 0, 'inp', inp);  % Parameters for Dynamic Time Warping (DTW)
parPimw = st('nItMa', 50, 'th', 0, 'lA', 1, 'lB', 1);  % Parameters for Pairwise IMW (PIMW)
parCca = st('d', 2, 'lams', 0);  % Parameters for Canonical Correlation Analysis (CCA)
parFtw = st('nItMa', 2, 'th', 0, 'lam', 0, 'nor', 'n', 'qp', qp, 'inp', inp);  % Parameters for Fast Time Warping (FTW)



%% src
wsSrc = toyAliSrc(tag, l, m, 'svL', 1);
[Xs, aliT] = stFld(wsSrc, 'Xs', 'aliT');

%% basis
ns = cellDim(Xs, 2);
bas = baTems(l, ns, 'pol', [3 .4], 'tan', [3 .6 1]);

%% utw (initialization)
aliUtw = utw(Xs, bas, aliT);

%% pdtw
aliPdtw = pdtw(Xs, aliUtw, aliT, parDtw);

%% pddtw
aliPddtw = pddtw(Xs, aliUtw, aliT, parDtw);

%% pimw
aliPimw = pimw(Xs, aliUtw, aliT, parPimw, parDtw);

%% gtw
aliGtw = gtw(Xs, bas, aliUtw, aliT, parGtw, parCca, parFtw);

%% show result
shAliCmp(Xs, Xs, {aliPdtw, aliPddtw, aliPimw, aliGtw}, aliT, parCca, parDtw);

%% show basis
shGtwPs(bas{1}.P);

