function prSet(lMa)
% Set the promption level.
%
% Input
%   lMa     -  maximum level, 0 | 1 | 2 | ...

% variables
global lPr lMaPr;
global nmPrs ticPrs ticPr0s nRepPrs scaPrs;

% level
lPr = 1;         % Set the current level to 1
lMaPr = lMa;     % Set the maximum level to the input value

% list
nMa = 10;        % Maximum number of entries in lists/arrays
nmPrs = cell(1, nMa);  % Initialize a cell array for storing names/labels
ticPrs = zeros(1, nMa, 'uint64');  % Initialize an array for storing tic values
ticPr0s = zeros(1, nMa, 'uint64'); % Initialize an array for storing initial tic values
[nRepPrs, scaPrs] = zeross(1, nMa); % Initialize arrays for repetition counts and scale factors

