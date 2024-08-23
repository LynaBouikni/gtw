

% Generalized Time Warping (GTW) implementation function.

% This function implements Generalized Time Warping (GTW) for aligning multiple time series.
% It alternates between spatial transformation and temporal warping to align sequences optimally.
% The function is capable of visualizing the alignment process if debugging is enabled and stores 
% the alignment results in a structure (ali).
% The overall approach allows for a flexible and iterative alignment process,
% adapting both the spatial and temporal aspects of the sequences involved.


function ali = gtw(Xs, bas, ali0, aliT, parGtw, parCca, parFtw)


%---------------------------------------------------------


    % Input

    %   Xs       -  sequences, 1 x m (cell), di x ni

    %   Xs is a cell array where each cell contains a time series.

    %   bas      -  warping bases, 1 x m (cell)

    %   ali0     -  initial alignment

    %   aliT     -  ground-truth alignment (can be [])

    %   parGtw   -  parameter for GTW

    %     th     -  stop threshold, {.01}

    %     nItMa  -  maximum iteration number, {100}

    %     debg   -  debug flag, 'y' | {'n'}

    %     fig    -  figure used for debugging, {11}

    %   parCca   -  parameter for CCA. See function cca for more details

    %   parFtw   -  parameter for FTW. See function ftw for more details
    

    % Output

    %   ali      -  alignment structure
    %     objs   -  objective values across iterations, 1 x nIt
    %     its    -  iteration step ids, 1 x nIt
    %     P      -  warping path, l x m
    %     Vs     -  transformation matrices, 1 x m (cell), di x b
    %


    %    ali: The structure that contains the final alignment information.
    %    objs: The objective values at each iteration, used for tracking optimization progress.
    %    its: The iteration steps corresponding to the optimization.
    %    P: The warping paths determined for each sequence.
    %    Vs: The spatial transformation matrices applied to each sequence.


%---------------------------------------------------------

% function parameter
% Extract parameters from parGtw, with default values


%th: Stopping threshold, controls when the optimization process should terminate.

th = ps(parGtw, 'th', .01);

%nItMa: Maximum number of iterations the optimization will run for.

nItMa = ps(parGtw, 'nItMa', 100);

%isDebg: A flag to enable or disable debugging visualizations.

isDebg = psY(parGtw, 'debg', 'n');

%fig: Figure number for any debugging plots.

fig = ps(parGtw, 'fig', 11);


%---------------------------------------------------------

% dimension

m = length(Xs);                             % Number of sequences
dims = cellDim(Xs, 1);                      % Dimensions of each sequence
dimsStr = vec2str(dims, '%d', 'delim', ' ');
ns = cellDim(Xs, 2);                        % Number of time points in each sequence
nsStr = vec2str(ns, '%d', 'delim', ' ');
ks = zeros(1, m);
for j = 1 : m
    ks(j) = size(bas{j}.P, 2);               % Warping basis size for each sequence
end
ksStr = vec2str(ks, '%d', 'delim', ' ');
prIn('gtw', 'dims %s, ns %s, ks %s', dimsStr, nsStr, ksStr); % Print input information


% m: Number of sequences to align.
%dims: Dimensions of each sequence, typically the number of channels or features.
%ns: Number of time points in each sequence.
%ks: Warping basis dimensions, used to adjust the temporal alignment for each sequence.


%---------------------------------------------------------

% homogeneous coordinate

% Add a homogeneous coordinate to the sequences
Xs = homoX(Xs);

%homoX(Xs): This function adds a homogeneous coordinate 
%(usually appending a row of ones) to each sequence in Xs.
% This is often used in transformation operations.

%---------------------------------------------------------

% Initialize axes for debugging if the debug flag is set
if isDebg
    Ax1 = iniAx(fig, 2, 5, [250 * 4, 250 * 5], 'pos', [0 .5 1 .5]);
    Ax2 = iniAx(0, 2, m, [], 'pos', [0 0 1 .5]);
else
    Ax1 = [];
    Ax2 = [];
end

%Ax1, Ax2: Handles for the axes used for plotting debug information. 
% If isDebg is true, these are initialized to show intermediate results.

%---------------------------------------------------------

% Perform coordinate-descent optimization

ali = ali0;                            % Start with the initial alignment
[objs, its] = zeross(1, nItMa);        % Initialize storage for objectives and iteration steps
prCIn('EM', nItMa, .1);                % Print initial status

for nIt = 1 : 2 : nItMa
    prC(nIt);                          % Print current iteration number

    % Spatial transformation
    Ys = seqInp(Xs, ali0.P, parFtw);   % Transform sequences using current alignment
    ali.Vs = mcca(Ys, parCca);         % Apply CCA to the transformed sequences

    objs(nIt) = gtwObj(Xs, ali, parCca, parFtw); % Compute the objective function
    its(nIt) = a2it('spa');            % Record the iteration step as spatial transformation
    if nIt == 1
        parCca.d = size(ali.Vs{1}, 2) - 1; % Set the dimensionality for CCA
    end
    
    % Debugging visualization
    if isDebg
        debg(Ax1, Ax2, nIt, objs, its, Xs, bas, aliT, ali, parCca, parFtw);
    end

    % Temporal warping
    Ys = cellTim(cellTra(ali.Vs), Xs); % Transform sequences with the spatial transformations
    aliF = ftw(Ys, bas, ali0, aliT, parFtw); % Apply Fast Time Warping (FTW)
    ali.P = aliF.P;                    % Update the warping path

    objs(nIt + 1) = gtwObj(Xs, ali, parCca, parFtw); % Update the objective value after warping
    its(nIt + 1) = a2it('tem');        % Record the iteration step as temporal warping
    
    % Debugging visualization
    if isDebg
        debg(Ax1, Ax2, nIt + 1, objs, its, Xs, bas, aliT, ali, parCca, parFtw);
    end

    % Stopping condition: check if the warping path has converged

    if pDif(ali.P, ali0.P) <= th
        break;
    end
    ali0 = ali;                        % Update alignment for the next iteration
end
prCOut(nIt);                           % Print completion status

%Coordinate Descent: The optimization alternates between spatial transformation
%(aligning the sequences spatially) and temporal warping (aligning them temporally).
%Stopping Condition: The loop stops if the alignment change between iterations is less
% than the threshold (th).

%---------------------------------------------------------

% Store final alignment information
ali.alg = 'gtw';                       % Store algorithm name
ali.obj = objs(nIt + 1);               % Store final objective value
ali.objs = objs(1 : nIt);              % Store objective values for all iterations
ali.its = its(1 : nIt);                % Store iteration steps
if ~isempty(aliT)
    ali.dif = pDif(ali.P, aliT.P);     % If ground-truth is available, calculate the difference
end

prOut;                                 % Print final status

%---------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function debg(isDebg, Ax1, Ax2, nIt, objs, its, Xs, bas, aliT, ali, parCca, parDtw)
% Show current status in the optimization of GTW.

if ~isDebg
    return;
end

% Visualization parameters
parPca = st('d', 3, 'homo', 'y');
parMk = st('mkSiz', 3, 'lnWid', 1, 'ln', '-');
parChan = st('nms', {'x', 'y', 'z'}, 'gapWid', 1);
parAx = st('mar', .1, 'ang', [30 80]);

% Plot objective values
shIt(objs(1 : nIt), its(1 : nIt), 'ax', Ax1{1, 1}, 'mkSiz', 7, 'itNms', {'space', 'time'});
title('objective');

%Debug Function: If debugging is enabled, this function visualizes the current state of the
% optimization process, including the objective values, sequence transformations, and warping paths.

%---------------------------------------------------------

% dimension

m = length(Xs);     % Number of sequences to be aligned (number of cells in Xs)
ds = cellDim(Xs, 1);  % Get the dimensionality of each sequence (number of features/channels)
ns = cellDim(Xs, 2);  % Get the number of time points in each sequence
Ps = cellFld(bas, 'cell', 'P');  % Extract the warping paths from the bases
l = size(Ps{1}, 1);  % Get the length of the warping path (number of points in the warping basis)


%---------------------------------------------------------

% original

if nIt == 1
    col = 2;  % Set the column index for visualization
    XXs = pcas(Xs, stAdd(parPca, 'cat', 'n'));  % Perform PCA on the original sequences
    
    % Plot the original sequences in 3D space
    shs(XXs, parMk, parAx, 'ax', Ax1{1, col});
    view([30 80]);
    title('original sequence (x vs y)');
    
    % Plot the original sequences with respect to time
    shChans(XXs, parMk, parChan, 'ax', Ax1{2, col});
    title('original sequence (x y vs time)');
end

%First Iteration (nIt == 1): This block runs only in the first iteration to visualize the original 
%sequences.
%XXs = pcas(Xs, stAdd(parPca, 'cat', 'n')): 
%This applies PCA to reduce the dimensionality of the sequences for visualization.
%shs: Plots the sequences in a 3D space where each sequence is projected onto its principal components.
%shChans: Plots the sequences over time, showing how the original sequences change with respect to time.

%---------------------------------------------------------

% truth
if nIt == 1 && ~isempty(aliT)
    col = 3;                                    % Set the column index for visualization
    Ys = gtwTra(Xs, aliT, parCca, parDtw);      % Apply ground-truth alignment transformation
    YYs = pcas(Ys, parPca);                     % Perform PCA on the aligned sequences
    
    % Plot the sequences after applying the ground-truth alignment
    shs(YYs, parMk, parAx, 'ax', Ax1{1, col});
    title('true sequence (x vs y)');
    
    % Plot the sequences after alignment with respect to time
    shChans(YYs, parMk, stAdd(parChan, 'P', aliT.P), 'ax', Ax1{2, col});
    title('true sequence (x y vs time)');
end

%Ground-Truth Alignment: If ground-truth alignment aliT is available, this block visualizes 
%how the sequences should look after applying this correct alignment.
gtwTra: This function applies the ground-truth temporal and spatial transformations to the sequences.
%YYs = pcas(Ys, parPca): Applies PCA to the sequences after ground-truth transformation for visualization.
%Visualization: Similar to the previous block, the sequences are plotted in 3D space and over time, but this time showing the "ideal" alignment.


%---------------------------------------------------------

% current
[Ys, ~, Us] = gtwTra(Xs, ali, parCca, parDtw);      % Apply the current alignment transformation

if nIt == 1
    set(Ax1{1, 4}, 'UserData', ali);                % Store the current alignment in the plot
    ali0 = ali;                                     % Save the initial alignment for reference
    col = 4;                                        % Set the column index for visualization
else
    ali0 = get(Ax1{1, 4}, 'UserData');              % Retrieve the stored alignment
    col = 5;                                        % Set the column index for the next visualization
                                                    % step
end

YYs = pcas(Ys, parPca);                             % Perform PCA on the currently aligned sequences
shs(YYs, parMk, parAx, 'ax', Ax1{1, col});          % Plot the aligned sequences in 3D space
title(sprintf('step %d sequence (x vs y)', nIt));   % Title indicating the iteration step

shChans(YYs, parMk, stAdd(parChan, 'P', ali.P), 'ax', Ax1{2, col});  % Plot aligned sequences over time
title(sprintf('step %d sequence (x y vs time)', nIt));  % Title indicating the iteration step


%Current Alignment Visualization: Shows how the sequences look after applying the current alignment
% during this iteration.
%gtwTra(Xs, ali, parCca, parDtw): Applies the spatial and temporal transformations according to the 
%current alignment (ali).
%set and get: Stores and retrieves the alignment state between iterations, ensuring that the plots 
%update correctly.
%Visualization: The sequences are plotted similarly to the original and ground-truth cases but now 
%reflect the alignment computed at this specific iteration.


%---------------------------------------------------------

% transformation
for i = 1 : length(Xs)
    shM(Us{i}, 'ax', Ax2{1, i}, 'eq', 'y', 'clMap', 'jet');  % Plot the spatial transformation matrix
    title(['Spatial transformation for sequence ' num2str(i)]);  % Title indicating the sequence number
end

%Spatial Transformation Visualization: For each sequence, this block visualizes the spatial 
%transformation matrix (Us{i}) applied during the current iteration.

%shM: Function to visualize the transformation matrix with a color map (jet).

%Us{i}: Transformation matrix for the i-th sequence, showing how the data is spatially aligned.

%---------------------------------------------------------

% bases & weight
for i = 1 : m
    shP(Ps{i}, 'ax', Ax2{2, i}, 'n', ns(i), 'mkSiz', 0);  % Plot the warping basis for the i-th sequence
    plot(1 : l, ali0.P(:, i), '--k', 'LineWidth', 1);  % Plot the initial warping path
    plot(1 : l, ali.P(:, i), '-r', 'LineWidth', 2);  % Plot the current warping path
    axis square;  % Set the axis to be square
    title(['Temporal warping for sequence ' num2str(i)]);  % Title indicating the sequence number
end

%Temporal Warping Visualization: Visualizes the warping paths (temporal alignments) for each sequence.
%shP(Ps{i}, ...): Plots the warping basis for the i-th sequence.
%plot(1 : l, ali0.P(:, i), '--k', 'LineWidth', 1): Plots the initial warping path as a dashed line.
%plot(1 : l, ali.P(:, i), '-r', 'LineWidth', 2): Plots the current warping path as a solid red line.
%Comparison: This visualization helps to compare how the warping path has changed from the initial 
%to the current iteration.

    %Summary:
    %This part of the code handles the visualization of the alignment process. 
    %It shows how the sequences look before and after alignment, compares them to the ground-truth 
    %alignment (if available), and visualizes both the spatial transformations and temporal warping 
    %paths applied to each sequence. These visualizations are crucial for debugging and understanding
    %how the alignment algorithm progresses and how the sequences are being aligned at each step.