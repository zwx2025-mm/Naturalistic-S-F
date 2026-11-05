%RUN_ANALYSIS Run the three core analyses on user-provided data.
% Prepare data/input.mat with variables SC, X, networkLabels and events.
% Edit only the file locations and optional parameters below.
clear; clc;
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'functions'));

inputFile  = fullfile(rootDir, 'data', 'input.mat');
outputFile = fullfile(rootDir, 'results', 'core_results.mat');

% Optional settings. Set to zero to use the event indices as supplied.
eventLagTR = 0;
nrdOptions = struct('eventLagTR', eventLagTR, ...
                    'minNonEventTR', 2, 'minValidROI', []);

if ~exist(inputFile, 'file')
    error('Input file missing: %s. See README.md for the required variables.', inputFile);
end
S = load(inputFile, 'SC', 'X', 'networkLabels', 'events');
required = {'SC', 'X', 'networkLabels', 'events'};
for i = 1:numel(required)
    if ~isfield(S, required{i})
        error('Missing variable "%s" in %s.', required{i}, inputFile);
    end
end

[Xc, Xd, harmonics] = decompose_harmonics(S.SC, S.X);
energy = compute_event_energy(Xc, Xd, S.networkLabels, S.events, eventLagTR);
nrd = compute_nrd(Xc, Xd, S.networkLabels, S.events, nrdOptions);

outDir = fileparts(outputFile);
if ~exist(outDir, 'dir'), mkdir(outDir); end
save(outputFile, 'Xc', 'Xd', 'harmonics', 'energy', 'nrd', '-v7.3');
fprintf('Saved core results to: %s\n', outputFile);
