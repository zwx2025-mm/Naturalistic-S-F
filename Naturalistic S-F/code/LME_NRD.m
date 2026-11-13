%% LME_NRD_MainAnalysis.m
% Linear mixed-effects analysis of network-level NRD.
%
% Required workspace inputs (from NRD_MainAnalysis_Struct.m):
%   NRD_Network  Struct with six fields, each Subject x Network:
%                tom_coup, tom_dec, pain_coup, pain_dec,
%                non_coup, non_dec
%   netNames     Network names in the same order as NRD_Network columns
%
% Model (unchanged from the original MSD analysis):
%   Y ~ Component * EventType * Network + (1|SubID) + (1|SubID:Network)
% Fitted with REML; ANOVA uses Satterthwaite degrees of freedom.
% Outputs: lme_NRD, anova_NRD, EMCI_NRD, long_NRD (.mat only).

%% 1. Convert NRD_Network into a long-format table
[nSub, nNet] = size(NRD_Network.tom_coup);

% Each row contains one participant-network observation in the order:
% ToM coupled, ToM decoupled, Pain coupled, Pain decoupled,
% NonEvent coupled, NonEvent decoupled.
Ymatrix = [NRD_Network.tom_coup(:), NRD_Network.tom_dec(:), ...
           NRD_Network.pain_coup(:), NRD_Network.pain_dec(:), ...
           NRD_Network.non_coup(:), NRD_Network.non_dec(:)];

SubID = repelem(repmat((1:nSub)', nNet, 1), 6);
Network = repelem(netNames(:), nSub * 6);
Component = repmat({'coupling'; 'decoupling'; 'coupling'; ...
                    'decoupling'; 'coupling'; 'decoupling'}, nSub * nNet, 1);
EventType = repmat({'ToM'; 'ToM'; 'Pain'; 'Pain'; ...
                    'NonEvent'; 'NonEvent'}, nSub * nNet, 1);
Y = reshape(Ymatrix.', [], 1);

long_NRD = table(SubID, Network, Component, EventType, Y);
long_NRD.SubID = categorical(long_NRD.SubID);
long_NRD.Network = categorical(long_NRD.Network);
long_NRD.Component = categorical(long_NRD.Component);
long_NRD.EventType = categorical(long_NRD.EventType);

% Preserve the original Network reference level / category order.
net_order = {'DefaultMode', 'Frontoparietal', 'Limbic', ...
             'Somatomotor', 'DorsalAttention', 'VentralAttention', 'Visual'};
long_NRD.Network = reordercats(long_NRD.Network, net_order);

%% 2. Fit the NRD mixed-effects model
formula = 'Y ~ Component * EventType * Network + (1|SubID) + (1|SubID:Network)';
lme_NRD = fitlme(long_NRD, formula, 'FitMethod', 'REML');
anova_NRD = anova(lme_NRD, 'DFMethod', 'Satterthwaite');

%% 3. Model-estimated marginal means and 95% confidence intervals
% Fixed-effects predictions only, as in the original code.
net_levels = categories(long_NRD.Network);
comp_levels = categories(long_NRD.Component);
evt_levels = categories(long_NRD.EventType);

[iNet, iComp, iEvt] = ndgrid(1:numel(net_levels), ...
                             1:numel(comp_levels), 1:numel(evt_levels));
nPred = numel(iNet);

predTbl = table();
predTbl.SubID = repmat(long_NRD.SubID(1), nPred, 1);
predTbl.Network = categorical(net_levels(iNet(:)), net_levels);
predTbl.Component = categorical(comp_levels(iComp(:)), comp_levels);
predTbl.EventType = categorical(evt_levels(iEvt(:)), evt_levels);
predTbl.Y = nan(nPred, 1);

[estimatedMean, ci, DF] = predict(lme_NRD, predTbl, ...
    'Conditional', false, 'Prediction', 'curve', ...
    'Alpha', 0.05, 'DFMethod', 'satterthwaite');

if isscalar(DF)
    DF = repmat(DF, nPred, 1);
end

EMCI_NRD = table(cellstr(predTbl.Network), cellstr(predTbl.Component), ...
    cellstr(predTbl.EventType), estimatedMean, ci(:,1), ci(:,2), DF, ...
    'VariableNames', {'Network', 'Component', 'EventType', ...
                      'EstimatedMean', 'CI_Lower', 'CI_Upper', 'DF'});

% Retain the original display order of the estimated-means table.
EMCI_NRD.Network = categorical(EMCI_NRD.Network, net_order, 'Ordinal', true);
EMCI_NRD.Component = categorical(EMCI_NRD.Component, ...
    {'coupling', 'decoupling'}, 'Ordinal', true);
EMCI_NRD.EventType = categorical(EMCI_NRD.EventType, ...
    {'NonEvent', 'ToM', 'Pain'}, 'Ordinal', true);
EMCI_NRD = sortrows(EMCI_NRD, {'EventType', 'Component', 'Network'});

%% 4. Save MATLAB results
%save('LME_NRD_results.mat', 'lme_NRD', 'anova_NRD','EMCI_NRD', 'long_NRD', '-v7.3');

fprintf('NRD mixed-effects model completed');
