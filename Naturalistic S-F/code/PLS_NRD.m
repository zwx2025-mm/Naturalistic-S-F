%% PLS_NRD_MainAnalysis.m
% Cross-sectional PLS analysis of network-level NRD and behavioral measures.
% The PLS method and parameters follow PlS_longitudinal_plot.m.
% Longitudinal projection, plotting and CSV export are omitted.
%
% Required MAT files:
%   NRD_results.mat      NRD_Network (struct), netNames (1 x 7 cell)
%       Each NRD_Network field is Subject x Network:
%       tom_coup, tom_dec, pain_coup, pain_dec, non_coup, non_dec.
%   LME_NRD_results.mat  EMCI_NRD (table from LME_NRD_MainAnalysis.m)
%       Required columns: Network, Component, EventType, EstimatedMean.
%   BehaviorData.mat     Behav (table with SubID, Age, Sex, ToM_score,
%                        Empathy_Cognitive, Empathy_Affective)

%% 1. Input and parameters

nComp = 2;
nPerm = 5000;
nBoot = 2000;
rng(20260327);

behaviorVarNames = {'Age', 'ToM_score', 'Empathy_Cognitive', 'Empathy_Affective'};
covariateNames = {'Sex'};
corrType = 'Pearson';

%% 2. Automatically select 14 brain features from relative NRD
% Relative NRD = estimated mean - average estimated mean across 7 networks
% within the SAME component and event condition (not raw NRD dominance).
eventNames = {'ToM', 'Pain'};
eventFields = {'tom', 'pain'};
componentNames = {'coupling', 'decoupling'};
componentFields = {'coup', 'dec'};
netAbbrev = {'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN'};

nNet = numel(netNames);
nSub = size(NRD_Network.tom_coup, 1);
modelMeans = nan(nNet, 2, 2);     % Network x Component x Event
relativeNRD = nan(nNet, 2, 2);

for ei = 1:2
    for ci = 1:2
        for ni = 1:nNet
            idx = strcmp(cellstr(EMCI_NRD.Network), netNames{ni}) & ...
                  strcmp(cellstr(EMCI_NRD.Component), componentNames{ci}) & ...
                  strcmp(cellstr(EMCI_NRD.EventType), eventNames{ei});
            modelMeans(ni, ci, ei) = EMCI_NRD.EstimatedMean(idx);
        end
        relativeNRD(:, ci, ei) = modelMeans(:, ci, ei) - mean(modelMeans(:, ci, ei));
    end
end

Xraw = nan(nSub, nNet * 2);
brainVarNames = cell(nNet * 2, 1);
selectedComponent = cell(nNet * 2, 1);
selectedNetwork = cell(nNet * 2, 1);
selectedEvent = cell(nNet * 2, 1);
selectedCoupRelative = nan(nNet * 2, 1);
selectedDecRelative = nan(nNet * 2, 1);

% Preserve the original variable order: two events per network.
for ni = 1:nNet
    for ei = 1:2
        col = 2 * (ni - 1) + ei;
        if relativeNRD(ni, 1, ei) >= relativeNRD(ni, 2, ei)
            ci = 1;
        else
            ci = 2;
        end
        fieldName = [eventFields{ei}, '_', componentFields{ci}];
        dataMatrix = NRD_Network.(fieldName);
        Xraw(:, col) = dataMatrix(:, ni);
        brainVarNames{col} = [netAbbrev{ni}, '_', componentFields{ci}, '_', eventNames{ei}];
        selectedNetwork{col} = netNames{ni};
        selectedEvent{col} = eventNames{ei};
        selectedComponent{col} = componentNames{ci};
        selectedCoupRelative(col) = relativeNRD(ni, 1, ei);
        selectedDecRelative(col) = relativeNRD(ni, 2, ei);
    end
end

brainSelection = table(selectedNetwork, selectedEvent, selectedComponent, ...
    selectedCoupRelative, selectedDecRelative, brainVarNames, ...
    'VariableNames', {'Network', 'EventType', 'SelectedComponent', ...
                      'CoupledRelativeNRD', 'DecoupledRelativeNRD', 'BrainVariable'});
disp(brainSelection);

%% 3. Match brain and behavioral data by SubID
BrainWide = array2table(Xraw, 'VariableNames', brainVarNames);
BrainWide.SubID = (1:nSub)';
T = outerjoin(Behav, BrainWide, 'Keys', 'SubID', 'MergeKeys', true);
T = sortrows(T, 'SubID');

X = table2array(T(:, brainVarNames));
Y = table2array(T(:, behaviorVarNames));
C = table2array(T(:, covariateNames));
keep = all(~isnan(X), 2) & all(~isnan(Y), 2) & all(~isnan(C), 2);

X = X(keep, :);
Y = Y(keep, :);
C = C(keep, :);
SubID = T.SubID(keep);

fprintf('PLS sample: %d participants; %d brain features; %d behavioral measures.\n', ...
    size(X, 1), size(X, 2), size(Y, 2));

%% 4. Regress out sex from both brain and behavioral variables
% Equivalent to the original column-by-column least-squares residuals.
Xdesign = [ones(size(C, 1), 1), C];
XcovBeta = Xdesign \ X;
YcovBeta = Xdesign \ Y;
X = X - Xdesign * XcovBeta;
Y = Y - Xdesign * YcovBeta;

%% 5. Z-score all variables (same normalization as the original code)
Xmu = nanmean(X, 1);
Xsd = nanstd(X, 0, 1);
Xsd(Xsd == 0) = 1;
Xz = bsxfun(@rdivide, bsxfun(@minus, X, Xmu), Xsd);

Ymu = nanmean(Y, 1);
Ysd = nanstd(Y, 0, 1);
Ysd(Ysd == 0) = 1;
Yz = bsxfun(@rdivide, bsxfun(@minus, Y, Ymu), Ysd);

%% 6. Cross-sectional PLS
[XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(Xz, Yz, nComp);

% Orient LV1 so that its brain score is positively correlated with age.
idxAge = find(strcmp(behaviorVarNames, 'Age'), 1);
if corr(XS(:, 1), Y(:, idxAge), 'type', 'Pearson', 'rows', 'complete') < 0
    XS(:, 1) = -XS(:, 1);
    YS(:, 1) = -YS(:, 1);
    XL(:, 1) = -XL(:, 1);
    YL(:, 1) = -YL(:, 1);
    stats.W(:, 1) = -stats.W(:, 1);
end

scoreCorr = nan(nComp, 1);
scoreP = nan(nComp, 1);
for k = 1:nComp
    [scoreCorr(k), scoreP(k)] = corr(XS(:, k), YS(:, k), ...
        'type', 'Pearson', 'rows', 'complete');
end
cumPctX = cumsum(100 * PCTVAR(1, :));
cumPctY = cumsum(100 * PCTVAR(2, :));

%% 7. Permutation test (5000; unchanged from source)
permNull = nan(nPerm, nComp);
for ip = 1:nPerm
    Yperm = Yz(randperm(size(Xz, 1)), :);
    [~, ~, XS_p, YS_p] = plsregress(Xz, Yperm, nComp);
    for k = 1:nComp
        permNull(ip, k) = abs(corr(XS_p(:, k), YS_p(:, k), ...
            'type', 'Pearson', 'rows', 'complete'));
    end
end
permP = nan(nComp, 1);
for k = 1:nComp
    permP(k) = (sum(permNull(:, k) >= abs(scoreCorr(k))) + 1) / (nPerm + 1);
end

%% 8. Bootstrap loadings / saliences (2000; unchanged from source)
pX = size(Xz, 2);
pY = size(Yz, 2);
W_x_boot = nan(nBoot, pX, nComp);
W_y_boot = nan(nBoot, pY, nComp);
L_x_boot = nan(nBoot, pX, nComp);
L_y_boot = nan(nBoot, pY, nComp);

for ib = 1:nBoot
    idxBoot = randsample(size(Xz, 1), size(Xz, 1), true);
    [XL_b, YL_b, ~, ~, ~, ~, ~, stats_b] = ...
        plsregress(Xz(idxBoot, :), Yz(idxBoot, :), nComp);

    WbX = stats_b.W;
    WbY = YL_b;         % Original code uses Y loading as Y weight proxy.
    LbX = XL_b;
    LbY = YL_b;
    for k = 1:nComp
        if corr(WbX(:, k), stats.W(:, k), 'rows', 'complete') < 0
            WbX(:, k) = -WbX(:, k);
            WbY(:, k) = -WbY(:, k);
            LbX(:, k) = -LbX(:, k);
            LbY(:, k) = -LbY(:, k);
        end
    end
    W_x_boot(ib, :, :) = WbX;
    W_y_boot(ib, :, :) = WbY;
    L_x_boot(ib, :, :) = LbX;
    L_y_boot(ib, :, :) = LbY;
end

bootRes.W_x_mean = squeeze(mean(W_x_boot, 1));
bootRes.W_x_se   = squeeze(std(W_x_boot, 0, 1));
bootRes.W_x_z    = bootRes.W_x_mean ./ max(bootRes.W_x_se, eps);
bootRes.W_y_mean = squeeze(mean(W_y_boot, 1));
bootRes.W_y_se   = squeeze(std(W_y_boot, 0, 1));
bootRes.W_y_z    = bootRes.W_y_mean ./ max(bootRes.W_y_se, eps);
bootRes.L_x_mean = squeeze(mean(L_x_boot, 1));
bootRes.L_x_se   = squeeze(std(L_x_boot, 0, 1));
bootRes.L_x_z    = bootRes.L_x_mean ./ max(bootRes.L_x_se, eps);
bootRes.L_y_mean = squeeze(mean(L_y_boot, 1));
bootRes.L_y_se   = squeeze(std(L_y_boot, 0, 1));
bootRes.L_y_z    = bootRes.L_y_mean ./ max(bootRes.L_y_se, eps);

%% 9. LV summary and associations with behavioral variables
LV = (1:nComp)';
LVSummary = table(LV, scoreCorr, scoreP, cumPctX(:), cumPctY(:), permP, ...
    'VariableNames', {'LV', 'ScoreCorr', 'ScoreCorr_p', ...
                      'CumExplainedX_pct', 'CumExplainedY_pct', 'Permutation_p'});

scoreBehavR = nan(nComp, pY);
scoreBehavP = nan(nComp, pY);
for k = 1:nComp
    for j = 1:pY
        [scoreBehavR(k, j), scoreBehavP(k, j)] = corr(XS(:, k), Y(:, j), ...
            'type', corrType, 'rows', 'complete');
    end
end

%% 10. Save results 
PLS = struct();
PLS.SubID = SubID;
PLS.brainVarNames = brainVarNames;
PLS.brainSelection = brainSelection;
PLS.relativeNRD = relativeNRD;
PLS.behaviorVarNames = behaviorVarNames;
PLS.covariateNames = covariateNames;
PLS.netNames = netNames;
PLS.nComp = nComp;
PLS.nPerm = nPerm;
PLS.nBoot = nBoot;
PLS.rngSeed = 20260327;
PLS.XcovBeta = XcovBeta;
PLS.YcovBeta = YcovBeta;
PLS.Xmu = Xmu;
PLS.Xsd = Xsd;
PLS.Ymu = Ymu;
PLS.Ysd = Ysd;
PLS.XL = XL;
PLS.YL = YL;
PLS.XS = XS;
PLS.YS = YS;
PLS.BETA = BETA;
PLS.PCTVAR = PCTVAR;
PLS.MSE = MSE;
PLS.W = stats.W;
PLS.LVSummary = LVSummary;
PLS.permNull = permNull;
PLS.bootRes = bootRes;
PLS.scoreBehavR = scoreBehavR;
PLS.scoreBehavP = scoreBehavP;

save('PLS_NRD_results.mat', 'PLS', '-v7.3');
disp(LVSummary);
fprintf('Cross-sectional PLS completed: PLS_NRD_results.mat\n');
