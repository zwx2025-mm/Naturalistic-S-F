
%% ============================================================
% Event and whole-time energy quantification
%   E = ||F(T)||_2^2 / |T| = mean(F(T).^2)
% ============================================================

%% 1. network assignments
netassign = netassignments(:);

%% 2. Parameters

params = struct();
params.eventLagTR = 0;

net_ids = 1:7;
net_names = {'Visual', 'Somatomotor', 'DorsalAttention', ...
             'VentralAttention', 'Limbic', ...
             'Frontoparietal', 'DefaultMode'};

nNet = numel(net_ids);

%% 3. Validate input dimensions

assert(exist('X_c','var') == 1 && exist('X_d','var') == 1, ...
       'X_c and X_d must exist in workspace.');

[nROI, nTime, nSubj] = size(X_c);

assert(isequal(size(X_c), size(X_d)), ...
       'X_c and X_d must have identical dimensions.');

assert(numel(netassign) == nROI, ...
       'Network assignment length must match nROI.');

assert(nTime == 145, ...
       'Expected 145 TRs. Check event indices and input data.');

assert(all(ismember(netassign, net_ids)), ...
       'Some ROI network labels are outside net_ids.');

%% 4. Define event time points

event_idx = struct();

event_idx.pain = [ ...
    10 11, 22:24, 32 33, 39 40, 43 44, ...
    54:58, 67 68, 83 84, 90:92, 95 96, ...
    100:103, 119 120, 124:126, 142 143];

event_idx.tom = [ ...
    3 4, 29 30, 46:49, 59:61, 74:80, ...
    97:99, 106:110, 113:115, ...
    132:134, 136:138];

% Non-event time points (before HRF shifting)
event_all = unique([event_idx.pain, event_idx.tom]);
event_idx.nonevent = setdiff(1:nTime, event_all);

% Full-time condition
event_idx.fulltime = 1:nTime;

condition_names = {'pain', 'tom', 'nonevent', 'fulltime'};
nCond = numel(condition_names);

%% 5. Prepare network ROI indices

roi_idx = cell(nNet, 1);

for k = 1:nNet
    roi_idx{k} = find(netassign == net_ids(k));

    assert(~isempty(roi_idx{k}), ...
           'Network %d contains no ROIs.', net_ids(k));
end

%% 6. Prepare output arrays

% ROI x Subject x Condition
Ec_roi = nan(nROI, nSubj, nCond);
Ed_roi = nan(nROI, nSubj, nCond);

% Network x Subject x Condition
Ec_net = nan(nNet, nSubj, nCond);
Ed_net = nan(nNet, nSubj, nCond);

% 1 x Subject x Condition
Ec_brain = nan(1, nSubj, nCond);
Ed_brain = nan(1, nSubj, nCond);

% Record actual TR indices
condition_idx = struct();
condition_nTR = zeros(1, nCond);

%% 7. Calculate duration-normalized squared L2 energy

for e = 1:nCond

    cname = condition_names{e};
    idx = event_idx.(cname);

    % HRF lag applies only to event-related conditions.
    % Full-time always includes all available TRs.
    if ~strcmp(cname, 'fulltime')
        idx = idx + params.eventLagTR;
    end

    % Remove out-of-range and duplicate indices
    idx = unique(idx(idx >= 1 & idx <= nTime));

    assert(~isempty(idx), ...
           'Condition %s has no valid TRs.', cname);

    % Save actual time indices
    condition_idx.(cname) = idx;
    condition_nTR(e) = numel(idx);

    % ---- ROI-level energy ----
    % Duration-normalized squared L2 norm

    Ec_roi(:,:,e) = reshape( ...
        mean(X_c(:,idx,:).^2, 2, 'omitnan'), ...
        nROI, nSubj);

    Ed_roi(:,:,e) = reshape( ...
        mean(X_d(:,idx,:).^2, 2, 'omitnan'), ...
        nROI, nSubj);

    % ---- Network-level energy ----
    % Average ROI energy within each network

    for k = 1:nNet

        rois = roi_idx{k};

        Ec_net(k,:,e) = mean( ...
            Ec_roi(rois,:,e), 1, 'omitnan');

        Ed_net(k,:,e) = mean( ...
            Ed_roi(rois,:,e), 1, 'omitnan');
    end

    % ---- Whole-brain energy ----
    % Average energy across all ROIs

    Ec_brain(1,:,e) = mean( ...
        Ec_roi(:,:,e), 1, 'omitnan');

    Ed_brain(1,:,e) = mean( ...
        Ed_roi(:,:,e), 1, 'omitnan');

end

%% 8. Pack results

% ROI-level results
roiResults = struct();

roiResults.Ec = Ec_roi;
roiResults.Ed = Ed_roi;
roiResults.condition_names = condition_names;
roiResults.condition_idx = condition_idx;
roiResults.condition_nTR = condition_nTR;

% Network-level results
netResults = struct();

netResults.Ec = Ec_net;
netResults.Ed = Ed_net;
netResults.net_ids = net_ids;
netResults.net_names = net_names;
netResults.roi_idx = roi_idx;
netResults.condition_names = condition_names;

% Whole-brain results
brainResults = struct();

brainResults.Ec = Ec_brain;
brainResults.Ed = Ed_brain;
brainResults.condition_names = condition_names;

% Common metadata
energyInfo = struct();

energyInfo.params = params;
energyInfo.definition = 'Duration-normalized squared L2 norm';
energyInfo.condition_names = condition_names;
energyInfo.condition_idx = condition_idx;
energyInfo.condition_nTR = condition_nTR;
energyInfo.netassign = netassign;



