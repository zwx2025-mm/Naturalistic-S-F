%% NRD_MainAnalysis_Struct.m
% Neural representational deviation (NRD; mean squared deviation).
% This script preserves the calculation in Static_MainAnalysis.m.
%
% Required workspace inputs (same ROI and subject order):
%   Xc      ROI x Time x Subject, structure-coupled fMRI signals
%   Xd      ROI x Time x Subject, structure-decoupled fMRI signals
%
% Outputs (structs, each field is a 2D matrix):
%   NRD_ROI.<condition>_<component>     Subject x ROI
%   NRD_Network.<condition>_<component> Subject x Network
%   Conditions: tom, pain, non; components: coup, dec

%% 1. Event indices (original TR indexing)
event_idx.pain = [10:11, 22:24, 32:33, 39:40, 43:44, 54:58, 67:68, 83:84, ...
                  90:92, 95:96, 100:103, 119:120, 124:126, 142:143];
event_idx.tom  = [3:4, 29:30, 46:49, 59:61, 74:80, 97:99, 106:110, ...
                  113:115, 132:134, 136:138];

[nROI, T_len, nSub] = size(X_c);

%% 2. Yeo7 network mapping
netNames = {'Visual', 'Somatomotor', 'DorsalAttention', ...
            'VentralAttention', 'Limbic', 'Frontoparietal', 'DefaultMode'};
nNet = numel(netNames);
netMasks = cell(nNet, 1);
minValidNet = zeros(1, nNet);
for ni = 1:nNet
    netMasks{ni} = find(netassignments(:) == ni)';
    nRoiUse = numel(netMasks{ni});
    minValidNet(ni) = max(5, min(30, round(0.2 * nRoiUse)));
end

%% 3. Consecutive event and non-event segments
all_event_idx = unique([event_idx.tom(:); event_idx.pain(:)])';
non_idx = setdiff(1:T_len, all_event_idx);
indexSets = {event_idx.tom, event_idx.pain, non_idx};
segments = cell(1, 3);

for ci = 1:3
    idx = unique(indexSets{ci}(:)');
    cuts = [1, find(diff(idx) > 1) + 1, numel(idx) + 1];
    segs = cell(1, numel(cuts) - 1);
    for e = 1:numel(segs)
        segs{e} = idx(cuts(e):cuts(e+1)-1);
    end
    if ci == 3  % Keep only non-event segments of at least 2 TRs.
        segs = segs(cellfun(@numel, segs) >= 2);
    end
    segments{ci} = segs;
end

%% 4. ROI-level and network-level NRD
components = {X_c, X_d};
conditionNames = {'tom', 'pain', 'non'};
componentNames = {'coup', 'dec'};
NRD_ROI = struct();
NRD_Network = struct();

for k = 1:2  % Coupled / decoupled
    X = components{k};

    for ci = 1:3  % ToM / Social-Pain / Non-event
        segs = segments{ci};
        nSeg = numel(segs);
        roiSeg = nan(nSub, nROI, nSeg);
        netSeg = nan(nSub, nNet, nSeg);

        for e = 1:nSeg
            trSet = segs{e};
            P = reshape(mean(X(:, trSet, :), 2), [nROI, nSub]);
            sumP = sum(P, 2, 'omitnan');

            for s = 1:nSub
                ps = P(:, s);
                pt = (sumP - ps) ./ (nSub - 1);  % Original LOO reference
                valid = ~isnan(ps) & ~isnan(pt);

                % Squared deviation of each ROI in this segment.
                roiDev = (ps - pt).^2;
                roiDev(~valid) = NaN;
                roiSeg(s, :, e) = roiDev';

                % Network mean per segment (original valid-ROI threshold).
                for ni = 1:nNet
                    roiSet = netMasks{ni};
                    roiValid = valid(roiSet);
                    if sum(roiValid) < minValidNet(ni)
                        continue;
                    end
                    netSeg(s, ni, e) = mean(roiDev(roiSet(roiValid)));
                end
            end
        end

        % Equal-weight average across segments, as in the original code.
        fieldName = [conditionNames{ci}, '_', componentNames{k}];
        NRD_ROI.(fieldName) = mean(roiSeg, 3, 'omitnan');
        NRD_Network.(fieldName) = mean(netSeg, 3, 'omitnan');
    end
end

% Save results from the main pipeline if needed:
% save('NRD_results.mat', 'NRD_ROI', 'NRD_Network', 'event_idx', 'netNames');
