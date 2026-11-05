function nrd = compute_nrd(Xc, Xd, networkLabels, events, opts)
%COMPUTE_NRD Network-level neural response deviation (NRD = LOO MSD).
%   nrd = compute_nrd(Xc, Xd, networkLabels, events, opts)
%
%   For each contiguous event segment, average each subject's time series
%   across its TRs, yielding a spatial ROI pattern. Compare this pattern
%   with the mean of all OTHER subjects, then average squared differences
%   over ROIs. Average segment-wise MSD values equally within condition.
%
%   opts.eventLagTR    : integer TR shift (default 0).
%   opts.minNonEventTR : minimum consecutive non-event TRs (default 2).
%   opts.minValidROI   : scalar threshold, or [] for original adaptive rule:
%                        max(5, min(30, round(0.2 * number_of_network_ROIs))).
%
%   The source calculation uses segment-mean patterns, not sliding windows.

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isstruct(opts), error('opts must be a struct.'); end
if ~isfield(opts, 'eventLagTR'), opts.eventLagTR = 0; end
if ~isfield(opts, 'minNonEventTR'), opts.minNonEventTR = 2; end
if ~isfield(opts, 'minValidROI'), opts.minValidROI = []; end

if ~isnumeric(Xc) || ~isreal(Xc) || ~isnumeric(Xd) || ~isreal(Xd) || ...
        ~isequal(size(Xc), size(Xd)) || any(isinf(Xc(:))) || any(isinf(Xd(:)))
    error('Xc and Xd must be equally sized real numeric arrays without Inf.');
end
[nROI, nTime, nSub] = size(Xc);
if nSub < 2
    error('NRD requires at least two subjects for leave-one-out templates.');
end
if ~isscalar(opts.minNonEventTR) || ~isfinite(opts.minNonEventTR) || ...
        opts.minNonEventTR < 1 || opts.minNonEventTR ~= round(opts.minNonEventTR)
    error('opts.minNonEventTR must be a positive integer.');
end
if ~isempty(opts.minValidROI) && (~isscalar(opts.minValidROI) || ...
        ~isfinite(opts.minValidROI) || opts.minValidROI < 1 || ...
        opts.minValidROI ~= round(opts.minValidROI))
    error('opts.minValidROI must be [] or a positive integer.');
end

[networkIDs, roiGroups] = network_groups(networkLabels, nROI);
[names, indices] = build_conditions(events, nTime, opts.eventLagTR, false);
nNet = numel(networkIDs);
nCond = numel(names);
segments = cell(nCond, 1);
for c = 1:nCond
    segments{c} = contiguous_segments(indices{c});
    if strcmp(names{c}, 'nonevent')
        keep = cellfun(@numel, segments{c}) >= opts.minNonEventTR;
        segments{c} = segments{c}(keep);
    end
end

nrd.Ec = nan(nSub, nNet, nCond);  % Subject x Network x Condition
nrd.Ed = nan(nSub, nNet, nCond);
nrd.segmentEc = cell(nNet, nCond);  % Each cell: Subject x Segment
nrd.segmentEd = cell(nNet, nCond);
nrd.minValidROI = zeros(nNet, 1);

for k = 1:nNet
    rois = roiGroups{k};
    if isempty(opts.minValidROI)
        minValid = max(5, min(30, round(0.2 * numel(rois))));
    else
        minValid = opts.minValidROI;
    end
    nrd.minValidROI(k) = minValid;
    for c = 1:nCond
        nrd.segmentEc{k,c} = segment_msd(Xc, rois, segments{c}, minValid);
        nrd.segmentEd{k,c} = segment_msd(Xd, rois, segments{c}, minValid);
        nrd.Ec(:,k,c) = mean(nrd.segmentEc{k,c}, 2, 'omitnan');
        nrd.Ed(:,k,c) = mean(nrd.segmentEd{k,c}, 2, 'omitnan');
    end
end

nrd.networkIDs = networkIDs;
nrd.conditionNames = names;
nrd.conditionIndices = indices;
nrd.segmentIndices = segments;
nrd.eventLagTR = opts.eventLagTR;
nrd.minNonEventTR = opts.minNonEventTR;
nrd.definition = 'Unweighted segment mean of leave-one-out spatial MSD';
end

function segments = contiguous_segments(idx)
% Split a TR vector into maximal runs of consecutive TRs.
idx = unique(idx(:).');
cuts = [1, find(diff(idx) > 1) + 1, numel(idx) + 1];
segments = cell(numel(cuts)-1, 1);
for j = 1:numel(segments)
    segments{j} = idx(cuts(j):cuts(j+1)-1);
end
end

function msd = segment_msd(X, rois, segments, minValid)
% Calculate subject-wise leave-one-out MSD for every segment.
nSub = size(X, 3);
nSeg = numel(segments);
msd = nan(nSub, nSeg);
for j = 1:nSeg
    P = reshape(mean(X(rois,segments{j},:), 2, 'omitnan'), numel(rois), nSub);
    isValid = isfinite(P);
    sums = sum(P, 2, 'omitnan');
    counts = sum(isValid, 2);
    for s = 1:nSub
        peerCount = counts - double(isValid(:,s));
        valid = isValid(:,s) & peerCount > 0;
        if sum(valid) < minValid, continue; end
        ref = (sums(valid) - P(valid,s)) ./ peerCount(valid);
        residual = P(valid,s) - ref;
        msd(s,j) = mean(residual.^2);
    end
end
end
