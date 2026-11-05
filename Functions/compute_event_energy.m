function energy = compute_event_energy(Xc, Xd, networkLabels, events, lagTR)
%COMPUTE_EVENT_ENERGY Duration-normalized squared L2 energy.
%   energy = compute_event_energy(Xc, Xd, networkLabels, events, lagTR)
%
%   Xc, Xd       : ROI x Time x Subject arrays.
%   networkLabels: one numeric network label per ROI.
%   events       : struct of 1-based event TR index vectors.
%   lagTR        : integer shift in TRs (default 0).
%
%   E(r,s,c) = mean_t(X(r,t,s)^2), with t pooled within condition c.
%   Network and brain energy are unweighted means of ROI energies.

if nargin < 5 || isempty(lagTR), lagTR = 0; end
if ~isnumeric(Xc) || ~isreal(Xc) || ~isnumeric(Xd) || ~isreal(Xd) || ...
        ~isequal(size(Xc), size(Xd)) || any(isinf(Xc(:))) || any(isinf(Xd(:)))
    error('Xc and Xd must be equally sized real numeric arrays without Inf.');
end
[nROI, nTime, nSub] = size(Xc);
[networkIDs, roiGroups] = network_groups(networkLabels, nROI);
[names, indices] = build_conditions(events, nTime, lagTR, true);
nNet = numel(networkIDs);
nCond = numel(names);

EcROI = nan(nROI, nSub, nCond);
EdROI = nan(nROI, nSub, nCond);
EcNet = nan(nNet, nSub, nCond);
EdNet = nan(nNet, nSub, nCond);
EcBrain = nan(1, nSub, nCond);
EdBrain = nan(1, nSub, nCond);

for c = 1:nCond
    idx = indices{c};
    EcROI(:,:,c) = reshape(mean(Xc(:,idx,:).^2, 2, 'omitnan'), nROI, nSub);
    EdROI(:,:,c) = reshape(mean(Xd(:,idx,:).^2, 2, 'omitnan'), nROI, nSub);
    for k = 1:nNet
        rois = roiGroups{k};
        EcNet(k,:,c) = mean(EcROI(rois,:,c), 1, 'omitnan');
        EdNet(k,:,c) = mean(EdROI(rois,:,c), 1, 'omitnan');
    end
    EcBrain(1,:,c) = mean(EcROI(:,:,c), 1, 'omitnan');
    EdBrain(1,:,c) = mean(EdROI(:,:,c), 1, 'omitnan');
end

energy.roi.Ec = EcROI;       % ROI x Subject x Condition
energy.roi.Ed = EdROI;
energy.network.Ec = EcNet;   % Network x Subject x Condition
energy.network.Ed = EdNet;
energy.brain.Ec = EcBrain;   % 1 x Subject x Condition
energy.brain.Ed = EdBrain;
energy.networkIDs = networkIDs;
energy.conditionNames = names;
energy.conditionIndices = indices;
energy.conditionNTR = cellfun(@numel, indices);
energy.eventLagTR = lagTR;
end
