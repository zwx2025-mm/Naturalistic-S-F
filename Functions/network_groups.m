function [networkIDs, roiGroups] = network_groups(networkLabels, nROI)
%NETWORK_GROUPS Create ROI index sets for arbitrary numeric network labels.
if ~isnumeric(networkLabels) || ~isvector(networkLabels) || ...
        numel(networkLabels) ~= nROI || any(~isfinite(networkLabels(:)))
    error('networkLabels must be one finite numeric label per ROI.');
end
labels = networkLabels(:);
networkIDs = unique(labels, 'stable');
roiGroups = cell(numel(networkIDs), 1);
for k = 1:numel(networkIDs)
    roiGroups{k} = find(labels == networkIDs(k));
end
end
