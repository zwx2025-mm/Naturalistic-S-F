function [names, indices] = build_conditions(events, nTime, lagTR, includeFulltime)
%BUILD_CONDITIONS Create event, non-event, and optional full-time TR sets.
%   Input events must be a struct of 1-based TR index vectors.
%   Non-event TRs are defined before shifting. The lag is then applied to
%   both event and non-event sets, matching the source energy calculation.

if nargin < 3 || isempty(lagTR), lagTR = 0; end
if nargin < 4, includeFulltime = false; end
if ~isstruct(events) || numel(events) ~= 1 || isempty(fieldnames(events))
    error('events must be a scalar struct with at least one event category.');
end
if ~isscalar(lagTR) || ~isfinite(lagTR) || lagTR ~= round(lagTR)
    error('lagTR must be an integer number of TRs.');
end
names = fieldnames(events);
if any(strcmpi(names, 'nonevent')) || any(strcmpi(names, 'fulltime'))
    error('Event categories may not be named nonevent or fulltime.');
end

indices = cell(numel(names), 1);
allEvents = [];
for c = 1:numel(names)
    tr = events.(names{c});
    if ~isnumeric(tr) || ~isvector(tr) || isempty(tr) || ...
            any(~isfinite(tr(:))) || any(tr(:) ~= round(tr(:))) || ...
            any(tr(:) < 1) || any(tr(:) > nTime)
        error('Event "%s" needs valid 1-based TR indices.', names{c});
    end
    indices{c} = unique(tr(:).');
    allEvents = [allEvents, indices{c}]; %#ok<AGROW>
end

names{end+1,1} = 'nonevent';
indices{end+1,1} = setdiff(1:nTime, unique(allEvents));

for c = 1:numel(indices)
    tr = indices{c} + lagTR;
    indices{c} = unique(tr(tr >= 1 & tr <= nTime));
    if isempty(indices{c})
        error('Condition "%s" has no valid TRs after lag adjustment.', names{c});
    end
end

if includeFulltime
    names{end+1,1} = 'fulltime';
    indices{end+1,1} = 1:nTime;
end
end
