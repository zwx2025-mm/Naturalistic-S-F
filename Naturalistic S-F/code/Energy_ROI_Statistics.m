%% ROI_Statistics.m
% ================================================================
% ROI-level statistics for precomputed duration-normalized energy.
% ================================================================

%% 1. User settings
cfg = struct();
cfg.alpha_q = 0.05;
cfg.nROI = 360;                  % use first 360 ROIs, as original script
cfg.do_zscore_maps = true;
cfg.zscore_method = 'subject_first'; % 'subject_first' | 'group_mean'
cfg.minN = 5;                    % original minimum for ROI t-tests

%% 2. PRECOMPUTED ROI energy 
R = roiResults;
[nROI, nSub, nCond] = size(R.Ec);
roi_use = 1:cfg.nROI;

idxPain = find(strcmp(R.condition_names, 'pain'));
idxTom  = find(strcmp(R.condition_names, 'tom'));
idxNon  = find(strcmp(R.condition_names, 'nonevent'));
idxFull = find(strcmp(R.condition_names, 'fulltime'));
hasEvents = ~isempty(idxPain) && ~isempty(idxTom) && ~isempty(idxNon);

Maps = struct();
Maps.cfg = cfg;
Maps.event_display = struct('pain', 'Pain-Events', ...
    'tom', 'ToM-Events');
Maps.roi_index = (1:cfg.nROI)';
Maps.condition_names = R.condition_names;

fprintf('[ROI] %d ROIs, %d subjects, %d conditions\n', ...
        cfg.nROI, nSub, nCond);

%% 3. Full-time analysis (not part of event-vs-nonevent contrasts)
if ~isempty(idxFull)
    Ec_full = reshape(R.Ec(roi_use,:,idxFull), cfg.nROI, nSub);
    Ed_full = reshape(R.Ed(roi_use,:,idxFull), cfg.nROI, nSub);

    Ec_mean = mean(Ec_full, 2, 'omitnan');
    Ed_mean = mean(Ed_full, 2, 'omitnan');

    Maps.FullTime.GroupMean.Ec_fulltime_Energy = Ec_mean;
    Maps.FullTime.GroupMean.Ed_fulltime_Energy = Ed_mean;
    if cfg.do_zscore_maps
        [Ec_z, Ed_z] = group_energy_zmaps(Ec_full, Ed_full, cfg.zscore_method);
        Maps.FullTime.GroupMeanZ.Method = cfg.zscore_method;
        Maps.FullTime.GroupMeanZ.Ec_fulltime_Energy = Ec_z;
        Maps.FullTime.GroupMeanZ.Ed_fulltime_Energy = Ed_z;
    end
    % Positive difference / t / dz denotes Ed > Ec.
    fullStat = paired_roi_map(Ed_full - Ec_full, cfg.alpha_q, cfg.minN);
    Maps.FullTime.Diff.Ed_vs_Ec_fulltime_Energy = fullStat;
    fprintf('[ROI] Full-time analysis completed.\n');
else
    warning('No fulltime condition found; Full-time analysis skipped.');
end

%% 4. Original event statistics, WITHOUT fulltime in event contrasts
if hasEvents
    runEvents = {'pain', 'tom'};
    Ec_non = reshape(R.Ec(roi_use,:,idxNon), cfg.nROI, nSub);
    Ed_non = reshape(R.Ed(roi_use,:,idxNon), cfg.nROI, nSub);

    for ee = 1:numel(runEvents)
        ev = runEvents{ee};
        switch ev
            case 'pain'
                Ec_in = reshape(R.Ec(roi_use,:,idxPain), cfg.nROI, nSub);
                Ed_in = reshape(R.Ed(roi_use,:,idxPain), cfg.nROI, nSub);
            case 'tom'
                Ec_in = reshape(R.Ec(roi_use,:,idxTom), cfg.nROI, nSub);
                Ed_in = reshape(R.Ed(roi_use,:,idxTom), cfg.nROI, nSub);
        end

        % Group mean maps (ROI x 1)
        mapVals = struct();
        mapVals.Ec_event = mean(Ec_in, 2, 'omitnan');
        mapVals.Ec_nonevent = mean(Ec_non, 2, 'omitnan');
        mapVals.Ed_event = mean(Ed_in, 2, 'omitnan');
        mapVals.Ed_nonevent = mean(Ed_non, 2, 'omitnan');
        fields = fieldnames(mapVals);
        for j = 1:numel(fields)
            nm = fields{j};
            Maps.(ev).GroupMean.([nm '_Energy']) = mapVals.(nm);
        end

        % ROI-wise spatial z-score (descriptive maps only)
        if cfg.do_zscore_maps
            [ecInZ,edInZ] = group_energy_zmaps(Ec_in, Ed_in, cfg.zscore_method);
            [ecNonZ,edNonZ] = group_energy_zmaps(Ec_non, Ed_non, cfg.zscore_method);
            Maps.(ev).GroupMeanZ.Method = cfg.zscore_method;
            zVals = struct('Ec_event',ecInZ, 'Ec_nonevent',ecNonZ, ...
                           'Ed_event',edInZ, 'Ed_nonevent',edNonZ);
            zFields = fieldnames(zVals);
            for j = 1:numel(zFields)
                nm = zFields{j};
                Maps.(ev).GroupMeanZ.([nm '_Energy']) = zVals.(nm);
            end
        end

        % Three ORIGINAL paired contrasts, independent BH-FDR per map
        names = {'Ec_event_vs_nonevent', 'Ed_event_vs_nonevent', ...
                 'Ed_vs_Ec_event'};
        deltas = {Ec_in-Ec_non, Ed_in-Ed_non, Ed_in-Ec_in};
        csvBases = {sprintf('Ec_%s_vs_nonevent',ev), ...
                    sprintf('Ed_%s_vs_nonevent',ev), ...
                    sprintf('Ed_vs_Ec_%s',ev)};
        for j = 1:numel(names)
            stat = paired_roi_map(deltas{j}, cfg.alpha_q, cfg.minN);
            Maps.(ev).Diff.([names{j} '_Energy']) = stat;
        end
        fprintf('[ROI] Event %s completed.\n', ev);
    end
else
    warning('Pain/ToM/Non-event conditions incomplete; event statistics skipped.');
end


%% ======================= LOCAL FUNCTIONS =======================
function stat = paired_roi_map(Delta, alpha_q, minN)
% Identical test family to original: one-sample t-test on paired
% differences, tested independently at each ROI.
[nR, ~] = size(Delta);
stat.t = nan(nR, 1);
stat.p = nan(nR, 1);
stat.q = nan(nR, 1);
stat.meanDiff = nan(nR, 1);
stat.dz = nan(nR, 1);
stat.n = zeros(nR, 1);

for r = 1:nR
    v = Delta(r,:).';
    v = v(isfinite(v));
    stat.n(r) = numel(v);
    if isempty(v), continue; end
    stat.meanDiff(r) = mean(v);
    if numel(v) >= 2 && std(v) > eps
        stat.dz(r) = mean(v) / std(v);
    end
    if numel(v) < minN, continue; end
    [~, p, ~, stats] = ttest(v);
    stat.p(r) = p;
    stat.t(r) = stats.tstat;
end
stat.q = bh_fdr_adjust(stat.p);
stat.sig = stat.q < alpha_q;
end

function q = bh_fdr_adjust(p)
% Benjamini-Hochberg FDR adjusted p-values; finite entries only.
origSize = size(p);
p = p(:);
q = nan(size(p));
valid = isfinite(p);
pv = p(valid);
m = numel(pv);
if m == 0
    q = reshape(q, origSize);
    return;
end
[ps, order] = sort(pv);
qs = ps .* (m ./ (1:m)');
for i = m-1:-1:1
    qs(i) = min(qs(i), qs(i+1));
end
qs = min(qs, 1);
qv = nan(m,1);
qv(order) = qs;
q(valid) = qv;
q = reshape(q, origSize);
end

function [ecZ, edZ] = group_energy_zmaps(Ec, Ed, method)
switch lower(method)
    case 'subject_first'
        ecZ = mean(zscore_each_subject(Ec), 2, 'omitnan');
        edZ = mean(zscore_each_subject(Ed), 2, 'omitnan');
    case 'group_mean'
        ecZ = zscore_roi_vector(mean(Ec, 2, 'omitnan'));
        edZ = zscore_roi_vector(mean(Ed, 2, 'omitnan'));
    otherwise
        error('zscore_method must be subject_first or group_mean.');
end
end

function Z = zscore_each_subject(M)
Z = nan(size(M));
for s = 1:size(M,2)
    Z(:,s) = zscore_roi_vector(M(:,s));
end
end

function z = zscore_roi_vector(x)
x = x(:);
z = nan(size(x));
ok = isfinite(x);
if nnz(ok) < 2, return; end
mu = mean(x(ok));
sd = std(x(ok));
if sd < eps
    z(ok) = 0;
else
    z(ok) = (x(ok)-mu)/sd;
end
end

