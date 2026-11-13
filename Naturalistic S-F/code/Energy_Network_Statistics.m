%% Network_Statistics.m
% ================================================================
% Network-level statistics using PRECOMPUTED squared-L2 energy.
% ================================================================

%% 1. User settings
alpha = 0.05;
%% 2. Load PRECOMPUTED network energy
N = netResults;

[nNet,nSub,nCond] = size(N.Ec);
net_ids = N.net_ids;
net_names = {'Visual','Somatomotor','DorsalAttention', ...
        'VentralAttention','Limbic','Frontoparietal','DefaultMode'};

idxPain = find(strcmp(N.condition_names,'pain'));
idxTom  = find(strcmp(N.condition_names,'tom'));
idxNon  = find(strcmp(N.condition_names,'nonevent'));
idxFull = find(strcmp(N.condition_names,'fulltime'));

%% 3. ORIGINAL event comparisons (Full-time is excluded)
if hasEvents
    event_names = {'pain','tom','nonevent'};
    eventColumns = [idxPain idxTom idxNon];
    E_c = N.Ec(:,:,eventColumns);
    E_d = N.Ed(:,:,eventColumns);

    % Field names identical to the original event_contrast(2).m.
    contrast_list = { ...
        'pain_vs_tom', 'pain', 'tom'; ...
        'pain_vs_non', 'pain', 'nonevent'; ...
        'tom_vs_non', 'tom', 'nonevent'};

    results_Ec = struct();
    results_Ed = struct();
    for c = 1:size(contrast_list,1)
        cname = contrast_list{c,1};
        ia = find(strcmp(event_names, contrast_list{c,2}));
        ib = find(strcmp(event_names, contrast_list{c,3}));
        results_Ec.(cname) = paired_network_test( ...
            E_c(:,:,ia), E_c(:,:,ib), alpha);
        results_Ed.(cname) = paired_network_test( ...
            E_d(:,:,ia), E_d(:,:,ib), alpha);
    end

    %% 4. ORIGINAL INOUT network-by-subject quantities
    Ec_pain = E_c(:,:,1);
    Ec_tom  = E_c(:,:,2);
    Ec_non  = E_c(:,:,3);
    Ed_pain = E_d(:,:,1);
    Ed_tom  = E_d(:,:,2);
    Ed_non  = E_d(:,:,3);

    % Intentional preservation of original arithmetic average.
    Ec_all = 0.5*(Ec_pain + Ec_tom);
    Ed_all = 0.5*(Ed_pain + Ed_tom);

    INOUT = struct();
    INOUT.net_names = net_names;
    INOUT.nNet = nNet;
    INOUT.nSub = nSub;
    INOUT.event_names = event_names;
    INOUT.Ec = struct('pain',Ec_pain,'tom',Ec_tom, ...
                      'all',Ec_all,'nonevent',Ec_non);
    INOUT.Ed = struct('pain',Ed_pain,'tom',Ed_tom, ...
                      'all',Ed_all,'nonevent',Ed_non);
    INOUT.Gain.Ec = struct( ...
        'pain',Ec_pain-Ec_non, 'tom',Ec_tom-Ec_non, ...
        'all',Ec_all-Ec_non);
    INOUT.Gain.Ed = struct( ...
        'pain',Ed_pain-Ed_non, 'tom',Ed_tom-Ed_non, ...
        'all',Ed_all-Ed_non);
    INOUT.stats_Ec = results_Ec;
    INOUT.stats_Ed = results_Ed;
else
    fprintf('[Network] Event conditions unavailable; event analysis skipped.\n');
end

%% 5. NEW FULL-TIME analysis -- completely separate from event tests
if ~isempty(idxFull)
    Ec_full = reshape(N.Ec(:,:,idxFull),nNet,nSub);
    Ed_full = reshape(N.Ed(:,:,idxFull),nNet,nSub);

    % Positive difference / t / dz denotes Ed > Ec.
    format long g
    fullStat = paired_network_test(Ed_full, Ec_full, alpha);

    FullTime = struct();
    FullTime.Ec = Ec_full;
    FullTime.Ed = Ed_full;
    FullTime.Delta_Ed_vs_Ec = Ed_full-Ec_full;
    FullTime.stats.Ed_vs_Ec = fullStat;
    FullTime.net_names = net_names;
    FullTime.net_ids = net_ids;
    FullTime.nSub = nSub;
    FullTime.alpha = alpha;

    % Standalone paired wide table (nSubject * nNet rows).
    nRow = nSub*nNet;
    SubID_col = zeros(nRow,1);
    NetID_col = zeros(nRow,1);
    NetName_col = cell(nRow,1);
    Ec_col = nan(nRow,1);
    Ed_col = nan(nRow,1);
    Diff_col = nan(nRow,1);
    rr = 0;
    for s = 1:nSub
        for k = 1:nNet
            rr = rr+1;
            SubID_col(rr) = s;
            NetID_col(rr) = k;
            NetName_col{rr} = net_names{k};
            Ec_col(rr) = Ec_full(k,s);
            Ed_col(rr) = Ed_full(k,s);
            Diff_col(rr) = Ed_full(k,s)-Ec_full(k,s);
        end
    end

else
    warning('No fulltime condition; Full-time analysis skipped.');
end


%% ======================= LOCAL FUNCTIONS =======================
function st = paired_network_test(A,B,alpha)
% A,B are Network x Subject. Statistics follow original network code:
% paired, two-sided t-tests + BH-FDR ACROSS NETWORKS, separately
% for every contrast/component.
[nNet, ~] = size(A);
assert(isequal(size(A),size(B)), 'Paired matrix dimensions differ.');
st.t = nan(nNet,1);
st.p = nan(nNet,1);
st.p_fdr = nan(nNet,1);
st.effect = nan(nNet,1); % mean(A - B)
st.dz = nan(nNet,1);
st.n = zeros(nNet,1);
for k = 1:nNet
    x = A(k,:).';
    y = B(k,:).';
    ok = isfinite(x) & isfinite(y);
    d = x(ok)-y(ok);
    st.n(k) = numel(d);
    if isempty(d), continue; end
    st.effect(k) = mean(d);
    if numel(d)>=2 && std(d)>eps
        st.dz(k) = mean(d)/std(d);
    end
    if numel(d)<2, continue; end

    [~,p,~,stats] = ttest(x(ok),y(ok));
    
    t = stats.tstat;
    df = stats.df;
    
    % More numerically stable two-sided p-value
    if isfinite(t)
        p = betainc(df / (df + t^2), df/2, 0.5);
    end
    
    st.p(k) = p;
    st.t(k) = t;
end
st.p_fdr = bh_fdr_adjust(st.p);
st.sig = st.p_fdr < alpha;
end

function q = bh_fdr_adjust(p)
% Equivalent to mafdr(p,'BHFDR',true) on finite p values.
origSize = size(p);
p = p(:);
q = nan(size(p));
valid = isfinite(p);
pv = p(valid);
m = numel(pv);
if m==0
    q = reshape(q,origSize);
    return;
end
[ps,order] = sort(pv);
qs = ps.*(m./(1:m)');
for k = m-1:-1:1
    qs(k) = min(qs(k),qs(k+1));
end
qs = min(qs,1);
qv = nan(m,1);
qv(order) = qs;
q(valid) = qv;
q = reshape(q,origSize);
end

