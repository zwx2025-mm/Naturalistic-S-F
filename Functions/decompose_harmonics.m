function [Xc, Xd, info] = decompose_harmonics(SC, X)
%DECOMPOSE_HARMONICS Individualized graph-spectral decomposition of fMRI.
%   [Xc, Xd, info] = decompose_harmonics(SC, X)
%
%   SC : ROI x ROI x Subject structural connectomes, or one shared ROI x ROI matrix.
%   X  : ROI x Time x Subject preprocessed fMRI time series.
%   Xc : Low-frequency (structure-coupled) reconstructed signals.
%   Xd : High-frequency (structure-decoupled) reconstructed signals.
%
%   For each subject: L = I - D^(-1/2)*SC*D^(-1/2).
%   The split is the first mode where cumulative trapezoidal PSD area
%   reaches 50% of its total, as in the supplied analysis implementation.
%
%   Method: Preti & Van De Ville, Nat Commun (2019), doi:10.1038/s41467-019-12765-7.

if ~isnumeric(X) || ~isreal(X) || isempty(X) || any(~isfinite(X(:)))
    error('X must be a nonempty, real, finite numeric array.');
end
if ~isnumeric(SC) || ~isreal(SC) || isempty(SC)
    error('SC must be a real numeric array.');
end

X = double(X);
SC = double(SC);
[nROI, nTime, nSub] = size(X);
if nROI < 2 || nTime < 2
    error('X must have at least 2 ROIs and 2 time points.');
end
if size(SC,1) ~= nROI || size(SC,2) ~= nROI
    error('SC and X must have the same ROI count and ordering.');
end
if size(SC,3) ~= 1 && size(SC,3) ~= nSub
    error('SC must have one matrix or one matrix per subject.');
end

Xc = zeros(nROI, nTime, nSub);
Xd = zeros(nROI, nTime, nSub);
info.eigenvectors = zeros(nROI, nROI, nSub);
info.eigenvalues = zeros(nROI, nSub);
info.PSD = zeros(nROI, nSub);
info.cutoff = zeros(nSub, 1);
info.laplacian = 'I - D^(-1/2) * SC * D^(-1/2)';
info.cutoffRule = '50% cumulative trapezoidal integral of PSD across harmonic index';

for s = 1:nSub
    W = SC(:,:,min(s, size(SC,3)));
    if any(~isfinite(W(:))) || any(W(:) < 0)
        error('SC for subject %d must have finite, nonnegative weights.', s);
    end
    if max(max(abs(W - W.'))) > 1e-8 * max(1, max(abs(W(:))))
        error('SC for subject %d must be symmetric.', s);
    end

    degree = sum(W, 2);
    if ~any(degree > 0)
        error('SC for subject %d has no positive connections.', s);
    end
    invSqrtDegree = zeros(nROI, 1);
    positive = degree > 0;
    invSqrtDegree(positive) = 1 ./ sqrt(degree(positive));
    Wnormalized = bsxfun(@times, bsxfun(@times, W, invSqrtDegree), invSqrtDegree.');
    L = eye(nROI) - Wnormalized;

    [U, eigenvalueMatrix] = eig(L);
    [lambda, order] = sort(diag(eigenvalueMatrix), 'ascend');
    U = U(:, order);
    coefficients = U.' * X(:,:,s);
    psd = mean(abs(coefficients).^2, 2);

    totalArea = trapz(psd);
    if ~isfinite(totalArea) || totalArea <= 0
        error('Subject %d has zero or invalid spectral power.', s);
    end
    cutoff = find(cumtrapz(psd) >= totalArea / 2, 1, 'first');

    Xc(:,:,s) = U(:,1:cutoff) * coefficients(1:cutoff,:);
    if cutoff < nROI
        Xd(:,:,s) = U(:,cutoff+1:end) * coefficients(cutoff+1:end,:);
    end

    info.eigenvectors(:,:,s) = U;
    info.eigenvalues(:,s) = lambda;
    info.PSD(:,s) = psd;
    info.cutoff(s) = cutoff;
end
end
