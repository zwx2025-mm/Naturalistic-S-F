% GSP_LAPLACIAN_INDIVIDUAL Individual connectome harmonic decomposition.
%
% Inputs
%   fmri_z       - nROI x nTime x nSubj regional fMRI time series.
%                  Must be preprocessed and z-scored across time for each
%                  ROI and subject before calling this function.
%   SC_volnorm   - nROI x nROI x nSubj individual structural connectomes.
%                  Must be nonnegative, symmetric, and volume-normalized
%                  (e.g., streamlines / (volume_i + volume_j)) beforehand.
%                  Graph degree normalization is performed below.
%                  Subject and ROI orders must match fmri_z.
%
% Output (struct GSP)
%   U_all        - nROI x nROI x nSubj: harmonics, low to high frequency
%   Lambda_all   - nROI x nSubj: Laplacian eigenvalues
%   PSD_all      - nROI x nSubj: mean squared harmonic coefficients
%   NN_half      - nSubj x 1: subject-specific low-frequency cutoff index
%   Vlow_all     - nROI x nROI x nSubj: zero-padded low-frequency basis
%   Vhigh_all    - nROI x nROI x nSubj: zero-padded high-frequency basis
%
% The Laplacian, PSD cutoff, and basis construction follow the original code.

[nROI, nTime, nSubj] = size(fmri_z);
[nSC1, nSC2, nSCSubj] = size(SC_volnorm);
assert(nSC1 == nROI && nSC2 == nROI && nSCSubj == nSubj, ...
    'SC_volnorm must be nROI x nROI x nSubj and match fmri_z.');
assert(isnumeric(fmri_z) && isnumeric(SC_volnorm) && ...
    all(isfinite(fmri_z(:))) && all(isfinite(SC_volnorm(:))), ...
    'Inputs must be numeric and contain only finite values.');

fprintf('Data: %d ROIs x %d time points x %d subjects\n', nROI, nTime, nSubj);

% Preallocate subject-level results.
GSP.U_all      = zeros(nROI, nROI, nSubj);
GSP.Lambda_all = zeros(nROI, nSubj);
GSP.PSD_all    = zeros(nROI, nSubj);
GSP.NN_half    = zeros(nSubj, 1);
GSP.Vlow_all   = zeros(nROI, nROI, nSubj);
GSP.Vhigh_all  = zeros(nROI, nROI, nSubj);

for s = 1:nSubj
    W = SC_volnorm(:,:,s);

    % Symmetric degree normalization: Wnorm = D^(-1/2) * W * D^(-1/2).
    deg = sum(W, 2);
    inv_sqrt_deg = zeros(size(deg));
    idx = deg > 0;
    inv_sqrt_deg(idx) = 1 ./ sqrt(deg(idx));
    D_inv_sqrt = diag(inv_sqrt_deg);
    Wnorm = D_inv_sqrt * W * D_inv_sqrt;

    % Normalized graph Laplacian and eigenmodes.
    L = eye(nROI) - Wnorm;
    [U, LambdaL] = eig(L);
    [lambda_vec, indL] = sort(diag(LambdaL));
    U = U(:, indL);

    GSP.U_all(:,:,s)     = U;
    GSP.Lambda_all(:,s)  = lambda_vec;

    % Project z-scored fMRI onto individual structural harmonics.
    X_hat_L = U' * fmri_z(:,:,s);
    PSDs = mean(abs(X_hat_L).^2, 2);
    GSP.PSD_all(:,s) = PSDs;

    % Subject-specific half-energy cutoff (original trapz procedure).
    AUCTOT = trapz(PSDs(1:nROI));
    i = 0;
    AUC = 0;
    while AUC < AUCTOT/2 && i < nROI
        i = i + 1;
        AUC = trapz(PSDs(1:i));
    end
    NN = i;
    GSP.NN_half(s) = NN;

    % Preserve the original reversed and zero-padded harmonic bases.
    M = fliplr(U);
    NNL = nROI - NN;
    Vlow = zeros(size(M));
    Vhigh = zeros(size(M));
    Vhigh(:,1:NNL) = M(:,1:NNL);
    Vlow(:,end-NN+1:end) = M(:,end-NN+1:end);

    GSP.Vlow_all(:,:,s)  = Vlow;
    GSP.Vhigh_all(:,:,s) = Vhigh;

    fprintf('Subject %d/%d: %d low-frequency, %d high-frequency modes\n', ...
        s, nSubj, NN, NNL);
end

