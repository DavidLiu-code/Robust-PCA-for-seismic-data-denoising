%% rpca_tapered_reproduction.m
% Single-file MATLAB reproduction of Cheng, Chen & Sacchi (2015):
% f-x-y RPCA / PCP denoising for seismic erratic-noise attenuation.
%
% Tapered reproduction workflow:
%   1) Replaces hard frequency cutoffs with a raised-cosine spectral taper.
%      This strongly reduces time-domain Gibbs/ringing artifacts.
%   2) Uses time zero-padding before FFT and crops after inverse FFT.
%      This gives smoother frequency sampling and less circular leakage.
%   3) Reports Q against both the original clean data and the tapered clean
%      reference. The latter is the fair reference if a band-limited/tapered
%      reconstruction is used.
%   4) Keeps the RPCA + optional rank projection workflow, but applies the
%      same smooth spectral taper to bandpass, PCA, RPCA and sparse output.
%
% Recommended first run:
%   run('src/rpca_tapered_reproduction.m')
%
% For real data, set USER_DATA_MAT and DATA_VAR_NAME below.
%
% Author: generated for Dawei Liu
% Date: 2026-06-05

clear; close all; clc;

%% ---------------- User parameters ----------------
USER_DATA_MAT = '';       % e.g., 'my_seismic_patch.mat'; leave empty for synthetic
DATA_VAR_NAME = 'data';   % variable in USER_DATA_MAT, size [nt,nx] or [nt,nx,ny]

dt = 0.002;

% Smooth frequency taper. Avoid hard fmin/fmax cutoffs.
% Passband is [flow_pass, fhigh_pass]; transition zones are
% [flow_stop, flow_pass] and [fhigh_pass, fhigh_stop].
freq_taper.flow_stop  = 0;    % Hz. Keep 0 to avoid sharp low-cut ringing.
freq_taper.flow_pass  = 2;    % Hz. Frequencies above this are fully kept.
freq_taper.fhigh_pass = 45;   % Hz. Fully kept up to this frequency.
freq_taper.fhigh_stop = 75;   % Hz. Smoothly tapered to zero here.
freq_taper.min_weight_to_process = 0.02;  % skip tiny transition weights

% FFT padding. 2 means roughly double length before FFT.
pad_factor = 2;

% f-x-y PCA/eigenimage rank filtering baseline
rank_pca = 3;

% RPCA/PCP parameters
% lambda = lambda_scale / sqrt(max(nx,ny)). Smaller lambda sends more energy
% to sparse S; larger lambda preserves more energy in low-rank L.
lambda_scale = 0.70;
maxIter      = 800;
tol          = 1e-6;
verbose      = false;

% Hybrid improvement: RPCA removes sparse erratic traces first, then rank
% projection suppresses residual Gaussian/no incoherent energy.
% Set rank_cap = [] to disable rank projection.
rank_cap = 4;

% Synthetic sweep. Clean data are available only for synthetic tests.
DO_SYNTHETIC_SWEEP = true;
lambda_grid = [0.35 0.50 0.70 0.90 1.10];
rank_grid   = {2,3,4,5,6,8,[]};

% Synthetic noise controls
snr_gaussian_db  = 12;
erratic_fraction = 0.08;
erratic_amp      = 3.0;
rng_seed         = 2025;

%% ---------------- Load or synthesize clean data ----------------
using_synthetic = isempty(USER_DATA_MAT);
if ~using_synthetic
    tmp = load(USER_DATA_MAT);
    if ~isfield(tmp, DATA_VAR_NAME)
        error('Variable "%s" not found in %s.', DATA_VAR_NAME, USER_DATA_MAT);
    end
    clean = double(tmp.(DATA_VAR_NAME));
    if ismatrix(clean)
        clean = reshape(clean, size(clean,1), size(clean,2), 1);
    end
    clean = clean ./ (max(abs(clean(:))) + eps);
    fprintf('Loaded user data: nt=%d, nx=%d, ny=%d\n', size(clean,1), size(clean,2), size(clean,3));
else
    rng(rng_seed);
    nt = round(1.0/dt) + 1;
    nx = 30;
    ny = 30;
    clean = make_lowrank_fxy_synthetic(nt, nx, ny, dt);
    clean = clean ./ (max(abs(clean(:))) + eps);
    fprintf('Created synthetic 3D data: nt=%d, nx=%d, ny=%d\n', nt, nx, ny);
end

[nt,nx,ny] = size(clean);
rng(rng_seed);
[noisy, gaussian_noise, erratic_noise] = add_gaussian_and_erratic_noise( ...
    clean, snr_gaussian_db, erratic_fraction, erratic_amp);

% Tapered clean reference. This is useful because any band-limited/tapered
% reconstruction cannot exactly equal the unfiltered clean trace.
clean_tapered = fxy_taper_reconstruct(clean, dt, freq_taper, pad_factor);

%% ---------------- Denoising ----------------
fprintf('\nRunning tapered bandpass-only reconstruction...\n');
bandpass_only = fxy_taper_reconstruct(noisy, dt, freq_taper, pad_factor);

fprintf('Running tapered f-x-y eigenimage filtering, rank=%d...\n', rank_pca);
pca_denoised = fxy_rank_filter_tapered(noisy, dt, freq_taper, pad_factor, rank_pca);

fprintf('Running tapered enhanced f-x-y RPCA/PCP denoising...\n');
opts = struct('lambda_scale',lambda_scale,'maxIter',maxIter,'tol',tol, ...
              'verbose',verbose,'rank_cap',rank_cap);
[rpca_denoised, rpca_sparse_noise, info] = fxy_rpca_pcp_denoise_tapered(noisy, dt, freq_taper, pad_factor, opts);

%% ---------------- Optional synthetic parameter sweep ----------------
sweep_table = [];
best = struct('Q',-Inf,'lambda_scale',lambda_scale,'rank_cap',rank_cap, ...
              'denoised',rpca_denoised,'sparse',rpca_sparse_noise,'info',info);

if using_synthetic && DO_SYNTHETIC_SWEEP
    fprintf('\nRunning synthetic parameter sweep using tapered-clean reference...\n');
    row = 0;
    for il = 1:numel(lambda_grid)
        for ir = 1:numel(rank_grid)
            ls = lambda_grid(il);
            rc = rank_grid{ir};
            opts_tmp = struct('lambda_scale',ls,'maxIter',maxIter,'tol',tol, ...
                              'verbose',false,'rank_cap',rc);
            [tmp_den,tmp_sp,tmp_info] = fxy_rpca_pcp_denoise_tapered(noisy, dt, freq_taper, pad_factor, opts_tmp);
            q_tmp = quality_db(clean_tapered, tmp_den);
            row = row + 1;
            if isempty(rc), rc_print = NaN; else, rc_print = rc; end
            sweep_table(row,:) = [ls, rc_print, q_tmp, mean(tmp_info.ranks), mean(tmp_info.ranks_after)]; %#ok<SAGROW>
            fprintf('  lambda_scale=%4.2f, rank_cap=%4s, Q_tapered=%6.2f dB, mean rank %.2f -> %.2f\n', ...
                ls, rankcap_to_string(rc), q_tmp, mean(tmp_info.ranks), mean(tmp_info.ranks_after));
            if q_tmp > best.Q
                best.Q = q_tmp;
                best.lambda_scale = ls;
                best.rank_cap = rc;
                best.denoised = tmp_den;
                best.sparse = tmp_sp;
                best.info = tmp_info;
            end
        end
    end
    fprintf('Best synthetic setting: lambda_scale=%.2f, rank_cap=%s, Q_tapered=%.2f dB\n', ...
        best.lambda_scale, rankcap_to_string(best.rank_cap), best.Q);

    rpca_denoised = best.denoised;
    rpca_sparse_noise = best.sparse;
    info = best.info;
    lambda_scale = best.lambda_scale;
    rank_cap = best.rank_cap;
end

%% ---------------- Metrics ----------------
q_noisy_full = quality_db(clean, noisy);
q_bp_full    = quality_db(clean, bandpass_only);
q_pca_full   = quality_db(clean, pca_denoised);
q_rpca_full  = quality_db(clean, rpca_denoised);

q_noisy_tap = quality_db(clean_tapered, noisy);
q_bp_tap    = quality_db(clean_tapered, bandpass_only);
q_pca_tap   = quality_db(clean_tapered, pca_denoised);
q_rpca_tap  = quality_db(clean_tapered, rpca_denoised);

fprintf('\nQ_full = 10 log10(||clean||_F^2 / ||estimate-clean||_F^2)\n');
fprintf('Noisy data       : %.2f dB\n', q_noisy_full);
fprintf('Taper only       : %.2f dB\n', q_bp_full);
fprintf('Eigenimage PCA   : %.2f dB\n', q_pca_full);
fprintf('Enhanced RPCA    : %.2f dB\n', q_rpca_full);

fprintf('\nQ_tapered_ref = 10 log10(||clean_tapered||_F^2 / ||estimate-clean_tapered||_F^2)\n');
fprintf('Noisy data       : %.2f dB\n', q_noisy_tap);
fprintf('Taper only       : %.2f dB\n', q_bp_tap);
fprintf('Eigenimage PCA   : %.2f dB\n', q_pca_tap);
fprintf('Enhanced RPCA    : %.2f dB\n', q_rpca_tap);

fprintf('\nProcessed positive frequency bins: %d\n', info.nfreq_processed);
fprintf('Mean RPCA rank before rank projection: %.2f\n', mean(info.ranks));
fprintf('Mean RPCA rank after  rank projection: %.2f\n', mean(info.ranks_after));
fprintf('Final lambda_scale = %.2f, rank_cap = %s\n', lambda_scale, rankcap_to_string(rank_cap));
fprintf('Taper: %.1f-%.1f-%.1f-%.1f Hz, pad_factor=%g\n', ...
    freq_taper.flow_stop, freq_taper.flow_pass, freq_taper.fhigh_pass, freq_taper.fhigh_stop, pad_factor);

%% ---------------- Display one gather/slice ----------------
iy = max(1, round(ny/2));
t = (0:nt-1) * dt;
plot_scale = max(abs(noisy(:,:,iy)), [], 'all') + eps;

figure('Name','v4 tapered RPCA: common wiggle scale','Color','w','Position',[80 80 1300 760]);
subplot(2,3,1); wiggle_plot_fixed(clean(:,:,iy), t, plot_scale); title('(a) Clean gather');
subplot(2,3,2); wiggle_plot_fixed(noisy(:,:,iy), t, plot_scale); title(sprintf('(b) Noisy, Q=%.2f dB', q_noisy_full));
subplot(2,3,3); wiggle_plot_fixed(bandpass_only(:,:,iy), t, plot_scale); title(sprintf('(c) Taper only, Q=%.2f dB', q_bp_tap));
subplot(2,3,4); wiggle_plot_fixed(pca_denoised(:,:,iy), t, plot_scale); title(sprintf('(d) f-x-y PCA, Q=%.2f dB', q_pca_tap));
subplot(2,3,5); wiggle_plot_fixed(rpca_denoised(:,:,iy), t, plot_scale); title(sprintf('(e) Tapered RPCA, Q=%.2f dB', q_rpca_tap));
subplot(2,3,6); wiggle_plot_fixed(rpca_denoised(:,:,iy)-clean_tapered(:,:,iy), t, plot_scale); title('(f) RPCA error vs tapered clean');

figure('Name','Image display, same color scale','Color','w','Position',[160 160 1250 700]);
clim = [-1 1] * max(abs(noisy(:,:,iy)), [], 'all');
subplot(2,3,1); imagesc(1:nx,t,clean(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('Clean'); colorbar;
subplot(2,3,2); imagesc(1:nx,t,noisy(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('Noisy'); colorbar;
subplot(2,3,3); imagesc(1:nx,t,bandpass_only(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('Taper only'); colorbar;
subplot(2,3,4); imagesc(1:nx,t,pca_denoised(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('PCA'); colorbar;
subplot(2,3,5); imagesc(1:nx,t,rpca_denoised(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('Tapered RPCA'); colorbar;
subplot(2,3,6); imagesc(1:nx,t,rpca_denoised(:,:,iy)-clean_tapered(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('RPCA error vs tapered clean'); colorbar;

figure('Name','Frequency taper','Color','w','Position',[180 180 700 360]);
ntfft_plot = choose_ntfft(nt, pad_factor);
[freqs_plot,w_plot] = make_frequency_taper(ntfft_plot, dt, freq_taper);
plot(freqs_plot(1:floor(ntfft_plot/2)+1), w_plot(1:floor(ntfft_plot/2)+1), 'k', 'LineWidth', 1.5);
grid on; xlabel('Frequency (Hz)'); ylabel('Weight'); title('Raised-cosine spectral taper'); ylim([-0.05 1.05]);

%% ---------------- Save results ----------------
results.clean = clean;
results.clean_tapered = clean_tapered;
results.noisy = noisy;
results.gaussian_noise = gaussian_noise;
results.erratic_noise = erratic_noise;
results.bandpass_only = bandpass_only;
results.pca_denoised = pca_denoised;
results.rpca_denoised = rpca_denoised;
results.rpca_sparse_noise = rpca_sparse_noise;
results.info = info;
results.sweep_table = sweep_table;
results.freq_taper = freq_taper;
results.pad_factor = pad_factor;
results.metrics = struct('Q_noisy_full_dB',q_noisy_full,'Q_taper_full_dB',q_bp_full, ...
                         'Q_pca_full_dB',q_pca_full,'Q_rpca_full_dB',q_rpca_full, ...
                         'Q_noisy_tapered_dB',q_noisy_tap,'Q_taper_dB',q_bp_tap, ...
                         'Q_pca_tapered_dB',q_pca_tap,'Q_rpca_tapered_dB',q_rpca_tap, ...
                         'lambda_scale',lambda_scale,'rank_cap',rank_cap);
save('rpca_tapered_reproduction_results.mat','results','-v7.3');
fprintf('\nSaved results to rpca_tapered_reproduction_results.mat\n');

%% ========================================================================
%% Local functions
%% ========================================================================
function data = make_lowrank_fxy_synthetic(nt, nx, ny, dt)
    t = (0:nt-1)' * dt;
    [X,Y] = meshgrid(1:nx, 1:ny);
    X = X'; Y = Y';
    data = zeros(nt,nx,ny);

    events = [ ...
        0.16,  0.0035,  0.0010,  1.00, 25;
        0.34, -0.0020,  0.0018, -0.85, 22;
        0.58,  0.0015, -0.0025,  0.75, 18;
        0.72,  0.0030,  0.0005, -0.60, 20];

    x0 = (nx+1)/2; y0 = (ny+1)/2;
    for ie = 1:size(events,1)
        t0 = events(ie,1);
        sx = events(ie,2);
        sy = events(ie,3);
        amp = events(ie,4);
        fdom = events(ie,5);
        tau = t0 + sx*(X-x0) + sy*(Y-y0);
        for ix = 1:nx
            for iy = 1:ny
                data(:,ix,iy) = data(:,ix,iy) + amp * ricker_wavelet(t - tau(ix,iy), fdom);
            end
        end
    end
end

function w = ricker_wavelet(t, fdom)
    a = (pi*fdom*t).^2;
    w = (1 - 2*a) .* exp(-a);
end

function [noisy, gaussian_noise, erratic_noise] = add_gaussian_and_erratic_noise(clean, snr_db, erratic_fraction, erratic_amp)
    signal_power = mean(clean(:).^2);
    noise_power = signal_power / (10^(snr_db/10));
    gaussian_noise = sqrt(noise_power) * randn(size(clean));

    [nt,nx,ny] = size(clean);
    erratic_noise = zeros(size(clean));
    ntr = nx * ny;
    n_bad = max(1, round(erratic_fraction * ntr));
    bad_idx = randperm(ntr, n_bad);
    maxamp = max(abs(clean(:)));
    t = linspace(0,1,nt)';

    for k = 1:n_bad
        [ix,iy] = ind2sub([nx,ny], bad_idx(k));
        if rand < 0.7
            burst = randn(nt,1);
            burst = smooth_vector(burst, 5);
        else
            f0 = 12 + 18*rand;
            phase = 2*pi*rand;
            burst = cos(2*pi*f0*t + phase);
        end
        burst = burst ./ (max(abs(burst)) + eps);
        erratic_noise(:,ix,iy) = erratic_amp * maxamp * (0.7 + 0.6*rand) * burst;
    end

    noisy = clean + gaussian_noise + erratic_noise;
end

function y = smooth_vector(x, n)
    h = ones(n,1)/n;
    y = conv(x,h,'same');
end

function ntfft = choose_ntfft(nt, pad_factor)
    if nargin < 2 || isempty(pad_factor), pad_factor = 1; end
    ntfft = 2^nextpow2(max(nt, ceil(pad_factor*nt)));
end

function [freqs, w] = make_frequency_taper(ntfft, dt, taper)
% Raised-cosine taper for both positive and negative FFT frequencies.
    fs = 1/dt;
    fpos = (0:ntfft-1)'/(ntfft*dt);
    fnyq = fs/2;
    fabs = fpos;
    neg = fpos > fnyq;
    fabs(neg) = fs - fpos(neg);

    w = zeros(ntfft,1);
    f0 = taper.flow_stop;
    f1 = taper.flow_pass;
    f2 = taper.fhigh_pass;
    f3 = taper.fhigh_stop;

    w(fabs >= f1 & fabs <= f2) = 1;

    idx = fabs > f0 & fabs < f1;
    if f1 > f0
        x = (fabs(idx)-f0)/(f1-f0);
        w(idx) = 0.5 - 0.5*cos(pi*x);
    else
        w(fabs <= f1) = 1;
    end

    idx = fabs > f2 & fabs < f3;
    if f3 > f2
        x = (fabs(idx)-f2)/(f3-f2);
        w(idx) = 0.5 + 0.5*cos(pi*x);
    end

    freqs = fpos;
end

function Df = fft_time_padded(data, ntfft)
    [nt,nx,ny] = size(data);
    Dpad = zeros(ntfft,nx,ny);
    Dpad(1:nt,:,:) = data;
    Df = fft(Dpad, [], 1);
end

function data = ifft_time_crop(Df, nt)
    dpad = real(ifft(Df, [], 1));
    data = dpad(1:nt,:,:);
end

function den = fxy_taper_reconstruct(data, dt, taper, pad_factor)
    [nt,~,~] = size(data);
    ntfft = choose_ntfft(nt, pad_factor);
    Df = fft_time_padded(data, ntfft);
    [~,w] = make_frequency_taper(ntfft, dt, taper);
    den_f = bsxfun(@times, Df, reshape(w,[],1,1));
    den = ifft_time_crop(den_f, nt);
end

function den = fxy_rank_filter_tapered(data, dt, taper, pad_factor, rank_pca)
    [nt,~,~] = size(data);
    ntfft = choose_ntfft(nt, pad_factor);
    Df = fft_time_padded(data, ntfft);
    den_f = zeros(size(Df));
    [~,w] = make_frequency_taper(ntfft, dt, taper);
    npos = floor(ntfft/2)+1;
    pos_bins = find((1:ntfft)' <= npos & w > taper.min_weight_to_process);

    for ii = 1:numel(pos_bins)
        k = pos_bins(ii);
        M = squeeze(Df(k,:,:));
        [U,S,V] = svd(M,'econ');
        r = min([rank_pca, size(U,2), size(V,2)]);
        L = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';
        L = w(k) * L;
        den_f(k,:,:) = L;
        kc = ntfft - k + 2;
        if k > 1 && k < npos && kc <= ntfft
            den_f(kc,:,:) = conj(L);
        end
    end
    den = ifft_time_crop(den_f, nt);
end

function [den, sparse_noise_time, info] = fxy_rpca_pcp_denoise_tapered(data, dt, taper, pad_factor, opts)
    [nt,~,~] = size(data);
    ntfft = choose_ntfft(nt, pad_factor);
    Df = fft_time_padded(data, ntfft);
    den_f = zeros(size(Df));
    sparse_f = zeros(size(Df));
    [freqs,w] = make_frequency_taper(ntfft, dt, taper);
    npos = floor(ntfft/2)+1;
    pos_bins = find((1:ntfft)' <= npos & w > taper.min_weight_to_process);

    nfreq = numel(pos_bins);
    ranks = zeros(nfreq,1);
    ranks_after = zeros(nfreq,1);
    iters = zeros(nfreq,1);
    relres = zeros(nfreq,1);

    for ii = 1:nfreq
        k = pos_bins(ii);
        M = squeeze(Df(k,:,:));
        [L,S,rpca_info] = rpca_pcp_ialm(M, opts);
        ranks(ii) = rank_estimate(L);

        if isfield(opts,'rank_cap') && ~isempty(opts.rank_cap)
            L = rank_project(L, opts.rank_cap);
        end
        ranks_after(ii) = rank_estimate(L);

        Lw = w(k) * L;
        Sw = w(k) * S;
        den_f(k,:,:) = Lw;
        sparse_f(k,:,:) = Sw;
        kc = ntfft - k + 2;
        if k > 1 && k < npos && kc <= ntfft
            den_f(kc,:,:) = conj(Lw);
            sparse_f(kc,:,:) = conj(Sw);
        end
        iters(ii) = rpca_info.iter;
        relres(ii) = rpca_info.relres;
        if isfield(opts,'verbose') && opts.verbose
            fprintf('  f=%6.2f Hz, w=%.2f, iter=%4d, rank=%2d -> %2d, relres=%.2e\n', ...
                freqs(k), w(k), iters(ii), ranks(ii), ranks_after(ii), relres(ii));
        end
    end

    den = ifft_time_crop(den_f, nt);
    sparse_noise_time = ifft_time_crop(sparse_f, nt);
    info = struct('freq_bins',pos_bins,'frequencies',freqs(pos_bins), ...
                  'weights',w(pos_bins),'nfreq_processed',nfreq, ...
                  'ranks',ranks,'ranks_after',ranks_after, ...
                  'iters',iters,'relres',relres);
end

function [L,S,info] = rpca_pcp_ialm(D, opts)
% Inexact Augmented Lagrange Multiplier solver for PCP:
%   min ||L||_* + lambda ||S||_1, subject to D = L + S.
% Works for real or complex matrices.
    if nargin < 2, opts = struct(); end
    [m,n] = size(D);

    if isfield(opts,'lambda') && ~isempty(opts.lambda)
        lambda = opts.lambda;
    elseif isfield(opts,'lambda_scale') && ~isempty(opts.lambda_scale)
        lambda = opts.lambda_scale / sqrt(max(m,n));
    else
        lambda = 1 / sqrt(max(m,n));
    end
    if ~isfield(opts,'maxIter') || isempty(opts.maxIter), maxIter = 1000; else, maxIter = opts.maxIter; end
    if ~isfield(opts,'tol') || isempty(opts.tol), tol = 1e-7; else, tol = opts.tol; end

    normD = norm(D,'fro') + eps;
    norm2D = norm(D,2);
    normInfD = max(abs(D(:))) / lambda;
    dual_norm = max(norm2D, normInfD) + eps;
    Y = D / dual_norm;

    mu = 1.25 / (norm2D + eps);
    mu_bar = mu * 1e7;
    rho = 1.5;

    L = zeros(m,n,class(D));
    S = zeros(m,n,class(D));
    relres = Inf;

    for iter = 1:maxIter
        L = svd_threshold(D - S + (1/mu)*Y, 1/mu);
        S = soft_threshold_complex(D - L + (1/mu)*Y, lambda/mu);
        Z = D - L - S;
        relres = norm(Z,'fro') / normD;
        if relres < tol
            break;
        end
        Y = Y + mu * Z;
        mu = min(mu*rho, mu_bar);
    end

    info.iter = iter;
    info.relres = relres;
    info.lambda = lambda;
    info.rank = rank_estimate(L);
end

function X = soft_threshold_complex(Y, tau)
    mag = abs(Y);
    X = max(mag - tau, 0) ./ (mag + eps) .* Y;
end

function L = svd_threshold(M, tau)
    [U,S,V] = svd(M,'econ');
    s = diag(S);
    s2 = max(s - tau, 0);
    r = sum(s2 > 0);
    if r == 0
        L = zeros(size(M), class(M));
    else
        L = U(:,1:r) * diag(s2(1:r)) * V(:,1:r)';
    end
end

function Lr = rank_project(M, r)
    [U,S,V] = svd(M,'econ');
    r = min([r, size(U,2), size(V,2)]);
    if r <= 0
        Lr = zeros(size(M), class(M));
    else
        Lr = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';
    end
end

function r = rank_estimate(M)
    s = svd(M,'econ');
    if isempty(s) || s(1) == 0
        r = 0;
    else
        r = sum(s > max(size(M))*eps(s(1)));
    end
end

function q = quality_db(clean_ref, estimate)
    err = estimate - clean_ref;
    q = 10 * log10( sum(clean_ref(:).^2) / (sum(err(:).^2) + eps) );
end

function s = rankcap_to_string(r)
    if isempty(r)
        s = 'none';
    else
        s = sprintf('%d', r);
    end
end

function wiggle_plot_fixed(d, t, amp_ref)
    [nt,nx] = size(d);
    if nargin < 2 || isempty(t), t = (0:nt-1)'; end
    if nargin < 3 || isempty(amp_ref), amp_ref = max(abs(d(:))) + eps; end
    d = d ./ (amp_ref + eps);
    scale = 0.45;
    hold on;
    for ix = 1:nx
        tr = d(:,ix);
        x = ix + scale * tr;
        plot(x, t, 'k', 'LineWidth', 0.6);
        pos = tr;
        pos(pos < 0) = 0;
        xp = ix + scale * pos;
        fill([ix*ones(nt,1); flipud(xp)], [t(:); flipud(t(:))], 'k', ...
             'EdgeColor','none', 'FaceAlpha',0.25);
    end
    hold off;
    set(gca,'YDir','reverse','Box','on','FontName','Arial','FontSize',9);
    xlabel('trace number'); ylabel('time (s)');
    xlim([0 nx+1]); ylim([min(t) max(t)]);
end
