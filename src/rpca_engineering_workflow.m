%% rpca_engineering_workflow.m
% Engineering-style single-file MATLAB implementation of f-x-y RPCA/PCP
% seismic erratic-noise attenuation inspired by Cheng, Chen & Sacchi (2015).
%
% Engineering workflow additions:
%   1) Sliding time-window processing with overlap-add.
%   2) Adaptive frequency-dependent rank projection.
%   3) Adaptive frequency-dependent lambda for PCP/RPCA.
%   4) 2D/3D .mat input support and optional SEG-Y I/O hooks.
%   5) QC outputs: denoised data, removed noise, sparse noise, ranks,
%      iteration counts, Q metrics when a clean reference is available.
%
% Recommended first run:
%   run('src/rpca_engineering_workflow.m')
%
% For your own .mat data, set USER_DATA_MAT and DATA_VAR_NAME below.
% The data variable should be [nt,nx] or [nt,nx,ny].
%
% Author: generated for Dawei Liu
% Date: 2026-06-05

clear; close all; clc;

%% ========================= User parameters ==============================
% ---------------- Input mode ----------------
% 'synthetic' : generate clean/noisy synthetic example
% 'mat'       : read seismic matrix from USER_DATA_MAT
% 'segy'      : optional; requires MATLAB segyread/segywrite support
INPUT_MODE = 'synthetic';

USER_DATA_MAT = '';       % e.g., 'my_seismic_patch.mat'
DATA_VAR_NAME = 'data';   % variable name in mat file

USER_SEGY_IN  = '';       % optional SEG-Y input path if INPUT_MODE='segy'
USER_SEGY_OUT = 'rpca_denoised_output.sgy';

% ---------------- Sampling ----------------
dt = 0.002;               % seconds. Used for synthetic/mat input.

% ---------------- Smooth spectral taper ----------------
freq_taper.flow_stop  = 0;     % Hz
freq_taper.flow_pass  = 2;     % Hz
freq_taper.fhigh_pass = 45;    % Hz
freq_taper.fhigh_stop = 75;    % Hz
freq_taper.min_weight_to_process = 0.02;
pad_factor = 2;

% ---------------- Sliding time-window ----------------
USE_SLIDING_WINDOW = true;
win_length_s = 0.50;      % seconds. Try 0.5-1.0 s for field data.
win_overlap  = 0.50;      % 0.5 means 50 percent overlap.

% ---------------- Baseline PCA/eigenimage filtering ----------------
rank_pca = 3;

% ---------------- RPCA/PCP settings ----------------
opts = struct();
opts.maxIter = 800;
opts.tol = 1e-6;
opts.verbose = false;

% lambda(f) = lambda_scale/sqrt(max(nx,ny)) * (1 + lambda_freq_boost*f/fref)
% Smaller lambda sends more energy to sparse S. Larger lambda keeps more in L.
opts.lambda_scale = 0.70;
opts.lambda_adaptive = true;
opts.lambda_freq_boost = 0.40;     % 0 disables frequency dependence
opts.lambda_fref = freq_taper.fhigh_pass;

% Rank projection after RPCA. This is not pure PCP, but improves seismic QC.
% If rank_mode='fixed', rank_fixed is used for all frequencies.
% If rank_mode='adaptive', rank decreases with frequency.
% Set opts.use_rank_projection=false to disable.
opts.use_rank_projection = true;
opts.rank_mode = 'adaptive';       % 'adaptive' or 'fixed'
opts.rank_fixed = 4;
opts.rank_min = 2;
opts.rank_max = 7;
opts.rank_decay_hz = 35;           % larger = slower rank decay with f
opts.rank_energy_keep = [];        % e.g., 0.995. Empty = use rank formula only.

% ---------------- Optional synthetic parameter sweep ----------------
% Only used when INPUT_MODE='synthetic'. For field data, keep false.
DO_SYNTHETIC_SWEEP = true;
lambda_grid = [0.45 0.60 0.75 0.90];
rankmax_grid = [5 7 9];
rankmin_grid = [1 2 3];

% ---------------- Synthetic controls ----------------
snr_gaussian_db  = 12;
erratic_fraction = 0.08;
erratic_amp      = 3.0;
rng_seed         = 2025;

% ---------------- Display / save ----------------
DISPLAY_IY = [];          % [] = middle y slice
SAVE_MAT_RESULTS = true;
RESULT_MAT_NAME = 'rpca_engineering_workflow_results.mat';

%% ========================= Load or synthesize data =======================
rng(rng_seed);
[clean, noisy, input_meta, has_clean] = load_or_make_data(INPUT_MODE, USER_DATA_MAT, DATA_VAR_NAME, USER_SEGY_IN, dt, ...
    snr_gaussian_db, erratic_fraction, erratic_amp, rng_seed);

if isfield(input_meta,'dt') && ~isempty(input_meta.dt)
    dt = input_meta.dt;
end

[nt,nx,ny] = size(noisy);
fprintf('Input data size: nt=%d, nx=%d, ny=%d, dt=%.6f s\n', nt, nx, ny, dt);

%% ========================= Tapered clean reference =======================
if has_clean
    clean_tapered = process_by_windows(clean, dt, freq_taper, pad_factor, USE_SLIDING_WINDOW, win_length_s, win_overlap, 'taper', struct());
else
    clean_tapered = [];
end

%% ========================= Main processing ===============================
fprintf('\nRunning v5 processing...\n');
fprintf('  Sliding window: %d, window=%.3f s, overlap=%.0f%%\n', USE_SLIDING_WINDOW, win_length_s, 100*win_overlap);
fprintf('  Adaptive lambda: %d, adaptive rank: %s\n', opts.lambda_adaptive, opts.rank_mode);

fprintf('\n1) Smooth tapered reconstruction only...\n');
bandpass_only = process_by_windows(noisy, dt, freq_taper, pad_factor, USE_SLIDING_WINDOW, win_length_s, win_overlap, 'taper', struct());

fprintf('2) f-x-y PCA/eigenimage filtering, rank=%d...\n', rank_pca);
pca_opts = struct('rank_pca',rank_pca);
pca_denoised = process_by_windows(noisy, dt, freq_taper, pad_factor, USE_SLIDING_WINDOW, win_length_s, win_overlap, 'pca', pca_opts);

fprintf('3) Enhanced f-x-y RPCA/PCP denoising...\n');
[rpca_denoised, rpca_sparse_noise, info] = process_by_windows(noisy, dt, freq_taper, pad_factor, USE_SLIDING_WINDOW, win_length_s, win_overlap, 'rpca', opts);
removed_noise = noisy - rpca_denoised;

%% ========================= Optional synthetic sweep ======================
sweep_table = [];
if has_clean && DO_SYNTHETIC_SWEEP
    fprintf('\nRunning v5 synthetic sweep against tapered clean reference...\n');
    best.Q = -Inf;
    best.opts = opts;
    best.den = rpca_denoised;
    best.sparse = rpca_sparse_noise;
    best.info = info;
    row = 0;
    for il = 1:numel(lambda_grid)
        for imax = 1:numel(rankmax_grid)
            for imin = 1:numel(rankmin_grid)
                tmp_opts = opts;
                tmp_opts.lambda_scale = lambda_grid(il);
                tmp_opts.rank_max = rankmax_grid(imax);
                tmp_opts.rank_min = rankmin_grid(imin);
                if tmp_opts.rank_min > tmp_opts.rank_max
                    continue;
                end
                [tmp_den,tmp_sparse,tmp_info] = process_by_windows(noisy, dt, freq_taper, pad_factor, USE_SLIDING_WINDOW, win_length_s, win_overlap, 'rpca', tmp_opts);
                q = quality_db(clean_tapered, tmp_den);
                row = row + 1;
                sweep_table(row,:) = [tmp_opts.lambda_scale, tmp_opts.rank_min, tmp_opts.rank_max, q, mean_safe(tmp_info.ranks_before), mean_safe(tmp_info.ranks_after)]; %#ok<SAGROW>
                fprintf('  lambda=%.2f, rank_min=%d, rank_max=%d, Q_tapered=%6.2f dB, rank %.2f -> %.2f\n', ...
                    tmp_opts.lambda_scale, tmp_opts.rank_min, tmp_opts.rank_max, q, mean_safe(tmp_info.ranks_before), mean_safe(tmp_info.ranks_after));
                if q > best.Q
                    best.Q = q;
                    best.opts = tmp_opts;
                    best.den = tmp_den;
                    best.sparse = tmp_sparse;
                    best.info = tmp_info;
                end
            end
        end
    end
    opts = best.opts;
    rpca_denoised = best.den;
    rpca_sparse_noise = best.sparse;
    info = best.info;
    removed_noise = noisy - rpca_denoised;
    fprintf('Best v5 setting: lambda=%.2f, rank_min=%d, rank_max=%d, Q_tapered=%.2f dB\n', ...
        opts.lambda_scale, opts.rank_min, opts.rank_max, best.Q);
end

%% ========================= Metrics =======================================
metrics = struct();
if has_clean
    metrics.Q_noisy_full_dB = quality_db(clean, noisy);
    metrics.Q_taper_full_dB = quality_db(clean, bandpass_only);
    metrics.Q_pca_full_dB   = quality_db(clean, pca_denoised);
    metrics.Q_rpca_full_dB  = quality_db(clean, rpca_denoised);

    metrics.Q_noisy_tapered_dB = quality_db(clean_tapered, noisy);
    metrics.Q_taper_dB         = quality_db(clean_tapered, bandpass_only);
    metrics.Q_pca_tapered_dB   = quality_db(clean_tapered, pca_denoised);
    metrics.Q_rpca_tapered_dB  = quality_db(clean_tapered, rpca_denoised);

    fprintf('\nQ_full = 10log10(||clean||^2 / ||estimate-clean||^2)\n');
    fprintf('Noisy data       : %.2f dB\n', metrics.Q_noisy_full_dB);
    fprintf('Taper only       : %.2f dB\n', metrics.Q_taper_full_dB);
    fprintf('Eigenimage PCA   : %.2f dB\n', metrics.Q_pca_full_dB);
    fprintf('Enhanced RPCA    : %.2f dB\n', metrics.Q_rpca_full_dB);

    fprintf('\nQ_tapered_ref = 10log10(||clean_tapered||^2 / ||estimate-clean_tapered||^2)\n');
    fprintf('Noisy data       : %.2f dB\n', metrics.Q_noisy_tapered_dB);
    fprintf('Taper only       : %.2f dB\n', metrics.Q_taper_dB);
    fprintf('Eigenimage PCA   : %.2f dB\n', metrics.Q_pca_tapered_dB);
    fprintf('Enhanced RPCA    : %.2f dB\n', metrics.Q_rpca_tapered_dB);
else
    fprintf('\nNo clean reference was provided. Q metrics are skipped.\n');
end

fprintf('\nRPCA QC summary:\n');
fprintf('  Processed windows: %d\n', info.nwin);
fprintf('  Processed frequency slices total: %d\n', numel(info.ranks_before));
fprintf('  Mean rank before projection: %.2f\n', mean_safe(info.ranks_before));
fprintf('  Mean rank after  projection: %.2f\n', mean_safe(info.ranks_after));
fprintf('  Mean IALM iterations: %.1f\n', mean_safe(info.iters));
fprintf('  Mean final relative residual: %.2e\n', mean_safe(info.relres));
fprintf('  Final lambda_scale=%.2f, rank_mode=%s, rank_min=%d, rank_max=%d\n', ...
    opts.lambda_scale, opts.rank_mode, opts.rank_min, opts.rank_max);

%% ========================= Display =======================================
if isempty(DISPLAY_IY)
    iy = max(1, round(ny/2));
else
    iy = min(max(1, DISPLAY_IY), ny);
end
t = (0:nt-1)*dt;
plot_scale = max(abs(noisy(:,:,iy)), [], 'all') + eps;

figure('Name','v5 engineering RPCA: wiggle QC','Color','w','Position',[60 60 1400 800]);
if has_clean
    subplot(2,3,1); wiggle_plot_fixed(clean(:,:,iy), t, plot_scale); title('(a) Clean');
else
    subplot(2,3,1); wiggle_plot_fixed(noisy(:,:,iy), t, plot_scale); title('(a) Input noisy');
end
subplot(2,3,2); wiggle_plot_fixed(noisy(:,:,iy), t, plot_scale); title('(b) Noisy');
subplot(2,3,3); wiggle_plot_fixed(bandpass_only(:,:,iy), t, plot_scale); title('(c) Smooth taper only');
subplot(2,3,4); wiggle_plot_fixed(pca_denoised(:,:,iy), t, plot_scale); title('(d) f-x-y PCA');
subplot(2,3,5); wiggle_plot_fixed(rpca_denoised(:,:,iy), t, plot_scale); title('(e) v5 enhanced RPCA');
subplot(2,3,6); wiggle_plot_fixed(removed_noise(:,:,iy), t, plot_scale); title('(f) Removed noise');

figure('Name','v5 engineering RPCA: image QC','Color','w','Position',[120 120 1450 780]);
clim = [-1 1] * max(abs(noisy(:,:,iy)), [], 'all');
if has_clean
    subplot(2,3,1); imagesc(1:nx,t,clean(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('Clean'); colorbar;
else
    subplot(2,3,1); imagesc(1:nx,t,noisy(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('Input noisy'); colorbar;
end
subplot(2,3,2); imagesc(1:nx,t,noisy(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('Noisy'); colorbar;
subplot(2,3,3); imagesc(1:nx,t,bandpass_only(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('Smooth taper'); colorbar;
subplot(2,3,4); imagesc(1:nx,t,pca_denoised(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('PCA'); colorbar;
subplot(2,3,5); imagesc(1:nx,t,rpca_denoised(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('v5 RPCA'); colorbar;
subplot(2,3,6); imagesc(1:nx,t,removed_noise(:,:,iy),clim); axis tight; set(gca,'YDir','reverse'); title('Removed noise'); colorbar;

if has_clean
    figure('Name','v5 RPCA error QC','Color','w','Position',[200 200 1000 420]);
    err = rpca_denoised(:,:,iy) - clean_tapered(:,:,iy);
    err_clim = [-1 1] * max(abs(noisy(:,:,iy)), [], 'all');
    subplot(1,2,1); imagesc(1:nx,t,err,err_clim); axis tight; set(gca,'YDir','reverse'); title('RPCA error vs tapered clean'); colorbar;
    subplot(1,2,2); wiggle_plot_fixed(err, t, plot_scale); title('Error wiggle, common scale');
end

figure('Name','v5 spectral taper and adaptive rank/lambda','Color','w','Position',[240 240 1100 420]);
ntfft_plot = choose_ntfft(min(nt, max(16, round(win_length_s/dt))), pad_factor);
[freqs_plot,w_plot] = make_frequency_taper(ntfft_plot, dt, freq_taper);
pos = 1:floor(ntfft_plot/2)+1;
subplot(1,3,1); plot(freqs_plot(pos), w_plot(pos), 'k', 'LineWidth', 1.5); grid on; xlabel('Frequency (Hz)'); ylabel('Weight'); title('Raised-cosine taper'); ylim([-0.05 1.05]);
fpos = freqs_plot(pos);
rank_curve = arrayfun(@(f) choose_rank_cap(f, opts), fpos);
lambda_curve = arrayfun(@(f) choose_lambda(f, nx, ny, opts), fpos);
subplot(1,3,2); plot(fpos, rank_curve, 'k', 'LineWidth', 1.5); grid on; xlabel('Frequency (Hz)'); ylabel('Rank cap'); title('Adaptive rank');
subplot(1,3,3); plot(fpos, lambda_curve, 'k', 'LineWidth', 1.5); grid on; xlabel('Frequency (Hz)'); ylabel('\lambda'); title('Adaptive RPCA lambda');

%% ========================= Save results ==================================
if SAVE_MAT_RESULTS
    results = struct();
    results.clean = clean;
    results.clean_tapered = clean_tapered;
    results.noisy = noisy;
    results.bandpass_only = bandpass_only;
    results.pca_denoised = pca_denoised;
    results.rpca_denoised = rpca_denoised;
    results.rpca_sparse_noise = rpca_sparse_noise;
    results.removed_noise = removed_noise;
    results.metrics = metrics;
    results.info = info;
    results.opts = opts;
    results.freq_taper = freq_taper;
    results.pad_factor = pad_factor;
    results.sliding_window = struct('enabled',USE_SLIDING_WINDOW,'win_length_s',win_length_s,'win_overlap',win_overlap);
    results.sweep_table = sweep_table;
    results.input_meta = input_meta;
    save(RESULT_MAT_NAME, 'results', '-v7.3');
    fprintf('\nSaved results to %s\n', RESULT_MAT_NAME);
end

%% ========================= Optional SEG-Y output ==========================
if strcmpi(INPUT_MODE,'segy') && ~isempty(USER_SEGY_OUT)
    try_write_segy(USER_SEGY_IN, USER_SEGY_OUT, rpca_denoised, input_meta); %#ok<UNRCH> Reached when INPUT_MODE is set to 'segy'.
end

%% ========================================================================
%% Local functions
%% ========================================================================
function [clean, noisy, meta, has_clean] = load_or_make_data(mode, matfile, varname, segyin, dt, snr_db, err_frac, err_amp, seed)
    meta = struct('mode',mode,'dt',dt);
    switch lower(mode)
        case 'synthetic'
            rng(seed);
            nt = round(1.0/dt)+1; nx = 30; ny = 30;
            clean = make_lowrank_fxy_synthetic(nt,nx,ny,dt);
            clean = clean ./ (max(abs(clean(:))) + eps);
            [noisy, gaussian_noise, erratic_noise] = add_gaussian_and_erratic_noise(clean, snr_db, err_frac, err_amp);
            meta.gaussian_noise = gaussian_noise;
            meta.erratic_noise = erratic_noise;
            meta.dt = dt;
            has_clean = true;
            fprintf('Created synthetic 3D data.\n');
        case 'mat'
            if isempty(matfile), error('USER_DATA_MAT is empty.'); end
            tmp = load(matfile);
            if ~isfield(tmp,varname), error('Variable "%s" not found in %s.', varname, matfile); end
            noisy = double(tmp.(varname));
            if ismatrix(noisy), noisy = reshape(noisy,size(noisy,1),size(noisy,2),1); end
            noisy = noisy ./ (max(abs(noisy(:))) + eps);
            clean = [];
            if isfield(tmp,'clean')
                clean = double(tmp.clean);
                if ismatrix(clean), clean = reshape(clean,size(clean,1),size(clean,2),1); end
                clean = clean ./ (max(abs(clean(:))) + eps);
                has_clean = true;
            else
                has_clean = false;
            end
            if isfield(tmp,'dt'), meta.dt = tmp.dt; end
            meta.matfile = matfile;
            fprintf('Loaded MAT data: %s, variable=%s.\n', matfile, varname);
        case 'segy'
            [noisy, meta] = try_read_segy(segyin, dt);
            clean = [];
            has_clean = false;
        otherwise
            error('Unknown INPUT_MODE: %s', mode);
    end
end

function [data, meta] = try_read_segy(segyin, default_dt)
    if isempty(segyin), error('USER_SEGY_IN is empty.'); end
    meta = struct('mode','segy','dt',default_dt,'segy_in',segyin);
    if exist('segyread','file') ~= 2
        error(['MATLAB segyread was not found. Use INPUT_MODE=''mat'', or install a SEG-Y reader. ', ...
               'You can load SEG-Y with your own reader and save data as [nt,nx] or [nt,nx,ny] MAT.']);
    end
    S = segyread(segyin);
    if istable(S)
        % Some MATLAB versions return table-like trace data.
        vars = S.Properties.VariableNames;
        candidate = '';
        for i = 1:numel(vars)
            if isnumeric(S.(vars{i})) && size(S.(vars{i}),1) > 10
                candidate = vars{i}; break;
            end
        end
        if isempty(candidate), error('Could not identify trace samples from segyread output.'); end
        data = double(S.(candidate));
    elseif isnumeric(S)
        data = double(S);
    else
        error('Unsupported segyread output. Please convert SEG-Y to MAT first.');
    end
    if size(data,1) < size(data,2)
        % Usually nt should be rows. Leave as-is if user already knows otherwise.
    end
    if ismatrix(data), data = reshape(data,size(data,1),size(data,2),1); end
    data = data ./ (max(abs(data(:))) + eps);
end

function try_write_segy(~, segyout, denoised, ~)
    if exist('segywrite','file') ~= 2
        fprintf('segywrite was not found. Skipping SEG-Y output. MAT result has been saved.\n');
        return;
    end
    try
        % This is intentionally conservative because MATLAB SEG-Y APIs vary.
        % Users with a known SEG-Y toolbox should replace this block with their
        % standard header-preserving writer.
        data2d = reshape(denoised, size(denoised,1), []);
        segywrite(data2d, segyout);
        fprintf('Wrote denoised SEG-Y to %s. Header preservation depends on your segywrite implementation.\n', segyout);
    catch ME
        fprintf('Could not write SEG-Y output: %s\n', ME.message);
        fprintf('Denoised data are still available in the saved MAT results.\n');
    end
end

function [out, sparse_out, info_all] = process_by_windows(data, dt, taper, pad_factor, use_window, win_length_s, win_overlap, method, opts)
    [nt,nx,ny] = size(data);
    if nargin < 6 || isempty(use_window), use_window = false; end
    if ~use_window
        starts = 1; lens = nt;
    else
        win_len = max(16, round(win_length_s/dt));
        win_len = min(win_len, nt);
        hop = max(1, round(win_len*(1-win_overlap)));
        starts = 1:hop:(nt-win_len+1);
        if starts(end) ~= nt-win_len+1
            starts = [starts, nt-win_len+1];
        end
        lens = win_len;
    end

    out_acc = zeros(nt,nx,ny);
    sp_acc  = zeros(nt,nx,ny);
    w_acc   = zeros(nt,1);
    info_all = init_info();

    for iw = 1:numel(starts)
        i1 = starts(iw);
        if isscalar(lens), Lw = lens; else, Lw = lens(iw); end
        i2 = i1 + Lw - 1;
        xw = data(i1:i2,:,:);
        tw = local_sine_taper(Lw);
        xw_tapered = bsxfun(@times, xw, reshape(tw,[],1,1));

        switch lower(method)
            case 'taper'
                yw = fxy_taper_reconstruct(xw_tapered, dt, taper, pad_factor);
                sw = zeros(size(yw));
                infow = init_info();
            case 'pca'
                yw = fxy_rank_filter_tapered(xw_tapered, dt, taper, pad_factor, opts.rank_pca);
                sw = zeros(size(yw));
                infow = init_info();
            case 'rpca'
                [yw, sw, infow] = fxy_rpca_pcp_denoise_tapered(xw_tapered, dt, taper, pad_factor, opts);
            otherwise
                error('Unknown processing method: %s', method);
        end

        out_acc(i1:i2,:,:) = out_acc(i1:i2,:,:) + bsxfun(@times, yw, reshape(tw,[],1,1));
        sp_acc(i1:i2,:,:)  = sp_acc(i1:i2,:,:)  + bsxfun(@times, sw, reshape(tw,[],1,1));
        w_acc(i1:i2) = w_acc(i1:i2) + tw.^2;
        info_all = merge_info(info_all, infow, iw, i1, i2);
        if strcmpi(method,'rpca')
            fprintf('  window %d/%d, samples %d:%d, mean rank %.2f -> %.2f\n', ...
                iw, numel(starts), i1, i2, mean_safe(infow.ranks_before), mean_safe(infow.ranks_after));
        end
    end
    out = bsxfun(@rdivide, out_acc, reshape(w_acc+eps,[],1,1));
    sparse_out = bsxfun(@rdivide, sp_acc, reshape(w_acc+eps,[],1,1));
    info_all.nwin = numel(starts);
    info_all.window_starts = starts(:);
end

function info = init_info()
    info = struct('freq_bins',[],'frequencies',[],'weights',[], ...
        'ranks_before',[],'ranks_after',[],'rank_caps',[], ...
        'lambdas',[],'iters',[],'relres',[], ...
        'nwin',0,'window_index',[],'window_start',[],'window_end',[]);
end

function info = merge_info(info, infow, iw, i1, i2)
    fields = {'freq_bins','frequencies','weights','ranks_before','ranks_after','rank_caps','lambdas','iters','relres'};
    n = numel(infow.ranks_before);
    for k = 1:numel(fields)
        f = fields{k};
        if isfield(infow,f)
            info.(f) = [info.(f); infow.(f)(:)];
        end
    end
    info.window_index = [info.window_index; iw*ones(n,1)];
    info.window_start = [info.window_start; i1*ones(n,1)];
    info.window_end   = [info.window_end; i2*ones(n,1)];
end

function w = local_sine_taper(n)
    if n <= 1
        w = 1;
    else
        x = (0:n-1)'/(n-1);
        w = sin(pi*x);
        % Avoid exactly zero edges causing weak boundary normalization.
        w = max(w, 0.05);
    end
end

function den = fxy_taper_reconstruct(data, dt, taper, pad_factor)
    [nt,~,~] = size(data);
    ntfft = choose_ntfft(nt,pad_factor);
    Df = fft_time_padded(data,ntfft);
    [~,w] = make_frequency_taper(ntfft,dt,taper);
    den_f = bsxfun(@times,Df,reshape(w,[],1,1));
    den = ifft_time_crop(den_f,nt);
end

function den = fxy_rank_filter_tapered(data, dt, taper, pad_factor, rank_pca)
    [nt,~,~] = size(data);
    ntfft = choose_ntfft(nt,pad_factor);
    Df = fft_time_padded(data,ntfft);
    den_f = zeros(size(Df));
    [~,w] = make_frequency_taper(ntfft,dt,taper);
    npos = floor(ntfft/2)+1;
    pos_bins = find((1:ntfft)'<=npos & w>taper.min_weight_to_process);
    for ii = 1:numel(pos_bins)
        k = pos_bins(ii);
        M = squeeze(Df(k,:,:));
        [U,S,V] = svd(M,'econ');
        r = min([rank_pca,size(U,2),size(V,2)]);
        L = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
        Lw = w(k)*L;
        den_f(k,:,:) = Lw;
        kc = ntfft-k+2;
        if k>1 && k<npos && kc<=ntfft, den_f(kc,:,:) = conj(Lw); end
    end
    den = ifft_time_crop(den_f,nt);
end

function [den, sparse_noise_time, info] = fxy_rpca_pcp_denoise_tapered(data, dt, taper, pad_factor, opts)
    [nt,nx,ny] = size(data);
    ntfft = choose_ntfft(nt,pad_factor);
    Df = fft_time_padded(data,ntfft);
    den_f = zeros(size(Df)); sparse_f = zeros(size(Df));
    [freqs,w] = make_frequency_taper(ntfft,dt,taper);
    npos = floor(ntfft/2)+1;
    pos_bins = find((1:ntfft)'<=npos & w>taper.min_weight_to_process);
    nfreq = numel(pos_bins);

    ranks_before = zeros(nfreq,1); ranks_after = zeros(nfreq,1); rank_caps = zeros(nfreq,1);
    lambdas = zeros(nfreq,1); iters = zeros(nfreq,1); relres = zeros(nfreq,1);

    for ii = 1:nfreq
        k = pos_bins(ii);
        f = freqs(k);
        M = squeeze(Df(k,:,:));
        opts_f = opts;
        opts_f.lambda = choose_lambda(f,nx,ny,opts);
        lambdas(ii) = opts_f.lambda;
        [L,S,rpca_info] = rpca_pcp_ialm(M, opts_f);
        ranks_before(ii) = rank_estimate(L);

        rc = choose_rank_cap(f, opts);
        if isfield(opts,'use_rank_projection') && opts.use_rank_projection && ~isempty(rc)
            if isfield(opts,'rank_energy_keep') && ~isempty(opts.rank_energy_keep)
                L = rank_project_energy(L, opts.rank_energy_keep, rc);
            else
                L = rank_project(L, rc);
            end
        end
        ranks_after(ii) = rank_estimate(L);
        rank_caps(ii) = rc;

        Lw = w(k)*L; Sw = w(k)*S;
        den_f(k,:,:) = Lw; sparse_f(k,:,:) = Sw;
        kc = ntfft-k+2;
        if k>1 && k<npos && kc<=ntfft
            den_f(kc,:,:) = conj(Lw); sparse_f(kc,:,:) = conj(Sw);
        end
        iters(ii) = rpca_info.iter;
        relres(ii) = rpca_info.relres;
        if isfield(opts,'verbose') && opts.verbose
            fprintf('    f=%6.2f Hz, lambda=%.4g, rank=%d->%d cap=%d, iter=%d\n', ...
                f, opts_f.lambda, ranks_before(ii), ranks_after(ii), rc, iters(ii));
        end
    end
    den = ifft_time_crop(den_f,nt);
    sparse_noise_time = ifft_time_crop(sparse_f,nt);
    info = struct('freq_bins',pos_bins(:),'frequencies',freqs(pos_bins(:)), ...
        'weights',w(pos_bins(:)),'ranks_before',ranks_before,'ranks_after',ranks_after, ...
        'rank_caps',rank_caps,'lambdas',lambdas,'iters',iters,'relres',relres);
end

function lambda = choose_lambda(f, nx, ny, opts)
    base = opts.lambda_scale / sqrt(max(nx,ny));
    if isfield(opts,'lambda_adaptive') && opts.lambda_adaptive
        if isfield(opts,'lambda_fref') && opts.lambda_fref>0, fref = opts.lambda_fref; else, fref = max(f,1); end
        if isfield(opts,'lambda_freq_boost'), boost = opts.lambda_freq_boost; else, boost = 0; end
        lambda = base * (1 + boost * min(f/fref, 2));
    else
        lambda = base;
    end
end

function rc = choose_rank_cap(f, opts)
    if ~isfield(opts,'use_rank_projection') || ~opts.use_rank_projection
        rc = [];
        return;
    end
    if ~isfield(opts,'rank_mode'), opts.rank_mode = 'fixed'; end
    switch lower(opts.rank_mode)
        case 'adaptive'
            rmin = opts.rank_min; rmax = opts.rank_max;
            fd = max(opts.rank_decay_hz, eps);
            rc = round(rmin + (rmax-rmin)*exp(-f/fd));
            rc = max(rmin, min(rmax, rc));
        otherwise
            rc = opts.rank_fixed;
    end
end

function [L,S,info] = rpca_pcp_ialm(D, opts)
% Inexact ALM solver for Principal Component Pursuit:
% min ||L||_* + lambda ||S||_1 subject to D = L + S.
    if nargin<2, opts=struct(); end
    [m,n] = size(D);
    if isfield(opts,'lambda') && ~isempty(opts.lambda)
        lambda = opts.lambda;
    elseif isfield(opts,'lambda_scale') && ~isempty(opts.lambda_scale)
        lambda = opts.lambda_scale/sqrt(max(m,n));
    else
        lambda = 1/sqrt(max(m,n));
    end
    maxIter = getfield_default(opts,'maxIter',1000);
    tol = getfield_default(opts,'tol',1e-7);

    normD = norm(D,'fro') + eps;
    norm2D = norm(D,2);
    normInfD = max(abs(D(:))) / lambda;
    dual_norm = max(norm2D,normInfD) + eps;
    Y = D / dual_norm;
    mu = 1.25/(norm2D+eps);
    mu_bar = mu*1e7;
    rho = 1.5;
    L = zeros(m,n,class(D)); S = zeros(m,n,class(D)); relres = Inf;
    for iter = 1:maxIter
        L = svd_threshold(D-S+(1/mu)*Y,1/mu);
        S = soft_threshold_complex(D-L+(1/mu)*Y,lambda/mu);
        Z = D-L-S;
        relres = norm(Z,'fro')/normD;
        if relres < tol, break; end
        Y = Y + mu*Z;
        mu = min(mu*rho,mu_bar);
    end
    info.iter = iter; info.relres = relres; info.lambda = lambda; info.rank = rank_estimate(L);
end

function val = getfield_default(s, field, default)
    if isfield(s,field) && ~isempty(s.(field)), val = s.(field); else, val = default; end
end

function X = soft_threshold_complex(Y,tau)
    mag = abs(Y);
    X = max(mag-tau,0)./(mag+eps).*Y;
end

function L = svd_threshold(M,tau)
    [U,S,V] = svd(M,'econ');
    s = diag(S); s2 = max(s-tau,0); r = sum(s2>0);
    if r==0, L = zeros(size(M),class(M)); else, L = U(:,1:r)*diag(s2(1:r))*V(:,1:r)'; end
end

function Lr = rank_project(M,r)
    if isempty(r), Lr=M; return; end
    [U,S,V] = svd(M,'econ');
    r = min([r,size(U,2),size(V,2)]);
    if r<=0, Lr=zeros(size(M),class(M)); else, Lr=U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; end
end

function Lr = rank_project_energy(M, energy_keep, rmax)
    [U,S,V] = svd(M,'econ');
    s = diag(S);
    if isempty(s) || sum(s.^2)==0, Lr=zeros(size(M),class(M)); return; end
    e = cumsum(s.^2)/sum(s.^2);
    r = find(e>=energy_keep,1,'first');
    r = min([r,rmax,size(U,2),size(V,2)]);
    Lr = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
end

function r = rank_estimate(M)
    s = svd(M,'econ');
    if isempty(s) || s(1)==0, r=0; else, r=sum(s>max(size(M))*eps(s(1))); end
end

function ntfft = choose_ntfft(nt,pad_factor)
    if nargin<2 || isempty(pad_factor), pad_factor=1; end
    ntfft = 2^nextpow2(max(nt,ceil(pad_factor*nt)));
end

function [freqs,w] = make_frequency_taper(ntfft,dt,taper)
    fs = 1/dt; fpos = (0:ntfft-1)'/(ntfft*dt); fnyq = fs/2;
    fabs = fpos; neg = fpos>fnyq; fabs(neg)=fs-fpos(neg);
    w = zeros(ntfft,1);
    f0=taper.flow_stop; f1=taper.flow_pass; f2=taper.fhigh_pass; f3=taper.fhigh_stop;
    w(fabs>=f1 & fabs<=f2)=1;
    idx = fabs>f0 & fabs<f1;
    if f1>f0
        x=(fabs(idx)-f0)/(f1-f0); w(idx)=0.5-0.5*cos(pi*x);
    else
        w(fabs<=f1)=1;
    end
    idx = fabs>f2 & fabs<f3;
    if f3>f2
        x=(fabs(idx)-f2)/(f3-f2); w(idx)=0.5+0.5*cos(pi*x);
    end
    freqs = fpos;
end

function Df = fft_time_padded(data,ntfft)
    [nt,nx,ny] = size(data);
    Dpad = zeros(ntfft,nx,ny);
    Dpad(1:nt,:,:) = data;
    Df = fft(Dpad,[],1);
end

function data = ifft_time_crop(Df,nt)
    dpad = real(ifft(Df,[],1));
    data = dpad(1:nt,:,:);
end

function q = quality_db(clean_ref, estimate)
    err = estimate-clean_ref;
    q = 10*log10(sum(clean_ref(:).^2)/(sum(err(:).^2)+eps));
end

function m = mean_safe(x)
    if isempty(x), m = NaN; else, m = mean(x(:),'omitnan'); end
end

function data = make_lowrank_fxy_synthetic(nt,nx,ny,dt)
    t = (0:nt-1)'*dt;
    [X,Y] = meshgrid(1:nx,1:ny); X=X'; Y=Y';
    data = zeros(nt,nx,ny);
    events = [0.16, 0.0035, 0.0010, 1.00, 25; ...
              0.34,-0.0020, 0.0018,-0.85, 22; ...
              0.58, 0.0015,-0.0025, 0.75, 18; ...
              0.72, 0.0030, 0.0005,-0.60, 20];
    x0=(nx+1)/2; y0=(ny+1)/2;
    for ie=1:size(events,1)
        t0=events(ie,1); sx=events(ie,2); sy=events(ie,3); amp=events(ie,4); fdom=events(ie,5);
        tau = t0 + sx*(X-x0) + sy*(Y-y0);
        for ix=1:nx
            for iy=1:ny
                data(:,ix,iy) = data(:,ix,iy) + amp*ricker_wavelet(t-tau(ix,iy),fdom);
            end
        end
    end
end

function w = ricker_wavelet(t,fdom)
    a=(pi*fdom*t).^2;
    w=(1-2*a).*exp(-a);
end

function [noisy, gaussian_noise, erratic_noise] = add_gaussian_and_erratic_noise(clean,snr_db,erratic_fraction,erratic_amp)
    signal_power = mean(clean(:).^2);
    noise_power = signal_power/(10^(snr_db/10));
    gaussian_noise = sqrt(noise_power)*randn(size(clean));
    [nt,nx,ny] = size(clean);
    erratic_noise = zeros(size(clean));
    ntr = nx*ny; n_bad = max(1,round(erratic_fraction*ntr));
    bad_idx = randperm(ntr,n_bad);
    maxamp = max(abs(clean(:))); t = linspace(0,1,nt)';
    for k=1:n_bad
        [ix,iy] = ind2sub([nx,ny],bad_idx(k));
        if rand<0.7
            burst = conv(randn(nt,1),ones(5,1)/5,'same');
        else
            f0 = 12+18*rand; phase=2*pi*rand; burst=cos(2*pi*f0*t+phase);
        end
        burst = burst/(max(abs(burst))+eps);
        erratic_noise(:,ix,iy) = erratic_amp*maxamp*(0.7+0.6*rand)*burst;
    end
    noisy = clean + gaussian_noise + erratic_noise;
end

function wiggle_plot_fixed(d,t,amp_ref)
    [nt,nx] = size(d);
    if nargin<2 || isempty(t), t=(0:nt-1)'; end
    if nargin<3 || isempty(amp_ref), amp_ref=max(abs(d(:)))+eps; end
    d = d/(amp_ref+eps); scale=0.45; hold on;
    for ix=1:nx
        tr=d(:,ix); x=ix+scale*tr; plot(x,t,'k','LineWidth',0.6);
        pos=tr; pos(pos<0)=0; xp=ix+scale*pos;
        fill([ix*ones(nt,1); flipud(xp)],[t(:); flipud(t(:))],'k','EdgeColor','none','FaceAlpha',0.25);
    end
    hold off; set(gca,'YDir','reverse','Box','on','FontName','Arial','FontSize',9);
    xlabel('trace number'); ylabel('time (s)'); xlim([0 nx+1]); ylim([min(t) max(t)]);
end
