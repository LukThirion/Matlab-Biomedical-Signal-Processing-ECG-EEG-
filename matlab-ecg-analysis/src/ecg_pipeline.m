function ecg_pipeline
%% ECG Signal Processing Pipeline (toolbox-friendly)
% Author: Lukha Thirion
% Loads ECG (or makes synthetic), denoises, detects R-peaks, computes HR/HRV,
% and saves plots + results to ../figs/.
%
% Works with or without Signal Processing Toolbox:
% - If available: uses butter/filtfilt/findpeaks/iirnotch
% - If not: uses moving-average bandpass, forward-backward filtering,
%           and a simple local-maximum peak detector.

%% -------- Paths and folders --------
thisFile  = mfilename('fullpath');
projRoot  = fileparts(fileparts(thisFile));   % go up from /src
dataFolder = fullfile(projRoot, 'data');
figFolder  = fullfile(projRoot, 'figs');
if ~exist(figFolder, 'dir'), mkdir(figFolder); end

%% -------- Settings you can tweak --------
Fs_default = 250;        % CSV sampling rate if unknown
bp_low  = 0.5;           % bandpass low cut (Hz)
bp_high = 40;            % bandpass high cut (Hz)
notchHz = [];            % set to 60 or 50 to enable; [] disables (safe w/out toolbox)
minPeakHeightFrac = 0.35; % threshold as fraction of max after filtering
minPeakDistance_s = 0.25; % min RR interval (sec) ~ 240 bpm upper bound

%% -------- 0) Capability checks (which functions exist?) --------
HAS_BUTTER   = exist('butter','file')   == 2;
HAS_FILTFILT = exist('filtfilt','file') == 2;
HAS_FINDPEAKS= exist('findpeaks','file')== 2;
HAS_IIRNOTCH = exist('iirnotch','file') == 2;

%% -------- 1) Load ECG or generate synthetic --------
[ecg, Fs] = tryLoadECG(dataFolder, Fs_default);
t = (0:numel(ecg)-1)'/Fs;

%% -------- 2) Detrend & Filter --------
% Detrend (remove slow drift)
ecg_dt = detrend(ecg, 'linear');

% --- Bandpass:
if HAS_BUTTER
    % 2nd-order Butterworth bandpass + zero-phase if available
    Wn = [bp_low bp_high]/(Fs/2);
    Wn(1) = max(Wn(1), 1e-6);          % guard rails
    [b,a] = butter(2, Wn, 'bandpass');
    if HAS_FILTFILT
        ecg_bp = filtfilt(b,a,ecg_dt);
    else
        ecg_bp = zfilter_fb(b,a,ecg_dt); % forward-backward fallback
    end
else
    % Toolbox-free movable-average "bandpass":
    % 1) Low-pass (short window) to remove high-frequency noise
    Nlp = max(3, round(Fs*0.04));             % ~40 ms window
    ecg_lp_short = movmean(ecg_dt, Nlp);
    % 2) Remove baseline (very long window) -> high-pass effect
    Nbase = max(Nlp+2, round(Fs*0.7));         % ~700 ms window
    baseline = movmean(ecg_dt, Nbase);
    ecg_bp = ecg_lp_short - baseline;          % crude bandpass
end

% --- Optional notch (60/50 Hz)
if ~isempty(notchHz)
    if HAS_IIRNOTCH
        wo = notchHz/(Fs/2);
        bw = wo/35;
        [bn,an] = iirnotch(wo,bw);
        if HAS_FILTFILT
            ecg_filt = filtfilt(bn,an,ecg_bp);
        else
            ecg_filt = zfilter_fb(bn,an,ecg_bp);
        end
    else
        % Simple manual notch using pole/zero pair (fallback)
        w0 = 2*pi*notchHz/Fs;                 % rad/sample
        r  = 0.97;                            % pole radius (notch width)
        b  = [1, -2*cos(w0), 1];              % zeros on unit circle
        a  = [1, -2*r*cos(w0), r^2];          % poles just inside
        ecg_filt = zfilter_fb(b,a,ecg_bp);
    end
else
    ecg_filt = ecg_bp;
end

%% -------- 3) R-peak detection --------
minPeakHeight  = minPeakHeightFrac * max(ecg_filt);
minPeakDistance = round(minPeakDistance_s * Fs);

if HAS_FINDPEAKS
    [pkVals, pkLocs] = findpeaks(ecg_filt, 'MinPeakHeight', minPeakHeight, ...
        'MinPeakDistance', minPeakDistance);
else
    [pkLocs, pkVals] = simple_findpeaks(ecg_filt, minPeakHeight, minPeakDistance);
end

R_times = pkLocs / Fs;
RR      = diff(R_times);             % seconds
HR_inst = 60 ./ RR;                  % bpm

% HRV metrics
if ~isempty(RR)
    meanHR   = 60 / mean(RR);
    SDNN_ms  = std(RR) * 1000;                  % ms
    RMSSD_ms = sqrt(mean(diff(RR).^2)) * 1000;  % ms
else
    meanHR = NaN; SDNN_ms = NaN; RMSSD_ms = NaN;
end

%% -------- 4) Report metrics --------
fprintf('\n--- Metrics ---\n');
fprintf('Beats detected: %d\n', numel(pkLocs));
if ~isnan(meanHR)
    fprintf('Mean HR: %.1f bpm | SDNN: %.1f ms | RMSSD: %.1f ms\n', ...
        meanHR, SDNN_ms, RMSSD_ms);
else
    fprintf('Not enough beats for HRV metrics.\n');
end

%% -------- 5) Plots (saved to figs/) --------
% A) Raw → Detrended → Filtered
figure('Name','ECG—Raw vs Detrended vs Filtered','Position',[100 100 900 520]);
subplot(3,1,1); plot(t, ecg);      title('Raw ECG');        xlabel('Time (s)'); ylabel('Amp');
subplot(3,1,2); plot(t, ecg_dt);   title('Detrended ECG');  xlabel('Time (s)'); ylabel('Amp');
subplot(3,1,3); plot(t, ecg_filt); title('Filtered ECG');   xlabel('Time (s)'); ylabel('Amp');
saveas(gcf, fullfile(figFolder, '01_ecg_raw_detrended_filtered.png'));

% B) Filtered with R-peaks
figure('Name','ECG—R-peaks','Position',[100 100 900 300]);
plot(t, ecg_filt); hold on;
plot(R_times, ecg_filt(pkLocs), 'rv', 'MarkerFaceColor','r');
xlabel('Time (s)'); ylabel('Amp'); title('Filtered ECG with R-peaks'); grid on;
saveas(gcf, fullfile(figFolder, '02_ecg_rpeaks.png'));

% C) Instantaneous HR (if enough beats)
if numel(HR_inst) >= 2
    tHR = (R_times(1:end-1) + R_times(2:end)) / 2;
    figure('Name','Instantaneous HR','Position',[100 100 900 300]);
    plot(tHR, HR_inst, '-o'); grid on;
    xlabel('Time (s)'); ylabel('HR (bpm)'); title('Instantaneous Heart Rate');
    saveas(gcf, fullfile(figFolder, '03_hr_instantaneous.png'));
end

%% -------- 6) Save results --------
results = struct('meanHR', meanHR, 'SDNN_ms', SDNN_ms, 'RMSSD_ms', RMSSD_ms, ...
                 'R_times', R_times, 'RR_s', RR, 'Fs', Fs);
save(fullfile(figFolder,'results.mat'), 'results');
fprintf('\nSaved figures in figs/ and metrics in figs/results.mat\n');

end % function ecg_pipeline

%% ---------- Helpers (no toolboxes needed) ----------

function y = zfilter_fb(b,a,x)
% Simple forward-backward IIR filtering (fallback if filtfilt is unavailable).
% Note: lacks filtfilt's edge correction but works well for this use.
y = filter(b,a,x);
y = flipud(filter(b,a,flipud(y)));
end

function [locs, vals] = simple_findpeaks(x, minH, minD)
% Very small peak finder:
% - local maxima greater than neighbors
% - above minH
% - at least minD samples apart
n = numel(x);
locs = [];
i = 2;
while i <= n-1
    if x(i) > x(i-1) && x(i) >= x(i+1) && x(i) >= minH
        % enforce min distance
        if isempty(locs) || (i - locs(end)) >= minD
            locs(end+1) = i; %#ok<AGROW>
        else
            % keep the higher one within the window
            if x(i) > x(locs(end))
                locs(end) = i;
            end
        end
        i = i + 1;
    end
    i = i + 1;
end
vals = x(locs);
end

function [ecg, Fs] = tryLoadECG(dataFolder, Fs_default)
% Load first .mat (ecg, Fs) or .csv (column 'ecg'); else generate synthetic ECG.
ecg = []; Fs = [];
matFiles = dir(fullfile(dataFolder, '*.mat'));
csvFiles = dir(fullfile(dataFolder, '*.csv'));
if ~isempty(matFiles)
    f = fullfile(matFiles(1).folder, matFiles(1).name);
    s = load(f);
    if isfield(s,'ecg'), ecg = s.ecg(:); end
    if isfield(s,'Fs'),  Fs  = s.Fs;     end
    if ~isempty(ecg) && ~isempty(Fs)
        fprintf('Loaded MAT: %s\n', matFiles(1).name); return;
    end
elseif ~isempty(csvFiles)
    f = fullfile(csvFiles(1).folder, csvFiles(1).name);
    T = readtable(f);
    if any(strcmpi(T.Properties.VariableNames, 'ecg'))
        ecg = T.ecg(:); Fs = Fs_default;
        fprintf('Loaded CSV: %s (Fs assumed %g Hz)\n', csvFiles(1).name, Fs); return;
    end
end
% Fallback: synthetic ECG
duration = 12;                 % seconds
Fs = 250;                      % Hz
t  = (0:1/Fs:duration-1/Fs)';
HR_bpm = 72; rr = 60/HR_bpm;
base = 0.15*sin(2*pi*0.3*t);   % mild baseline wander
qrs  = 0.8*exp(-((mod(t,rr)-0.02)/0.015).^2);
T    = 0.12*exp(-((mod(t,rr)-0.18)/0.03).^2);
ecg  = base + qrs + T + 0.05*randn(size(t));
fprintf('No data found. Generated synthetic ECG (%gs @ %g Hz).\n', duration, Fs);
end