% =============================================================
%  Chua Circuit Data Reader, Live Plotter & Analyzer
%  - Offline mode: reads saved CSV and runs analysis
%  - Live mode: reads CSV from Arduino over serial, plots in real time,
%               and optionally saves to CSV
% =============================================================

clear; clc; close all;

%% ---- 0. MODE SELECTION -------------------------------------
% Choose:
%   'offline' -> read saved CSV file and analyze
%   'live'    -> stream from Arduino, plot in real time, then analyze

mode = 'live';   % <-- change to 'offline' if reading a saved file

%% ---- 1. SETTINGS -------------------------------------------

% Offline file
filename = 'chua_data.csv';

% Live serial settings
port = 'COM3';          % e.g. 'COM3' on Windows, '/dev/ttyUSB0' on Linux
baud = 115200;
n_samples = 5000;       % how many samples to capture in live mode
save_live_csv = true;   % save streamed data to CSV
live_save_filename = 'chua_live_capture.csv';

% Live plotting settings
maxPoints = 2000;       % rolling window for time-series display

%% ---- 2. LOAD OR STREAM DATA --------------------------------

switch lower(mode)
    case 'offline'
        [t, x, y, z] = read_saved_csv(filename);

    case 'live'
        data = read_serial_live_plot(port, baud, n_samples, maxPoints, ...
            save_live_csv, live_save_filename);

        t = data(:,1) / 1000.0;   % ms -> s
        x = data(:,2);
        y = data(:,3);
        z = data(:,4);

    otherwise
        error('Unknown mode. Use ''offline'' or ''live''.');
end

%% ---- 3. BASIC INFO -----------------------------------------

N = length(t);
fprintf('Loaded %d samples\n', N);
fprintf('Duration: %.3f s\n', t(end) - t(1));

% Compute uniform sample interval (Arduino timing is not perfect)
dt_nominal = median(diff(t));
fs = 1 / dt_nominal;
fprintf('Median sample interval: %.4f ms  (fs ~ %.1f Hz)\n', ...
    dt_nominal*1000, fs);

%% ---- 4. RESAMPLE TO UNIFORM GRID ---------------------------
% Arduino micros() drift means timestamps are not perfectly uniform.
% Resample onto a uniform grid before FFT or further analysis.

t_uniform = linspace(t(1), t(end), N)';
x_u = interp1(t, x, t_uniform, 'pchip');
y_u = interp1(t, y, t_uniform, 'pchip');
z_u = interp1(t, z, t_uniform, 'pchip');

%% ---- 5. PHASE PORTRAITS ------------------------------------

figure('Name', 'Chua Attractor - Phase Portraits', ...
    'Position', [100 100 1200 400]);

subplot(1,3,1)
plot(x_u, y_u, 'b', 'LineWidth', 0.3)
xlabel('V_{C1}  (V)'); ylabel('V_{C2}  (V)')
title('x-y  (Vc1 vs Vc2)')
grid on; axis tight

subplot(1,3,2)
plot(x_u, z_u, 'r', 'LineWidth', 0.3)
xlabel('V_{C1}  (V)'); ylabel('V_{IL}  (V)')
title('x-z  (Vc1 vs IL)')
grid on; axis tight

subplot(1,3,3)
plot(y_u, z_u, 'm', 'LineWidth', 0.3)
xlabel('V_{C2}  (V)'); ylabel('V_{IL}  (V)')
title('y-z  (Vc2 vs IL)')
grid on; axis tight

sgtitle('Chua Double Scroll Attractor')

%% ---- 6. 3D ATTRACTOR ---------------------------------------

figure('Name', 'Chua 3D Attractor')
plot3(x_u, y_u, z_u, 'k', 'LineWidth', 0.2)
xlabel('V_{C1}  (V)')
ylabel('V_{C2}  (V)')
zlabel('V_{IL}  (V)')
title('Chua Double Scroll — 3D')
grid on; axis tight; view(35, 25)
rotate3d on   % enable mouse rotation

%% ---- 7. TIME SERIES ----------------------------------------

figure('Name', 'Time Series')
subplot(3,1,1)
plot(t_uniform, x_u, 'b', 'LineWidth', 0.5)
ylabel('V_{C1} (V)'); title('x — V_{C1}'); grid on; xlim([t(1) t(end)])

subplot(3,1,2)
plot(t_uniform, y_u, 'r', 'LineWidth', 0.5)
ylabel('V_{C2} (V)'); title('y — V_{C2}'); grid on; xlim([t(1) t(end)])

subplot(3,1,3)
plot(t_uniform, z_u, 'm', 'LineWidth', 0.5)
ylabel('V_{IL} (V)'); title('z — I_L sense'); grid on; xlim([t(1) t(end)])
xlabel('Time (s)')

%% ---- 8. POWER SPECTRAL DENSITY -----------------------------

figure('Name', 'Power Spectral Density')
signals = {x_u, y_u, z_u};
names   = {'V_{C1}', 'V_{C2}', 'V_{IL}'};
colors  = {'b', 'r', 'm'};

for k = 1:3
    [pxx, f] = pwelch(signals{k}, hamming(512), 256, 512, fs);
    semilogy(f, pxx, colors{k}, 'LineWidth', 1.2); hold on
end
xlabel('Frequency (Hz)'); ylabel('PSD (V^2/Hz)')
title('Power Spectral Density')
legend(names); grid on
% Broadband spectrum is a signature of chaos

%% ---- 9. BASIC STATISTICS -----------------------------------

fprintf('\n--- Signal Statistics ---\n')
fprintf('%8s  %8s  %8s  %8s  %8s\n', 'Channel', 'Mean', 'Std', 'Min', 'Max')
fprintf('%8s  %8.4f  %8.4f  %8.4f  %8.4f\n', 'x (Vc1)', mean(x_u), std(x_u), min(x_u), max(x_u))
fprintf('%8s  %8.4f  %8.4f  %8.4f  %8.4f\n', 'y (Vc2)', mean(y_u), std(y_u), min(y_u), max(y_u))
fprintf('%8s  %8.4f  %8.4f  %8.4f  %8.4f\n', 'z (IL) ', mean(z_u), std(z_u), min(z_u), max(z_u))

%% ---- FUNCTIONS ---------------------------------------------

function [t, x, y, z] = read_saved_csv(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    raw = textscan(fid, '%f %f %f %f', ...
        'Delimiter', ',', ...
        'CommentStyle', '#', ...
        'HeaderLines', 1);
    fclose(fid);

    t = raw{1} / 1000.0;   % ms -> seconds
    x = raw{2};
    y = raw{3};
    z = raw{4};
end

function data = read_serial_live_plot(port, baud, n_samples, maxPoints, save_csv, save_filename)
    % Live serial read + oscilloscope-style plotting + optional CSV save
    %
    % Expected Arduino format:
    %   # comments...
    %   time_ms,x_V,y_V,z_V
    %   12.34,-0.123,0.456,0.789

    s = serialport(port, baud);
    configureTerminator(s, 'LF');
    flush(s);
    pause(2);  % allow Arduino reset on serial connection

    if save_csv
        fid = fopen(save_filename, 'w');
        if fid == -1
            error('Cannot create save file: %s', save_filename);
        end
        fprintf(fid, '# Chua live capture\n');
        fprintf(fid, 'time_ms,x_V,y_V,z_V\n');
    else
        fid = -1;
    end

    % Create live figure
    fig = figure('Name', 'Chua Live Plot', ...
        'Position', [100 100 1200 700], ...
        'NumberTitle', 'off');

    tiledlayout(2,1);

    % Top: rolling time traces
    ax1 = nexttile;
    hx = animatedline('DisplayName', 'x = Vc1');
    hy = animatedline('DisplayName', 'y = Vc2');
    hz = animatedline('DisplayName', 'z = Il');
    grid(ax1, 'on');
    xlabel(ax1, 'Time (ms)');
    ylabel(ax1, 'Voltage (V)');
    title(ax1, 'Live Signals');
    legend(ax1, 'show');

    % Bottom: live XY portrait
    ax2 = nexttile;
    hxy = animatedline('DisplayName', 'x vs y');
    grid(ax2, 'on');
    xlabel(ax2, 'x = Vc1 (V)');
    ylabel(ax2, 'y = Vc2 (V)');
    title(ax2, 'Live Phase Portrait');

    drawnow;

    fprintf('Waiting for Arduino... sending ''s'' to start\n');
    writeline(s, 's');

    data = zeros(n_samples, 4);
    k = 0;

    timeBuf = [];
    xBuf = [];
    yBuf = [];
    zBuf = [];

    while k < n_samples && ishandle(fig)
        line = strtrim(readline(s));

        if isempty(line) || startsWith(line, "#") || startsWith(line, "time_ms")
            continue;
        end

        vals = sscanf(line, '%f,%f,%f,%f');
        if numel(vals) ~= 4
            continue;
        end

        k = k + 1;
        data(k,:) = vals(:)';

        t = vals(1);
        x = vals(2);
        y = vals(3);
        z = vals(4);

        if fid ~= -1
            fprintf(fid, '%.3f,%.4f,%.4f,%.4f\n', t, x, y, z);
        end

        timeBuf(end+1) = t; %#ok<AGROW>
        xBuf(end+1) = x; %#ok<AGROW>
        yBuf(end+1) = y; %#ok<AGROW>
        zBuf(end+1) = z; %#ok<AGROW>

        % Rolling window for time-domain display
        if numel(timeBuf) > maxPoints
            timeBuf = timeBuf(end-maxPoints+1:end);
            xBuf = xBuf(end-maxPoints+1:end);
            yBuf = yBuf(end-maxPoints+1:end);
            zBuf = zBuf(end-maxPoints+1:end);

            clearpoints(hx); clearpoints(hy); clearpoints(hz); clearpoints(hxy);
            addpoints(hx, timeBuf, xBuf);
            addpoints(hy, timeBuf, yBuf);
            addpoints(hz, timeBuf, zBuf);
            addpoints(hxy, xBuf, yBuf);
        else
            addpoints(hx, t, x);
            addpoints(hy, t, y);
            addpoints(hz, t, z);
            addpoints(hxy, x, y);
        end

        if mod(k, 50) == 0
            drawnow limitrate;
        end

        if mod(k, 500) == 0
            fprintf('  %d / %d samples\n', k, n_samples);
        end
    end

    % Stop Arduino
    try
        writeline(s, 'q');
    catch
    end

    if fid ~= -1
        fclose(fid);
        fprintf('Saved live capture to %s\n', save_filename);
    end

    data = data(1:k,:);
    clear s;
    fprintf('Done. Collected %d samples.\n', k);
end
