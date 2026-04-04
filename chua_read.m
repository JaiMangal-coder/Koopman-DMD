% =============================================================
%  Chua Circuit Data Reader & Plotter
%  Reads CSV from Arduino serial log, plots attractor, does
%  basic analysis (FFT, phase portraits, statistics)
% =============================================================

clear; clc; close all;

%% ---- 1. READ DATA ------------------------------------------

% Option A: Read from saved CSV file
%   - In Arduino Serial Monitor: click save or copy-paste to .csv
%   - In Python: python -m serial.tools.miniterm COM3 115200 > chua_data.csv
filename = 'chua_data.csv';

% Option B: Read live from serial port (uncomment to use)
% data = read_serial_live('COM3', 115200, 5000);  % see function below

% --- Read CSV (skip comment lines starting with #) -----------
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end

raw = textscan(fid, '%f %f %f %f', ...
    'Delimiter', ',', ...
    'CommentStyle', '#', ...
    'HeaderLines', 1);   % skip the header row
fclose(fid);

t = raw{1} / 1000.0;   % ms -> seconds
x = raw{2};            % V_C1
y = raw{3};            % V_C2
z = raw{4};            % I_L sense voltage

N = length(t);
fprintf('Loaded %d samples\n', N);
fprintf('Duration: %.3f s\n', t(end) - t(1));

% Compute uniform sample interval (Arduino timing is not perfect)
dt_nominal = median(diff(t));
fs = 1 / dt_nominal;
fprintf('Median sample interval: %.4f ms  (fs ~ %.1f Hz)\n', ...
    dt_nominal*1000, fs);

%% ---- 2. RESAMPLE TO UNIFORM GRID ---------------------------
% Arduino micros() drift means timestamps are not perfectly uniform.
% Resample onto a uniform grid before FFT or further analysis.

t_uniform = linspace(t(1), t(end), N)';
x_u = interp1(t, x, t_uniform, 'pchip');
y_u = interp1(t, y, t_uniform, 'pchip');
z_u = interp1(t, z, t_uniform, 'pchip');

%% ---- 3. PHASE PORTRAITS ------------------------------------

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

%% ---- 4. 3D ATTRACTOR ---------------------------------------

figure('Name', 'Chua 3D Attractor')
plot3(x_u, y_u, z_u, 'k', 'LineWidth', 0.2)
xlabel('V_{C1}  (V)')
ylabel('V_{C2}  (V)')
zlabel('V_{IL}  (V)')
title('Chua Double Scroll — 3D')
grid on; axis tight; view(35, 25)
rotate3d on   % enable mouse rotation

%% ---- 5. TIME SERIES ----------------------------------------

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

%% ---- 6. POWER SPECTRAL DENSITY -----------------------------

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

%% ---- 7. BASIC STATISTICS -----------------------------------

fprintf('\n--- Signal Statistics ---\n')
fprintf('%8s  %8s  %8s  %8s  %8s\n', 'Channel', 'Mean', 'Std', 'Min', 'Max')
fprintf('%8s  %8.4f  %8.4f  %8.4f  %8.4f\n', 'x (Vc1)', mean(x_u), std(x_u), min(x_u), max(x_u))
fprintf('%8s  %8.4f  %8.4f  %8.4f  %8.4f\n', 'y (Vc2)', mean(y_u), std(y_u), min(y_u), max(y_u))
fprintf('%8s  %8.4f  %8.4f  %8.4f  %8.4f\n', 'z (IL) ', mean(z_u), std(z_u), min(z_u), max(z_u))

%% ---- 8. LIVE SERIAL READ FUNCTION --------------------------
% Uncomment Option B above and use this to read directly from Arduino.
% Call: data = read_serial_live('COM3', 115200, 5000);
% On Mac/Linux use port like '/dev/ttyUSB0' or '/dev/tty.usbmodem14101'

function data = read_serial_live(port, baud, n_samples)
    s = serialport(port, baud);
    configureTerminator(s, 'LF');
    flush(s);

    fprintf('Waiting for Arduino... send ''s'' to start\n');
    writeline(s, 's');   % send start command

    data = zeros(n_samples, 4);
    k = 0;

    while k < n_samples
        line = readline(s);
        line = strtrim(line);
        if isempty(line) || line(1) == '#', continue; end
        vals = sscanf(line, '%f,%f,%f,%f');
        if numel(vals) == 4
            k = k + 1;
            data(k,:) = vals';
            if mod(k, 500) == 0
                fprintf('  %d / %d samples\n', k, n_samples);
            end
        end
    end

    writeline(s, 'q');   % send stop command
    clear s;
    fprintf('Done. Collected %d samples.\n', k);
end
