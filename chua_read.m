% =============================================================
%  Chua Circuit — Multi-Regime Data Capture & Analyzer
%
%  Modes:
%   'live'    — waits for you to type regime commands in the Serial
%               Monitor, captures each regime as you send them,
%               saves each to a labelled CSV, then runs analysis.
%   'offline' — loads previously saved regime CSVs and runs analysis.
%
%  Live workflow:
%   1. Open this script and set port, baud, and the regime order
%      you plan to collect (REGIME_IDS below).
%   2. Run the script.  It opens the serial port and waits.
%   3. In your circuit: dial in Regime 1 (fixed point).
%      The script will prompt you; press Enter in the Command Window
%      to send the regime number + 's' to the Arduino automatically.
%   4. Data streams in, is saved to chua_regime1_fixed_point.csv, and
%      the next prompt appears.
%   5. Repeat for all four regimes.
%   6. After all regimes are collected, analysis figures are generated.
%
%  Offline workflow:
%   Set mode = 'offline'.  The script looks for the four CSV files in
%   the working directory and loads them directly.
% =============================================================

clear; clc; close all;

%% ---- 0. USER SETTINGS ----------------------------------------

mode = 'live';          % 'live' or 'offline'

% Serial port settings (live mode only)
port = '/dev/ttyACM0';          % Windows: 'COM3'; Linux/Mac: '/dev/ttyUSB0'
baud = 2000000;         % R4 USB-CDC ignores this; set high so host doesn't throttle

% Regimes to collect, in order.
% Each entry is [regime_id, label, friendly_name]
REGIMES = {
    1, 'fixed_point',      'Fixed Point';
    2, 'limit_cycle',      'Limit Cycle';
    3, 'period_doubled',   'Period-Doubled';
    4, 'double_scroll',    'Double Scroll (Chaos)';
};

% CSV filenames (same convention in live and offline mode)
csv_filenames = cellfun(@(id, lbl) ...
    sprintf('chua_regime%d_%s.csv', id, lbl), ...
    REGIMES(:,1), REGIMES(:,2), 'UniformOutput', false);

% Live display: rolling window length (samples)
maxPoints = 2000;

%% ---- 1. ACQUIRE DATA ----------------------------------------

n_regimes = size(REGIMES, 1);
regime_data = cell(n_regimes, 1);   % will hold struct(t,x,y,z,label,name)

switch lower(mode)

    % ----------------------------------------------------------
    case 'live'

        s = serialport(port, baud);
        configureTerminator(s, 'LF');
        flush(s);
        pause(2);   % let Arduino reset after serial open

        % Drain the Arduino boot message
        pause(0.5);
        while s.NumBytesAvailable > 0
            readline(s);
        end

        for r = 1:n_regimes
            regime_id   = REGIMES{r,1};
            regime_lbl  = REGIMES{r,2};
            regime_name = REGIMES{r,3};
            out_file    = csv_filenames{r};

            fprintf('\n============================================\n');
            fprintf('  REGIME %d: %s\n', regime_id, regime_name);
            fprintf('============================================\n');
            fprintf('  1. Adjust your circuit to the %s regime.\n', regime_name);
            fprintf('  2. Press Enter here when ready to start sampling.\n');
            input('  [Press Enter to start] ', 's');

            % Send regime select + start to Arduino
            writeline(s, sprintf('%d', regime_id));
            pause(0.1);
            writeline(s, 's');

            % The R4 sketch buffers samples in RAM and dumps the full CSV
            % block AFTER capture completes (between DATA START / DATA END).
            % We wait for that block, then plot it all at once.
            fprintf('  Arduino sampling (silent during capture)...\n');
            fprintf('  Waiting for DATA START block...\n');

            tBuf=[]; xBuf=[]; yBuf=[]; zBuf=[];
            in_data_block = false;
            done          = false;
            timeout_s     = 60;   % max wait for capture + dump
            t_wait        = tic;

            while ~done && toc(t_wait) < timeout_s
                if s.NumBytesAvailable == 0
                    pause(0.01);
                    continue;
                end
                raw = strtrim(readline(s));
                if isempty(raw), continue; end

                if startsWith(raw, '#')
                    fprintf('  Arduino: %s\n', raw);
                    if contains(raw, 'DATA START')
                        in_data_block = true;
                    elseif contains(raw, 'DATA END') || contains(raw, 'ABORTED')
                        done = true;
                    end
                    continue;
                end

                if startsWith(raw, 'time_ms'), continue; end

                if ~in_data_block, continue; end   % ignore anything before block

                vals = sscanf(raw, '%f,%f,%f,%f');
                if numel(vals) ~= 4, continue; end

                tBuf(end+1) = vals(1); %#ok<AGROW>
                xBuf(end+1) = vals(2); %#ok<AGROW>
                yBuf(end+1) = vals(3); %#ok<AGROW>
                zBuf(end+1) = vals(4); %#ok<AGROW>

                if mod(numel(tBuf), 200) == 0
                    fprintf('  %d samples received\n', numel(tBuf));
                end
            end

            if toc(t_wait) >= timeout_s
                warning('Timeout waiting for regime %d data.', regime_id);
            end

            % Save to CSV
            fid = fopen(out_file, 'w');
            fprintf(fid, '# Chua circuit — regime %d: %s\n', regime_id, regime_lbl);
            fprintf(fid, 'time_ms,x_V,y_V,z_V\n');
            for ii = 1:numel(tBuf)
                fprintf(fid, '%.3f,%.4f,%.4f,%.4f\n', tBuf(ii),xBuf(ii),yBuf(ii),zBuf(ii));
            end
            fclose(fid);

            % Static plot of the captured data
            fig = figure('Name', sprintf('Regime %d: %s', regime_id, regime_name), ...
                'Position', [80 80 1100 650], 'NumberTitle', 'off');
            tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');

            ax1 = nexttile(tl);
            plot(ax1, tBuf, xBuf, 'Color', [0.20 0.45 0.80], 'LineWidth', 0.5); hold(ax1,'on');
            plot(ax1, tBuf, yBuf, 'Color', [0.80 0.20 0.20], 'LineWidth', 0.5);
            plot(ax1, tBuf, zBuf, 'Color', [0.55 0.20 0.75], 'LineWidth', 0.5);
            xlabel(ax1, 'Time (ms)'); ylabel(ax1, 'Voltage (V)');
            title(ax1, sprintf('Time series — %s', regime_name));
            legend(ax1, {'x=Vc1','y=Vc2','z=IL'}); grid(ax1, 'on');

            ax2 = nexttile(tl);
            plot(ax2, xBuf, yBuf, 'Color', [0.10 0.55 0.40], 'LineWidth', 0.4);
            xlabel(ax2, 'x = Vc1 (V)'); ylabel(ax2, 'y = Vc2 (V)');
            title(ax2, 'Phase portrait x–y'); grid(ax2, 'on');
            drawnow;
            fprintf('  Saved to %s  (%d samples)\n', out_file, numel(tBuf));

            % Store in cell array for later analysis
            t_s = tBuf(:) / 1000.0;
            regime_data{r} = struct('t', t_s, 'x', xBuf(:), 'y', yBuf(:), ...
                'z', zBuf(:), 'label', regime_lbl, 'name', regime_name, ...
                'id', regime_id);
        end

        % Close serial
        try; writeline(s, 'q'); catch; end
        clear s;
        fprintf('\nAll regimes captured. Running analysis...\n');

    % ----------------------------------------------------------
    case 'offline'

        fprintf('Loading saved CSV files...\n');
        for r = 1:n_regimes
            fname = csv_filenames{r};
            if ~isfile(fname)
                error('File not found: %s\nRun in live mode first, or check filenames.', fname);
            end
            fid = fopen(fname, 'r');
            raw = textscan(fid, '%f %f %f %f', 'Delimiter', ',', ...
                'CommentStyle', '#', 'HeaderLines', 1);
            fclose(fid);
            t_s = raw{1} / 1000.0;
            regime_data{r} = struct('t', t_s, 'x', raw{2}, 'y', raw{3}, ...
                'z', raw{4}, 'label', REGIMES{r,2}, 'name', REGIMES{r,3}, ...
                'id', REGIMES{r,1});
            fprintf('  Loaded %d samples from %s\n', numel(t_s), fname);
        end

    otherwise
        error('Unknown mode ''%s''. Use ''live'' or ''offline''.', mode);
end

%% ---- 2. RESAMPLE TO UNIFORM GRID ----------------------------
% Arduino micros() drift makes timestamps slightly non-uniform;
% resample to a uniform grid before spectral analysis.

for r = 1:n_regimes
    d  = regime_data{r};
    N  = numel(d.t);
    tu = linspace(d.t(1), d.t(end), N)';
    regime_data{r}.t_u = tu;
    regime_data{r}.x_u = interp1(d.t, d.x, tu, 'pchip');
    regime_data{r}.y_u = interp1(d.t, d.y, tu, 'pchip');
    regime_data{r}.z_u = interp1(d.t, d.z, tu, 'pchip');
    regime_data{r}.dt  = median(diff(d.t));
    regime_data{r}.fs  = 1 / regime_data{r}.dt;
end

%% ---- 3. ANALYSIS PLOTS  -------------------------------------

colors = { [0.20 0.45 0.80],   % blue  — regime 1
           [0.10 0.60 0.38],   % green — regime 2
           [0.85 0.58 0.10],   % amber — regime 3
           [0.75 0.18 0.18] }; % red   — regime 4

% 3a. Phase portraits (x vs y) — all four regimes in one figure
figure('Name', 'Phase Portraits — All Regimes', 'Position', [100 100 1300 350]);
for r = 1:n_regimes
    d = regime_data{r};
    subplot(1, 4, r)
    plot(d.x_u, d.y_u, 'Color', colors{r}, 'LineWidth', 0.3)
    xlabel('V_{C1} (V)'); ylabel('V_{C2} (V)')
    title(sprintf('R%d: %s', d.id, d.name), 'FontSize', 9)
    grid on; axis tight
end
sgtitle('Chua Circuit Phase Portraits — x vs y')

% 3b. 3-D attractors — separate figure per regime
for r = 1:n_regimes
    d = regime_data{r};
    figure('Name', sprintf('3D Attractor — %s', d.name), ...
        'Position', [120+r*30 120+r*20 600 500]);
    plot3(d.x_u, d.y_u, d.z_u, 'Color', colors{r}, 'LineWidth', 0.25)
    xlabel('V_{C1} (V)'); ylabel('V_{C2} (V)'); zlabel('I_L (V)')
    title(sprintf('R%d: %s — 3D attractor', d.id, d.name))
    grid on; axis tight; view(35, 25); rotate3d on
end

% 3c. Time series — all channels, one regime per figure
for r = 1:n_regimes
    d = regime_data{r};
    figure('Name', sprintf('Time Series — %s', d.name), ...
        'Position', [140+r*20 100 900 500]);
    sigs = {d.x_u, d.y_u, d.z_u};
    ylabels = {'V_{C1} (V)', 'V_{C2} (V)', 'I_L sense (V)'};
    chnames = {'x  (Vc1)', 'y  (Vc2)', 'z  (IL)'};
    for k = 1:3
        subplot(3,1,k)
        plot(d.t_u, sigs{k}, 'Color', colors{r}, 'LineWidth', 0.5)
        ylabel(ylabels{k}); title(chnames{k}); grid on
        xlim([d.t_u(1) d.t_u(end)])
    end
    xlabel('Time (s)')
    sgtitle(sprintf('R%d: %s — Time Series', d.id, d.name))
end

% 3d. Power spectral density — all regimes overlaid
figure('Name', 'PSD — All Regimes', 'Position', [100 100 900 500]);
hold on
legends_psd = cell(n_regimes, 1);
for r = 1:n_regimes
    d = regime_data{r};
    [pxx, f] = pwelch(d.x_u, hamming(512), 256, 512, d.fs);
    semilogy(f, pxx, 'Color', colors{r}, 'LineWidth', 1.3)
    legends_psd{r} = sprintf('R%d: %s', d.id, d.name);
end
xlabel('Frequency (Hz)'); ylabel('PSD (V^2/Hz)')
title('Power Spectral Density — V_{C1} channel, all regimes')
legend(legends_psd); grid on
% Discrete peaks -> periodic; broadband floor -> chaos

%% ---- 4. STATISTICS TABLE ------------------------------------

fprintf('\n%s\n', repmat('=', 1, 72));
fprintf('  Signal statistics\n');
fprintf('%s\n', repmat('=', 1, 72));
fprintf('%-28s %8s %8s %8s %8s %8s\n', ...
    'Regime / Channel', 'Mean', 'Std', 'Min', 'Max', 'RMS');
fprintf('%s\n', repmat('-', 1, 72));

for r = 1:n_regimes
    d = regime_data{r};
    sigs = {d.x_u, d.y_u, d.z_u};
    ch   = {'x (Vc1)', 'y (Vc2)', 'z (IL) '};
    for k = 1:3
        v   = sigs{k};
        rms = sqrt(mean(v.^2));
        label = '';
        if k == 1
            label = sprintf('R%d %-20s', d.id, d.name);
        end
        fprintf('%-28s %8.4f %8.4f %8.4f %8.4f %8.4f   [%s]\n', ...
            label, mean(v), std(v), min(v), max(v), rms, ch{k});
    end
    fprintf('%s\n', repmat('-', 1, 72));
end

fprintf('\nAll figures generated.  CSVs saved in working directory.\n');
