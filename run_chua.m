% RUN_CHUA  Apply EDMD to all four Chua circuit experimental regimes
%
% Data files (in chua/ subdirectory), each with columns: time_ms, x_V, y_V, z_V
%   chua_regime1_fixed_point.csv
%   chua_regime2_limit_cycle.csv
%   chua_regime3_period_doubled.csv
%   chua_regime4_double_scroll.csv

clear; clc;

%% Regime definitions
regimes = { ...
    'chua/chua_regime1_fixed_point.csv',    'Fixed Point'; ...
    'chua/chua_regime2_limit_cycle.csv',    'Limit Cycle'; ...
    'chua/chua_regime3_period_doubled.csv', 'Period-Doubled'; ...
    'chua/chua_regime4_double_scroll.csv',  'Double Scroll'; ...
};
n_regimes = size(regimes, 1);

%% Dictionary options
dict_opts.type      = 'full';
dict_opts.degree    = 3;
dict_opts.n_centers = 20;
dict_opts.sigma     = 1.0;

state_names  = {'x (V)', 'y (V)', 'z (V)'};
state_labels = {'x1', 'x2', 'x3'};

%% Pre-allocate storage for sweep and eigenvalue figure
X_datasets = cell(n_regimes, 1);
fig_eig    = figure('Name', 'Chua: Koopman Eigenvalues by Regime');

%% Process each regime
for r = 1:n_regimes
    csv_file   = regimes{r, 1};
    regime_name = regimes{r, 2};

    fprintf('\n=== Regime %d: %s ===\n', r, regime_name);

    %% Load CSV (skip leading comment lines starting with #)
    raw  = load_chua_csv(csv_file);
    t_ms = raw(:, 1);
    X_data = raw(:, 2:4)';      % (3 x N)  rows: x_V, y_V, z_V
    dt   = mean(diff(t_ms)) * 1e-3;   % convert ms -> s

    fprintf('  Samples: %d,  dt = %.4e s\n', size(X_data, 2), dt);

    %% Snapshot pairs
    X_all = X_data(:, 1:end-1);
    Y_all = X_data(:, 2:end);

    %% Lift and run EDMD
    [Psi_X, labels] = build_dictionary(X_all, dict_opts);
    [Psi_Y, ~]      = build_dictionary(Y_all, dict_opts);

    fprintf('  Dictionary: %d observables\n', numel(labels));

    [K, lambda, ~] = edmd(Psi_X, Psi_Y);

    %% One-step prediction error
    Psi_Y_pred = K * Psi_X;
    rel_err = norm(Psi_Y - Psi_Y_pred, 'fro') / norm(Psi_Y, 'fro');
    fprintf('  Relative one-step error (lifted): %.4e\n', rel_err);

    %% Top eigenvalues
    [~, ord] = sort(abs(lambda), 'descend');
    fprintf('  Top 5 Koopman eigenvalues (by magnitude):\n');
    disp(lambda(ord(1:min(5, end))));

    %% Eigenvalue subplot
    ax_eig = subplot(2, 2, r, 'Parent', fig_eig);
    plot_eigenvalues(ax_eig, lambda, regime_name);

    %% State indices in the dictionary
    state_idx = cellfun(@(lbl) find(strcmp(labels, lbl)), state_labels);

    %% Time axis for snapshot pairs
    M  = size(X_all, 2);
    t  = t_ms(1:M) * 1e-3;   % s

    %% Time-series plot
    figure('Name', sprintf('Regime %d: %s — time series', r, regime_name));
    for j = 1:3
        subplot(3, 1, j);
        plot(t, X_all(j, :), 'b-', t, Psi_Y_pred(state_idx(j), :), 'r--', ...
             'LineWidth', 1.2);
        xlabel('t (s)'); ylabel(state_names{j});
        legend('Measured', 'EDMD one-step', 'Location', 'best');
        if j == 1, title(sprintf('%s: one-step prediction', regime_name)); end
    end

    %% Phase portrait (x vs y)
    figure('Name', sprintf('Regime %d: %s — phase portrait', r, regime_name));
    plot(X_all(1, :), X_all(2, :), 'b-', 'LineWidth', 1.0);
    xlabel('x (V)'); ylabel('y (V)');
    title(sprintf('%s: phase portrait (measured)', regime_name));

    %% Save data for dictionary sweep
    X_datasets{r} = X_data;
end

% =========================================================================
%% Dictionary Sweep — compare types and sizes across all 4 regimes

sweep_cfgs = {
    struct('type','linear'),
    struct('type','pwl'),
    struct('type','poly',  'degree',1),
    struct('type','poly',  'degree',2),
    struct('type','rbf',   'n_centers',5,  'sigma',1.0),
    struct('type','poly',  'degree',3),
    struct('type','rbf',   'n_centers',10, 'sigma',1.0),
    struct('type','poly',  'degree',4),
    struct('type','rbf',   'n_centers',20, 'sigma',1.0),
    struct('type','poly',  'degree',5),
    struct('type','rbf',   'n_centers',50, 'sigma',1.0),
    % full: poly + kink only (no rbf) — isolates physics-informed terms
    struct('type','full',  'degree',2, 'n_centers',20, 'sigma',1.0, 'use_rbf',false),
    struct('type','full',  'degree',3, 'n_centers',20, 'sigma',1.0, 'use_rbf',false),
    struct('type','full',  'degree',4, 'n_centers',20, 'sigma',1.0, 'use_rbf',false),
    % full: poly + rbf only (no kink) — data-driven without Chua structure
    struct('type','full',  'degree',2, 'n_centers',10, 'sigma',1.0, 'use_pwl',false),
    struct('type','full',  'degree',3, 'n_centers',10, 'sigma',1.0, 'use_pwl',false),
    struct('type','full',  'degree',4, 'n_centers',10, 'sigma',1.0, 'use_pwl',false),
    % full: all three components at increasing richness
    struct('type','full',  'degree',2, 'n_centers',5,  'sigma',1.0),
    struct('type','full',  'degree',3, 'n_centers',10, 'sigma',1.0),
    struct('type','full',  'degree',3, 'n_centers',20, 'sigma',1.0),
    struct('type','full',  'degree',4, 'n_centers',20, 'sigma',1.0),
};
N_cfgs = numel(sweep_cfgs);

sweep_Npsi = zeros(n_regimes, N_cfgs);
sweep_err  = zeros(n_regimes, N_cfgs);
sweep_time = zeros(n_regimes, N_cfgs);   % milliseconds

fprintf('\n=== Dictionary Sweep ===\n');
for r = 1:n_regimes
    Xd = X_datasets{r};
    Xs = Xd(:, 1:end-1);
    Ys = Xd(:, 2:end);

    for ci = 1:N_cfgs
        tic;
        [Psi_X_sw, ~] = build_dictionary(Xs, sweep_cfgs{ci});
        [Psi_Y_sw, ~] = build_dictionary(Ys, sweep_cfgs{ci});
        [K_sw, ~, ~]  = edmd(Psi_X_sw, Psi_Y_sw);
        elapsed = toc * 1e3;

        Psi_Y_pred_sw = K_sw * Psi_X_sw;
        sweep_Npsi(r, ci) = size(Psi_X_sw, 1);
        sweep_err(r,  ci) = norm(Psi_Y_sw - Psi_Y_pred_sw, 'fro') / ...
                            norm(Psi_Y_sw, 'fro');
        sweep_time(r, ci) = elapsed;
    end
    fprintf('  Regime %d (%s): done\n', r, regimes{r, 2});
end

%% Sweep figure: 4x2 subplots — left col: accuracy, right col: timing
type_colors = containers.Map( ...
    {'linear', 'pwl', 'poly', 'rbf', 'full'}, ...
    {'k',      'g',   'b',    'm',   'r'});

fig_sw = figure('Name', 'Chua: Dictionary Sweep — Error & Timing');
for r = 1:n_regimes
    ax_err  = subplot(4, 2, 2*r - 1, 'Parent', fig_sw);
    ax_time = subplot(4, 2, 2*r,     'Parent', fig_sw);
    hold(ax_err,  'on');
    hold(ax_time, 'on');

    % Plot each configuration as a marker
    for ci = 1:N_cfgs
        col = type_colors(sweep_cfgs{ci}.type);
        np  = sweep_Npsi(r, ci);
        semilogy(ax_err,  np, sweep_err(r, ci),  [col 'o'], ...
            'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(ax_time, np, sweep_time(r, ci), [col 's'], ...
            'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end

    % Connect poly, rbf, and full families with lines (sorted by N_psi)
    for fam = {'poly', 'rbf', 'full'}
        idx = find(cellfun(@(c) strcmp(c.type, fam{1}), sweep_cfgs));
        col = type_colors(fam{1});
        [npsi_sorted, so] = sort(sweep_Npsi(r, idx));
        semilogy(ax_err,  npsi_sorted, sweep_err(r, idx(so)),  [col '-'], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
        plot(ax_time, npsi_sorted, sweep_time(r, idx(so)), [col '--'], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
    end

    ylabel(ax_err,  'Relative error (log scale)');
    xlabel(ax_err,  'N_{\psi}  (observables)');
    title(ax_err,   sprintf('%s: error', regimes{r, 2}));
    grid(ax_err,  'on');

    ylabel(ax_time, 'Computation time (ms)');
    xlabel(ax_time, 'N_{\psi}  (observables)');
    title(ax_time,  sprintf('%s: timing', regimes{r, 2}));
    grid(ax_time, 'on');
end

% Legend in the first error subplot using dummy lines
ax1 = subplot(4, 2, 1, 'Parent', fig_sw);
h = gobjects(5, 1);
fam_names  = {'linear', 'pwl', 'poly', 'rbf', 'full'};
for fi = 1:5
    col   = type_colors(fam_names{fi});
    h(fi) = plot(ax1, NaN, NaN, [col 'o-'], 'MarkerSize', 6, 'LineWidth', 1.2);
end
legend(ax1, h, fam_names, 'Location', 'best', 'FontSize', 7);

% =========================================================================
%% Double Scroll Extended Sweep — richer dictionaries for regime 4 only

ds_idx = 4;   % Double Scroll is regime 4
Xds = X_datasets{ds_idx};
Xs_ds = Xds(:, 1:end-1);
Ys_ds = Xds(:, 2:end);

sweep_cfgs_ds = {
    struct('type','rbf',  'n_centers',100, 'sigma',1.0),
    struct('type','rbf',  'n_centers',200, 'sigma',1.0),
    struct('type','rbf',  'n_centers',20,  'sigma',0.5),
    struct('type','rbf',  'n_centers',20,  'sigma',2.0),
    struct('type','poly', 'degree',6),
    struct('type','poly', 'degree',7),
    struct('type','full', 'degree',4, 'n_centers',50,  'sigma',1.0),
    struct('type','full', 'degree',3, 'n_centers',50,  'sigma',0.5),
};
N_cfgs_ds = numel(sweep_cfgs_ds);

ds_Npsi = zeros(1, N_cfgs_ds);
ds_err  = zeros(1, N_cfgs_ds);
ds_time = zeros(1, N_cfgs_ds);

fprintf('\n=== Double Scroll Extended Sweep ===\n');
for ci = 1:N_cfgs_ds
    tic;
    [Psi_X_ds, ~] = build_dictionary(Xs_ds, sweep_cfgs_ds{ci});
    [Psi_Y_ds, ~] = build_dictionary(Ys_ds, sweep_cfgs_ds{ci});
    [K_ds, ~, ~]  = edmd(Psi_X_ds, Psi_Y_ds);
    elapsed = toc * 1e3;

    Psi_Y_pred_ds = K_ds * Psi_X_ds;
    ds_Npsi(ci) = size(Psi_X_ds, 1);
    ds_err(ci)  = norm(Psi_Y_ds - Psi_Y_pred_ds, 'fro') / norm(Psi_Y_ds, 'fro');
    ds_time(ci) = elapsed;
    fprintf('  Config %d (%s): N_psi=%d, err=%.4e, t=%.1f ms\n', ...
        ci, sweep_cfgs_ds{ci}.type, ds_Npsi(ci), ds_err(ci), ds_time(ci));
end

function data = load_chua_csv(filepath)
% Read a Chua CSV file, skipping any lines that start with '#'.
fid  = fopen(filepath, 'r');
if fid == -1
    error('Cannot open file: %s', filepath);
end

lines = {};
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if ischar(line) && ~isempty(line) && line(1) ~= '#'
        lines{end+1} = line; %#ok<AGROW>
    end
end
fclose(fid);

% First remaining line is the header (time_ms,x_V,y_V,z_V) — skip it
lines = lines(2:end);

data = cellfun(@(l) sscanf(l, '%f,%f,%f,%f')', lines, 'UniformOutput', false);
data = vertcat(data{:});   % (N x 4)
end
