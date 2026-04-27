% Chua Circuit Analysis
% Box-Counting Dimension · Correlation Dimension (G-P) · Largest Lyapunov Exponent
% ==================================================================================
clear; close all; clc;

data_dir  = '../chua';
files     = {'chua_regime1_fixed_point.csv', 'chua_regime2_limit_cycle.csv', ...
             'chua_regime3_period_doubled.csv', 'chua_regime4_double_scroll.csv'};
labels    = {'Regime 1: Fixed Point', 'Regime 2: Limit Cycle', ...
             'Regime 3: Period-Doubled', 'Regime 4: Double Scroll'};
short_lbl = {'Fixed Point', 'Limit Cycle', 'Period-Doubled', 'Double Scroll'};
n_reg     = numel(files);
colors    = lines(n_reg);

%% ── Load Data ────────────────────────────────────────────────────────────────
xyz_all = cell(n_reg,1);
dt_s    = zeros(n_reg,1);

fprintf('Loading data...\n');
for i = 1:n_reg
    raw      = readmatrix(fullfile(data_dir, files{i}), 'NumHeaderLines', 2);
    t_ms     = raw(:,1);
    dt_s(i)  = mean(diff(t_ms)) * 1e-3;   % ms → s
    xyz_all{i} = raw(:,2:4);              % [x_V, y_V, z_V]
    fprintf('  %-42s %4d pts,  dt = %.4f ms\n', files{i}, size(raw,1), dt_s(i)*1e3);
end

%% ══════════════════════════════════════════════════════════════════════════════
%% 1.  BOX-COUNTING FRACTAL DIMENSION
%% ══════════════════════════════════════════════════════════════════════════════
% Idea: cover the normalised phase portrait with ε-boxes; count N(ε) occupied.
%       D_box = slope of log N vs log(1/ε).
% With N=1000 pts the reliable scaling region is ε ≫ N^{-1/D}, so we trim the
% 3 coarsest scales (power-law not yet established) and the 3 finest (data-
% limited saturation) before fitting.

fprintf('\n=== BOX-COUNTING FRACTAL DIMENSIONS ===\n');

num_scales = 14;
eps_vals   = 2.^(-(1:num_scales));
D_box      = zeros(n_reg,1);

fig1 = figure('Name','Box-Counting','Position',[50 50 1400 700]);

for reg = 1:n_reg
    xyz = xyz_all{reg};
    N   = size(xyz,1);

    % Normalise independently per axis → [0,1]^3 (preserves local geometry)
    mn = min(xyz);  mx = max(xyz);
    sp = mx - mn;   sp(sp == 0) = 1;
    xyz_n = (xyz - mn) ./ sp;

    N_boxes = zeros(1, num_scales);
    for k = 1:num_scales
        bin_idx    = floor(xyz_n / eps_vals(k));
        N_boxes(k) = size(unique(bin_idx,'rows'), 1);
    end

    log_e = log2(eps_vals);          % negative (ε < 1)
    log_N = log2(N_boxes);

    % Valid: box count is finite, > 1, and < N (not trivially saturated)
    valid = isfinite(log_N) & log_N > 0.5 & log_N < log2(N) - 0.5;
    vidx  = find(valid);

    % Trim 3 coarsest + 3 finest from the fit; need at least 2 pts left
    trim = 3;
    if numel(vidx) > 2*trim + 2
        fit_idx = vidx((trim+1) : (end-trim));
    elseif numel(vidx) >= 2
        fit_idx = vidx;
    else
        D_box(reg) = NaN;
        continue;
    end

    p          = polyfit(log_e(fit_idx), log_N(fit_idx), 1);
    D_box(reg) = -p(1);

    fprintf('  %-30s D_box = %.4f  (fit over %d scales)\n', ...
            labels{reg}, D_box(reg), numel(fit_idx));

    % log-log plot
    ax = subplot(2, n_reg, reg);
    loglog(eps_vals(valid), 2.^log_N(valid), 'o-', ...
           'Color', colors(reg,:), 'MarkerFaceColor', colors(reg,:), 'MarkerSize',5);
    hold on;
    loglog(eps_vals(fit_idx), 2.^polyval(p, log_e(fit_idx)), '--k', 'LineWidth',2);
    set(ax,'XDir','reverse');
    xlabel('\epsilon  (box size)');  ylabel('N(\epsilon)');
    title(sprintf('%s\nD_{box} = %.3f', short_lbl{reg}, D_box(reg)));
    legend('Data','Fit','Location','northwest');  grid on;

    % Phase portrait
    subplot(2, n_reg, reg + n_reg);
    plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.', ...
          'Color', colors(reg,:), 'MarkerSize', 3);
    xlabel('x (V)');  ylabel('y (V)');  zlabel('z (V)');
    title(short_lbl{reg});  grid on;  view(45,30);
end
sgtitle('Chua Circuit — Box-Counting Fractal Dimension','FontSize',13,'FontWeight','bold');
saveas(fig1, 'chua_boxcounting.png');

%% ══════════════════════════════════════════════════════════════════════════════
%% 2.  CORRELATION DIMENSION  (Grassberger & Procaccia 1983)
%% ══════════════════════════════════════════════════════════════════════════════
% C(r) = fraction of point-pairs closer than r.  For small r: C(r) ∝ r^D₂.
% Using pdist (upper-triangle pairwise distances) avoids the O(N²) loop.

fprintf('\n=== CORRELATION DIMENSION  (Grassberger-Procaccia) ===\n');

D_corr = zeros(n_reg,1);
fig2   = figure('Name','Correlation Dimension','Position',[50 400 1400 380]);

for reg = 1:n_reg
    xyz  = xyz_all{reg};
    N    = size(xyz,1);

    % All N*(N-1)/2 pairwise distances (fast, vectorised)
    D_vec = pdist(xyz);

    % Radius range: 0.5 % to 150 % of median pairwise distance
    r_med  = median(D_vec);
    r_vals = r_med * logspace(-2, 0.7, 50);

    C_r = arrayfun(@(r) 2*sum(D_vec < r)/(N*(N-1)), r_vals);

    % Fit the scaling region: C in (0.01, 0.90) to avoid log(0) and saturation
    in_range = C_r > 0.01 & C_r < 0.90;
    if sum(in_range) >= 4
        p         = polyfit(log(r_vals(in_range)), log(C_r(in_range)), 1);
        D_corr(reg) = p(1);
    else
        D_corr(reg) = NaN;
    end

    fprintf('  %-30s D_corr = %.4f\n', labels{reg}, D_corr(reg));

    subplot(1, n_reg, reg);
    loglog(r_vals, C_r, 'o-', 'Color', colors(reg,:), ...
           'MarkerFaceColor', colors(reg,:), 'MarkerSize', 4);
    hold on;
    if sum(in_range) >= 4
        loglog(r_vals(in_range), exp(polyval(p, log(r_vals(in_range)))), ...
               '--k', 'LineWidth', 2);
    end
    xlabel('r  (V)');  ylabel('C(r)');
    title(sprintf('%s\nD_{corr} = %.3f', short_lbl{reg}, D_corr(reg)));
    grid on;  legend('C(r)','Fit','Location','northwest');
end
sgtitle('Chua Circuit — Correlation Dimension (Grassberger-Procaccia)', ...
        'FontSize',13,'FontWeight','bold');
saveas(fig2, 'chua_correlation_dim.png');

%% ══════════════════════════════════════════════════════════════════════════════
%% 3.  LARGEST LYAPUNOV EXPONENT  (Rosenstein, Collins & De Luca 1993)
%% ══════════════════════════════════════════════════════════════════════════════
% For each reference point i, find nearest spatial neighbour j outside a
% Theiler window W (prevents temporal correlations from masquerading as
% spatial closeness).  Track divergence d(t) = ||x_{i+t} – x_{j+t}||.
% The average <ln d(t)/d_0> grows linearly for the chaotic regime;
% the slope equals λ₁ in units of per-sample.  Multiply by 1/dt to get 1/s.
%
% Theiler window W: set to 2× the first zero-crossing lag of the
% autocorrelation so that at least one orbital period is excluded.

fprintf('\n=== LARGEST LYAPUNOV EXPONENTS  (Rosenstein method) ===\n');

max_iter  = 120;   % steps to track divergence (must be < N/2)
lambda1   = zeros(n_reg,1);
lambda1_s = zeros(n_reg,1);

fig3 = figure('Name','Lyapunov Exponents','Position',[50 550 1400 420]);

for reg = 1:n_reg
    xyz = xyz_all{reg};
    N   = size(xyz,1);

    % ── Auto-estimate Theiler window from autocorrelation of x-channel ──────
    x0    = xyz(:,1) - mean(xyz(:,1));
    acf   = xcorr(x0, min(N-1,400), 'normalized');
    acf   = acf(ceil(end/2):end);   % keep lags 0, 1, 2, …
    zc    = find(acf(1:end-1) > 0 & acf(2:end) <= 0, 1);
    if isempty(zc), zc = 30; end
    W     = max(2*zc, 20);   % Theiler window ≥ one full half-period
    fprintf('  %-30s  W_Theiler = %3d samples\n', labels{reg}, W);

    % ── Rosenstein nearest-neighbour divergence ──────────────────────────────
    ref_idx   = (1 : N - max_iter)';
    sum_log_d = zeros(max_iter,1);
    cnt       = zeros(max_iter,1);

    for ii = 1:numel(ref_idx)
        i  = ref_idx(ii);
        xi = xyz(i,:);

        d2          = sum((xyz - xi).^2, 2);
        lo          = max(1, i-W);
        hi          = min(N, i+W);
        d2(lo:hi)   = inf;      % blank Theiler window

        [d0_sq, j]  = min(d2);
        if isinf(d0_sq) || d0_sq == 0, continue; end
        if j + max_iter > N,           continue; end

        d0 = sqrt(d0_sq);
        for t = 1:max_iter
            dt_dist = norm(xyz(i+t,:) - xyz(j+t,:));
            if dt_dist > 0
                sum_log_d(t) = sum_log_d(t) + log(dt_dist / d0);
                cnt(t)       = cnt(t) + 1;
            end
        end
    end

    valid_t   = cnt > 0;
    avg_log_d = sum_log_d(valid_t) ./ cnt(valid_t);
    t_vec     = find(valid_t);

    % Fit the initial linear growth region (first ~25 % of tracked steps)
    n_fit   = max(5, round(0.25 * numel(t_vec)));
    p_fit   = polyfit(t_vec(1:n_fit)', avg_log_d(1:n_fit)', 1);

    lambda1(reg)   = p_fit(1);
    lambda1_s(reg) = p_fit(1) / dt_s(reg);

    fprintf('  %-30s  λ₁ = %+.5f /sample   (%+.1f /s)\n', ...
            labels{reg}, lambda1(reg), lambda1_s(reg));

    subplot(1, n_reg, reg);
    plot(t_vec, avg_log_d, '-', 'Color', colors(reg,:), 'LineWidth',1.5);
    hold on;
    plot(t_vec(1:n_fit), polyval(p_fit, t_vec(1:n_fit)'), '--k', 'LineWidth',2);
    xlabel('Time steps k');
    ylabel('avg  ln( d_k / d_0 )');
    title(sprintf('%s\n\\lambda_1 = %+.4f /samp\n(%+.1f /s)', ...
          short_lbl{reg}, lambda1(reg), lambda1_s(reg)));
    grid on;
end
sgtitle('Chua Circuit — Largest Lyapunov Exponent (Rosenstein 1993)', ...
        'FontSize',13,'FontWeight','bold');
saveas(fig3, 'chua_lyapunov.png');

%% ══════════════════════════════════════════════════════════════════════════════
%% 4.  SUMMARY TABLE
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n');
fprintf('╔══════════════════════════╦═══════════╦═══════════╦══════════════════╦══════════════════╗\n');
fprintf('║ Regime                   ║  D_box    ║  D_corr   ║  λ₁ (/sample)   ║  λ₁ (/s)        ║\n');
fprintf('╠══════════════════════════╬═══════════╬═══════════╬══════════════════╬══════════════════╣\n');
for reg = 1:n_reg
    if     lambda1(reg) >  0.005,   tag = 'CHAOTIC ←';
    elseif lambda1(reg) > -0.005,   tag = 'periodic';
    else,                            tag = 'stable  ';
    end
    fprintf('║ %-24s ║  %6.3f   ║  %6.3f   ║  %+12.5f   ║  %+12.1f   ║  %s\n', ...
            short_lbl{reg}, D_box(reg), D_corr(reg), lambda1(reg), lambda1_s(reg), tag);
end
fprintf('╚══════════════════════════╩═══════════╩═══════════╩══════════════════╩══════════════════╝\n');
fprintf('\n');
fprintf('Interpretation guide:\n');
fprintf('  D_box / D_corr > 2  →  fractal / strange attractor\n');
fprintf('  D_box / D_corr ≈ 1  →  limit cycle (1-D curve)\n');
fprintf('  D_box / D_corr ≈ 0  →  fixed point (0-D)\n');
fprintf('  λ₁ > 0              →  chaos  (sensitive dependence on initial conditions)\n');
fprintf('  λ₁ ≈ 0              →  periodic orbit (neutral stability along flow)\n');
fprintf('  λ₁ < 0              →  stable attractor (all trajectories contract)\n');
fprintf('\nData-density note:\n');
fprintf('  With N=1000 pts the box/correlation dimensions of the double scroll\n');
fprintf('  are underestimated vs the known theoretical value D_F ≈ 2.05–2.3.\n');
fprintf('  Reliable box-counting requires ~10^D pts; a 50k-pt dataset would\n');
fprintf('  push D_box closer to 2.  The relative ordering across regimes is valid.\n');
fprintf('\nFigures saved: chua_boxcounting.png, chua_correlation_dim.png, chua_lyapunov.png\n');
