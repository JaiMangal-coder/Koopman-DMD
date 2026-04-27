% chua_lyapunov_spectrum.m
% Full 3-exponent Lyapunov spectrum for all Chua regimes.
%
% Method: Sano & Sawada (1985) local linear maps + QR orthogonalisation.
%   At each step t, a 3×3 Jacobian is estimated from k_nn nearest neighbours,
%   an orthonormal basis Q is propagated and re-orthogonalised via QR, and
%   log|diag(R)| accumulates the growth rates for all three directions.
%
% Time normalisation for comparison with ODE parameters (α=9, β=14.286):
%   Chua ODE uses dimensionless time  t' = t / (R_eff × C2)
%   With C2 = 100 nF and R_eff ≈ 1 kΩ at chaos onset → τ = 0.1 ms
%   λ_dim = λ_per_sample × (τ / dt_s)
%   Reference: λ₁_dim ≈ +0.23 for double-scroll (Matsumoto 1984)
%
clear; close all; clc;

% ── Circuit time normalisation ────────────────────────────────────────────────
C2_F  = 100e-9;   % C₂ = 100 nF  (large capacitor)
R_eff = 1000;     % Ω — pot position at double-scroll onset (~1 kΩ)
tau_s = R_eff * C2_F;          % 0.1 ms per dimensionless time unit

fprintf('τ = R_eff × C₂ = %.3f ms  →  1 dim-less unit = τ/dt\n', tau_s*1e3);
fprintf('Reference (α=9, β=14.286):  λ₁ ≈ +0.23, λ₂ ≈ 0, λ₃ ≪ 0\n\n');

% ── Load data ────────────────────────────────────────────────────────────────
data_dir  = '../chua';
files     = {'chua_regime1_fixed_point.csv', 'chua_regime2_limit_cycle.csv', ...
             'chua_regime3_period_doubled.csv', 'chua_regime4_double_scroll.csv'};
labels    = {'Fixed Point', 'Limit Cycle', 'Period-Doubled', 'Double Scroll'};
n_reg     = numel(files);
colors    = lines(n_reg);

xyz_all = cell(n_reg,1);
dt_s    = zeros(n_reg,1);
for i = 1:n_reg
    raw       = readmatrix(fullfile(data_dir, files{i}), 'NumHeaderLines', 2);
    dt_s(i)   = mean(diff(raw(:,1))) * 1e-3;   % ms → s
    xyz_all{i} = raw(:,2:4);                   % [x_V, y_V, z_V]
end

% ── Algorithm settings ───────────────────────────────────────────────────────
k_nn = 15;   % nearest neighbours for local Jacobian (need ≥ d+1 = 4)

% ── Compute spectrum for each regime ─────────────────────────────────────────
all_lam = zeros(n_reg, 3);   % final exponents (per sample), columns ordered by QR
all_lt  = cell(n_reg, 1);    % running averages (N×3)

fprintf('=== FULL LYAPUNOV SPECTRUM  (Sano-Sawada / QR) ===\n\n');

for reg = 1:n_reg
    xyz = xyz_all{reg};
    N   = size(xyz, 1);
    d   = 3;

    % ── Theiler window: 2 × first autocorrelation zero crossing ──────────────
    x0  = xyz(:,1) - mean(xyz(:,1));
    acf = xcorr(x0, min(N-1,400), 'normalized');
    acf = acf(ceil(end/2):end);
    zc  = find(acf(1:end-1)>0 & acf(2:end)<=0, 1);
    if isempty(zc), W = 30; else, W = max(2*zc, 20); end

    % ── QR / Gram-Schmidt evolution ───────────────────────────────────────────
    Q      = eye(d);             % initial orthonormal perturbation basis
    ly_sum = zeros(d, 1);        % cumulative log growth per direction
    cnt    = 0;
    lt     = zeros(N, d);        % running average at each step

    for t = 1:(N-1)
        xi = xyz(t, :);

        % Squared distances to all points (blank Theiler window + last point)
        d2 = sum((xyz - xi).^2, 2);
        d2(max(1,t-W) : min(N-1,t+W)) = inf;
        d2(N) = inf;

        [~, sidx] = sort(d2);
        nn = sidx(1 : min(k_nn*2, N));     % candidate pool
        nn = nn(d2(nn) < inf);              % finite distances only
        nn = nn(nn < N);                    % need successor x_{nn+1}
        nn = nn(1 : min(k_nn, numel(nn)));  % trim to k_nn

        if numel(nn) < d+1
            if cnt > 0, lt(t,:) = (ly_sum/cnt)'; end
            continue
        end

        % Centred input/output deviations
        dX = xyz(nn, :)   - xi;               % k×3
        dY = xyz(nn+1, :) - xyz(t+1, :);      % k×3

        % Local Jacobian via least squares: dY ≈ dX × A'  →  A = Jf (3×3)
        A = ((dX'*dX + 1e-8*eye(d)) \ (dX'*dY))';

        % Propagate orthonormal basis then QR re-orthogonalise
        Z      = A * Q;
        [Q, R] = qr(Z, 0);

        % Force positive diagonal in R so log is well-defined (Q stays ortho)
        sg     = sign(diag(R));  sg(sg==0) = 1;
        Q      = Q * diag(sg);
        R      = diag(sg) * R;

        log_diag = log(abs(diag(R)));
        log_diag(~isfinite(log_diag)) = -20;   % clamp -Inf (near-degenerate map)
        ly_sum = ly_sum + log_diag;
        cnt    = cnt + 1;
        lt(t,:) = (ly_sum / cnt)';      % running per-sample average
    end

    all_lam(reg,:) = (ly_sum / cnt)';
    all_lt{reg}    = lt;

    lam_samp = all_lam(reg,:);
    lam_dim  = lam_samp / dt_s(reg) * tau_s;
    DKY      = kaplan_yorke(sort(lam_samp, 'descend'));

    fprintf('%s  (W_Theiler = %d)\n', labels{reg}, W);
    fprintf('  per sample : [%+.5f,  %+.5f,  %+.5f]   sum = %+.5f\n', ...
            lam_samp(1), lam_samp(2), lam_samp(3), sum(lam_samp));
    fprintf('  dim-less   : [%+.4f,  %+.4f,  %+.4f]   D_KY = %.4f\n\n', ...
            lam_dim(1), lam_dim(2), lam_dim(3), DKY);
end

% ─────────────────────────────────────────────────────────────────────────────
%  FIGURE 1 — Convergence of λ₁, λ₂, λ₃  +  phase portraits
% ─────────────────────────────────────────────────────────────────────────────
lc = {[0.05 0.55 0.05], [0.45 0.05 0.75], [0.85 0.30 0.05]};  % green / purple / orange

fig1 = figure('Name','Lyapunov Spectrum Convergence','Position',[30 30 1400 700]);

for reg = 1:n_reg
    lt   = all_lt{reg};
    lf   = all_lam(reg,:);
    tv   = find(any(lt ~= 0, 2));   % steps where estimate is available

    % ── Convergence curves ────────────────────────────────────────────────────
    ax = subplot(2, n_reg, reg);
    hold on;
    for k = 1:3
        plot(tv, lt(tv,k), '-', 'Color', lc{k}, 'LineWidth', 1.3);
    end
    yline(0, 'k:', 'LineWidth', 1.0);
    xlabel('Step t');
    ylabel('\lambda_k  /sample');
    title(labels{reg}, 'FontSize', 9, 'FontWeight', 'bold');
    legend(sprintf('\\lambda_1 = %+.4f', lf(1)), ...
           sprintf('\\lambda_2 = %+.4f', lf(2)), ...
           sprintf('\\lambda_3 = %+.4f', lf(3)), ...
           'Location', 'northeast', 'FontSize', 7);
    grid on;

    % ── Phase portrait ────────────────────────────────────────────────────────
    xyz = xyz_all{reg};
    subplot(2, n_reg, reg + n_reg);
    plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.', ...
          'Color', colors(reg,:), 'MarkerSize', 2.5);
    xlabel('x (V)');  ylabel('y (V)');  zlabel('z (V)');
    title(labels{reg}, 'FontSize', 9);
    grid on;  view(45, 30);
end

sgtitle('Chua Circuit — Lyapunov Spectrum Convergence  (\lambda_1, \lambda_2, \lambda_3)', ...
        'FontSize', 12, 'FontWeight', 'bold');
saveas(fig1, 'chua_lyapunov_convergence.png');

% ─────────────────────────────────────────────────────────────────────────────
%  FIGURE 2 — Grouped bar chart: all regimes, all exponents (per sample)
% ─────────────────────────────────────────────────────────────────────────────
fig2 = figure('Name','Lyapunov Spectrum — All Regimes','Position',[30 500 1050 430]);

b = bar(1:n_reg, all_lam, 'grouped');
b(1).FaceColor = lc{1};
b(2).FaceColor = lc{2};
b(3).FaceColor = lc{3};

yline(0, 'k-', 'LineWidth', 1.5);
set(gca, 'XTick', 1:n_reg, 'XTickLabel', labels, 'FontSize', 9);
ylabel('\lambda_i  (per sample)');
title('Full Lyapunov Spectrum — All Chua Regimes', 'FontWeight', 'bold');
legend('\lambda_1  (max / unstable)', '\lambda_2  (neutral)', ...
       '\lambda_3  (min / stable)', 'Location', 'southwest');
grid on;
saveas(fig2, 'chua_lyapunov_bar.png');

% ─────────────────────────────────────────────────────────────────────────────
%  FIGURE 3 — Double Scroll spectrum in dimensionless ODE-time units
% ─────────────────────────────────────────────────────────────────────────────
lam4_dim = all_lam(4,:) / dt_s(4) * tau_s;      % dimensionless
lam4_sorted = sort(lam4_dim, 'descend');

fig3 = figure('Name','Double Scroll — Dimensionless Spectrum','Position',[30 720 680 440]);

bh = bar([1 2 3], lam4_dim, 'FaceColor', 'flat');
bh.CData = [0.10 0.40 0.80;   % λ₁  blue
            0.55 0.55 0.55;   % λ₂  grey
            0.80 0.20 0.20];  % λ₃  red
hold on;
yline(0, 'k-', 'LineWidth', 1.5);

% Reference line for theoretical λ₁ = 0.23
yline(0.23, '--', 'Color', [0 0.6 0], 'LineWidth', 2);
text(1.65, 0.26, 'theory: \lambda_1^{ODE} = 0.23', ...
     'Color', [0 0.6 0], 'FontSize', 10, 'FontWeight', 'bold');

set(gca, 'XTick', 1:3, 'XTickLabel', {'\lambda_1  (unstable)', ...
                                        '\lambda_2  (neutral)', ...
                                        '\lambda_3  (stable)'}, 'FontSize', 9);
ylabel('\lambda_i  (dimensionless,  ODE time units)');
title(sprintf(['Double Scroll — Lyapunov Spectrum\n' ...
               '\\tau = R_{eff} \\times C_2 = %.2f ms\n' ...
               '[%+.3f,   %+.3f,   %+.3f]  (dim-less)'], ...
               tau_s*1e3, lam4_dim(1), lam4_dim(2), lam4_dim(3)), ...
      'FontWeight', 'bold');
grid on;
saveas(fig3, 'chua_doublescroll_spectrum.png');

% ─────────────────────────────────────────────────────────────────────────────
%  SUMMARY TABLE
% ─────────────────────────────────────────────────────────────────────────────
fprintf('╔══════════════════╦═══════════════╦═══════════════╦═══════════════╦═════════╦════════╦══════════════════════════╗\n');
fprintf('║ Regime           ║ λ₁ /sample    ║ λ₂ /sample    ║ λ₃ /sample    ║   Σλ    ║ D_KY   ║ λ_dim [λ₁, λ₂, λ₃]      ║\n');
fprintf('╠══════════════════╬═══════════════╬═══════════════╬═══════════════╬═════════╬════════╬══════════════════════════╣\n');
for reg = 1:n_reg
    L     = sort(all_lam(reg,:), 'descend');
    DKY   = kaplan_yorke(L);
    Ldim  = L / dt_s(reg) * tau_s;
    fprintf('║ %-16s ║ %+11.5f   ║ %+11.5f   ║ %+11.5f   ║ %+7.4f ║ %6.3f ║ [%+.3f, %+.3f, %+.4f]  ║\n', ...
            labels{reg}, L(1), L(2), L(3), sum(L), DKY, Ldim(1), Ldim(2), Ldim(3));
end
fprintf('╚══════════════════╩═══════════════╩═══════════════╩═══════════════╩═════════╩════════╩══════════════════════════╝\n\n');

fprintf('Unit conversion: λ_dim = λ_per_sample / dt_s × τ\n');
fprintf('  dt_s = %.4f ms,  τ = %.3f ms,  scale = τ/dt = %.2f\n', ...
        dt_s(4)*1e3, tau_s*1e3, tau_s/dt_s(4));
lam4_sorted = sort(all_lam(4,:), 'descend');
lam4_1_dim  = lam4_sorted(1) / dt_s(4) * tau_s;
fprintf('Double-scroll λ₁_dim = %.4f  (theory: 0.23,  error = %.1f%%)\n', ...
        lam4_1_dim, 100*abs(lam4_1_dim - 0.23)/0.23);
fprintf('\nKaplan-Yorke dimension:  D_KY = j + Σ_{i≤j} λᵢ / |λ_{j+1}|\n');
fprintf('Figures: chua_lyapunov_convergence.png, chua_lyapunov_bar.png, chua_doublescroll_spectrum.png\n');

% ─────────────────────────────────────────────────────────────────────────────
%  LOCAL FUNCTION
% ─────────────────────────────────────────────────────────────────────────────
function D = kaplan_yorke(L)
% L must be sorted descending.
L = L(:);
j = find(cumsum(L) >= 0, 1, 'last');
if isempty(j)
    D = 0;
elseif j == numel(L)
    D = numel(L);
else
    D = j + sum(L(1:j)) / abs(L(j+1));
end
end
