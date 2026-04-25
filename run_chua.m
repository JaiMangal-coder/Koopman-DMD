% RUN_CHUA  Apply EDMD to all four Chua circuit experimental regimes
%
% Data files (in chua/ subdirectory), each with columns: time_ms, x_V, y_V, z_V
%   chua_regime1_fixed_point.csv
%   chua_regime2_limit_cycle.csv
%   chua_regime3_period_doubled.csv
%   chua_regime4_double_scroll.csv

clear; clc; close all;

%% Regime definitions
regimes = { ...
    'chua/chua_regime1_fixed_point.csv',    'Fixed Point'; ...
    'chua/chua_regime2_limit_cycle.csv',    'Limit Cycle'; ...
    'chua/chua_regime3_period_doubled.csv', 'Period-Doubled'; ...
    'chua/chua_regime4_double_scroll.csv',  'Double Scroll'; ...
};
n_regimes = size(regimes, 1);

%% Dictionary options
% Use the same dictionary for cross-regime comparisons such as the
% eigenvalue figure, one-step error, and baseline predictions.
dict_opts.type = 'poly';

state_names  = {'x (V)', 'y (V)', 'z (V)'};
state_labels = {'x1', 'x2', 'x3'};

%% EDMD options
edmd_opts.method = 'tsvd';
edmd_opts.tol    = 1e-6;

%% Pre-allocate storage for sweep and eigenvalue figure
X_datasets = cell(n_regimes, 1);
fig_eig    = figure('Name', 'Chua: Koopman Eigenvalues by Regime');
eig_title_set = false;

%% Phase portrait figure (2x2, filled inside the loop)
fig_phase = figure('Name', 'Chua: Phase Portraits — Measured Signals', ...
                   'Position', [150 150 900 700]);
tl_phase  = tiledlayout(fig_phase, 2, 2, ...
                        'Padding', 'compact', 'TileSpacing', 'compact');
title(tl_phase, 'Chua Circuit Phase Portraits — Measured Signals', ...
      'FontSize', 13, 'FontWeight', 'bold');

%% Process each regime
for r = 1:n_regimes
    csv_file    = regimes{r, 1};
    regime_name = regimes{r, 2};

    fprintf('\n=== Regime %d: %s ===\n', r, regime_name);

    %% Load CSV (skip leading comment lines starting with #)
    raw  = load_chua_csv(csv_file);
    t_ms = raw(:, 1);
    X_data = raw(:, 2:4)';      % (3 x N) rows: x_V, y_V, z_V
    dt   = mean(diff(t_ms)) * 1e-3;   % convert ms -> s

    fprintf('  Samples: %d,  dt = %.4e s\n', size(X_data, 2), dt);

    %% Snapshot pairs
    X_all = X_data(:, 1:end-1);
    Y_all = X_data(:, 2:end);

    %% Lift and run EDMD
    [Psi_X, labels] = build_dictionary(X_all, dict_opts);
    [Psi_Y, ~]      = build_dictionary(Y_all, dict_opts);

    fprintf('  Dictionary: %d observables\n', numel(labels));

    [K, lambda, Phi] = edmd(Psi_X, Psi_Y, edmd_opts);

    %% One-step prediction error
    Psi_Y_pred = K * Psi_X;
    rel_err = norm(Psi_Y - Psi_Y_pred, 'fro') / norm(Psi_Y, 'fro');
    fprintf('  Relative one-step error (lifted): %.4e\n', rel_err);

    %% Top eigenvalues
    [~, ord] = sort(abs(lambda), 'descend');
    fprintf('  Top 5 Koopman eigenvalues (common dictionary, by magnitude):\n');
    disp(lambda(ord(1:min(5, end))));

    %% Eigenvalue subplot
    ax_eig = subplot(2, 2, r, 'Parent', fig_eig);
    plot_eigenvalues(ax_eig, lambda, regime_name);
    if ~eig_title_set
        sgtitle(fig_eig, sprintf('Chua: Koopman Eigenvalues by Regime (N_{\\psi} = %d)', ...
            numel(labels)));
        eig_title_set = true;
    end

    %% State indices in the dictionary
    state_idx = zeros(1, 3);
    for j = 1:3
        idx = find(strcmp(labels, state_labels{j}), 1);
        if isempty(idx)
            error('State label "%s" not found in dictionary.', state_labels{j});
        end
        state_idx(j) = idx;
    end

    %% Time axis for snapshot pairs
    M  = size(X_all, 2);
    t  = t_ms(1:M) * 1e-3;   % s

    %% Time-series plot: one-step reconstruction in state coordinates
    figure('Name', sprintf('Regime %d: %s — time series', r, regime_name));
    for j = 1:3
        subplot(3, 1, j);
        plot(t, X_all(j, :), 'b-', ...
             t, real(Psi_Y_pred(state_idx(j), :)), 'r--', ...
             'LineWidth', 1.2);
        xlabel('t (s)');
        ylabel(state_names{j});
        legend('Measured', 'EDMD one-step', 'Location', 'best');
        if j == 1
            title(sprintf('%s: one-step prediction', regime_name));
        end
        grid on;
    end

    %% Phase portrait (x vs y) — add tile to the shared 2×2 figure
    ax_phase = nexttile(tl_phase, r);
    plot(ax_phase, X_all(1, :), X_all(2, :), 'b-', 'LineWidth', 0.8);
    xlabel(ax_phase, 'x (V)', 'FontSize', 10);
    ylabel(ax_phase, 'y (V)', 'FontSize', 10);
    title(ax_phase, regime_name, 'FontSize', 11, 'FontWeight', 'bold');
    grid(ax_phase, 'on');
    box(ax_phase, 'on');
    set(ax_phase, 'FontSize', 10);

    %% Double-scroll eigenfunction coloring (r=4)
    % Use a richer Chua-matched dictionary, then evaluate the Koopman
    % eigenfunctions xi_k(x_n) = w_k^* psi(x_n).  Coloring the attractor by
    % arg(xi_k) or sign(Re(xi_k)) can reveal hidden coordinates that are not
    % obvious in the raw state-space trace.
    if r == 4
        fprintf('  Building eigenfunction-colored double-scroll plots...\n');

        dict_color = struct('type', 'full', 'degree', 3, ...
                            'use_rbf', false, 'use_pwl', true);
        [Psi_X_color, ~] = build_dictionary(X_all, dict_color);
        [Psi_Y_color, ~] = build_dictionary(Y_all, dict_color);
        [~, lambda_color, Phi_color] = edmd(Psi_X_color, Psi_Y_color, edmd_opts);

        W_color  = pinv(Phi_color);
        Xi_color = W_color * Psi_X_color;   % each row = one Koopman eigenfunction over time

        k_color = pick_visual_eigenfunction(lambda_color, Xi_color);
        xi_vis  = Xi_color(k_color, :);
        xi_phase = angle(xi_vis);
        xi_real  = real(xi_vis);
        xi_sign  = xi_real >= 0;

        fk = abs(angle(lambda_color(k_color))) / (2*pi*dt);
        fprintf('  Visual eigenfunction mode %d: lambda = %.4f%+.4fi, |lambda| = %.4f, f = %.3f Hz\n', ...
            k_color, real(lambda_color(k_color)), imag(lambda_color(k_color)), ...
            abs(lambda_color(k_color)), fk);
        fprintf('  Sign split: %d positive, %d negative samples\n', ...
            nnz(xi_sign), nnz(~xi_sign));

        fig_color = figure('Name', sprintf('%s — Eigenfunction Coloring', regime_name), ...
                           'Position', [120 80 1350 900]);
        tl_color  = tiledlayout(fig_color, 2, 2, ...
                                'Padding', 'compact', 'TileSpacing', 'compact');
        title(tl_color, sprintf(['%s: attractor colored by Koopman eigenfunction ', ...
            'phase/sign (mode %d, |\\lambda|=%.4f, f=%.3f Hz)'], ...
            regime_name, k_color, abs(lambda_color(k_color)), fk), ...
            'FontSize', 12, 'FontWeight', 'bold');

        ax1 = nexttile(tl_color);
        scatter(ax1, X_all(1, :), X_all(2, :), 8, xi_phase, 'filled');
        xlabel(ax1, 'x (V)');
        ylabel(ax1, 'y (V)');
        title(ax1, 'x-y projection colored by eigenfunction phase');
        colormap(ax1, hsv);
        cb1 = colorbar(ax1);
        cb1.Label.String = 'arg(\xi_k)';
        grid(ax1, 'on');

        ax2 = nexttile(tl_color);
        scatter(ax2, X_all(1, :), X_all(3, :), 8, xi_phase, 'filled');
        xlabel(ax2, 'x (V)');
        ylabel(ax2, 'z (V)');
        title(ax2, 'x-z projection colored by eigenfunction phase');
        colormap(ax2, hsv);
        cb2 = colorbar(ax2);
        cb2.Label.String = 'arg(\xi_k)';
        grid(ax2, 'on');

        ax3 = nexttile(tl_color);
        plot3(ax3, X_all(1, xi_sign),  X_all(2, xi_sign),  X_all(3, xi_sign), ...
              '.', 'Color', [0.15 0.35 0.85], 'MarkerSize', 5);
        hold(ax3, 'on');
        plot3(ax3, X_all(1, ~xi_sign), X_all(2, ~xi_sign), X_all(3, ~xi_sign), ...
              '.', 'Color', [0.85 0.25 0.20], 'MarkerSize', 5);
        xlabel(ax3, 'x (V)');
        ylabel(ax3, 'y (V)');
        zlabel(ax3, 'z (V)');
        title(ax3, '3D attractor split by sign(Re(\xi_k))');
        legend(ax3, 'Re(\xi_k) > 0', 'Re(\xi_k) < 0', 'Location', 'best');
        grid(ax3, 'on');
        view(ax3, 35, 25);

        ax4 = nexttile(tl_color);
        plot(ax4, t, xi_real, 'b-', 'LineWidth', 0.9);
        hold(ax4, 'on');
        plot(ax4, t, imag(xi_vis), 'r-', 'LineWidth', 0.9);
        yline(ax4, 0, 'k:');
        xlabel(ax4, 't (s)');
        ylabel(ax4, '\xi_k');
        title(ax4, 'Eigenfunction time trace along the trajectory');
        legend(ax4, 'Re(\xi_k)', 'Im(\xi_k)', 'Location', 'best');
        grid(ax4, 'on');

        exportgraphics(fig_color, 'writeup/doubleScrollColoring.png', 'Resolution', 150);
    end

    %% Koopman mode reconstruction — limit cycle (r=2) and period-doubled (r=3)
    %
    % MATHEMATICAL BASIS
    % ------------------
    % EDMD finds K such that K*Psi(x_n) ≈ Psi(x_{n+1}).  Its eigendecomposition
    %   K * Phi = Phi * diag(lambda)
    % gives the spectral expansion
    %   Psi(x_n) ≈ K^n * psi0 = Phi * diag(lambda)^n * a,   a = pinv(Phi)*psi0
    %
    % Projecting onto the physical state with the selector C (3 x N_psi):
    %   x_n ≈ C * Psi(x_n) = sum_k  a_k * lambda_k^n * v_k
    %
    % where v_k = C * phi_k in R^3 is the k-th Koopman mode — a 3-vector whose
    % three entries give the relative amplitude and phase with which each of the
    % three state variables (x, y, z) participates in mode k.  The scalar
    % a_k * lambda_k^n drives the temporal evolution.
    %
    % For the limit cycle the dominant eigenvalues come in conjugate pairs at
    % the fundamental frequency f0 and its harmonics; including all bounded
    % modes (|lambda| <= 1.05) reconstructs the full waveform.
    %
    % For the period-doubled attractor the same formula holds, but the leading
    % pair sits at f0/2, and pairs at 3f0/2, 5f0/2, ... also appear — a direct
    % signature of the period-doubling bifurcation visible in all three states.
    if r == 2 || r == 3
        fprintf('  Building Koopman mode reconstruction for %s...\n', regime_name);

        % For period-doubled: the pwl dictionary (7 observables) cannot
        % represent the period-doubling Koopman eigenfunction — the function
        % that changes sign between the two half-cycles.  Re-run EDMD with a
        % degree-3 polynomial dictionary (20 observables) so the half-period
        % mode at ~f0/2 can appear.  The main-loop one-step prediction is
        % unchanged; only the reconstruction uses this richer decomposition.
        if r == 3
            dict_recon = struct('type', 'poly', 'degree', 3);
            [Psi_X_recon, labels_recon] = build_dictionary(X_all, dict_recon);
            [Psi_Y_recon, ~]            = build_dictionary(Y_all, dict_recon);
            [K_recon, lambda, Phi] = edmd(Psi_X_recon, Psi_Y_recon, edmd_opts);
            Psi_X = Psi_X_recon;
            err_recon = norm(Psi_Y_recon - K_recon*Psi_X_recon, 'fro') / norm(Psi_Y_recon, 'fro');
            fprintf('  Poly-3 dictionary (%d observables), one-step error: %.4e\n', ...
                size(Psi_X_recon, 1), err_recon);
            % Recompute state indices in the new dictionary
            for j = 1:3
                idx = find(strcmp(labels_recon, state_labels{j}), 1);
                if isempty(idx)
                    error('State label "%s" not found in poly-3 dictionary.', state_labels{j});
                end
                state_idx(j) = idx;
            end
        end

        % Modal coordinates: a = pinv(Phi) * psi0
        W    = pinv(Phi);
        psi0 = Psi_X(:, 1);
        a    = W * psi0;

        % Selector matrix C (3 x N_psi): extracts physical states from Psi
        Npsi = size(Psi_X, 1);
        C    = zeros(3, Npsi);
        for j = 1:3
            C(j, state_idx(j)) = 1;
        end

        % Koopman modes in physical coordinates: V(:,k) = C * phi_k in R^3
        Vkoop = C * Phi;

        % Retain only bounded/near-bounded modes
        stable_mask = abs(lambda) <= 1.05;
        fprintf('  Number of bounded modes |lambda| <= 1.05: %d / %d\n', ...
            nnz(stable_mask), numel(lambda));

        % state_modes(j, n, k) = Re[ v_k(j) * a_k * lambda_k^{n-1} ]
        % Each slice (:,:,k) is the 3-state time series contributed by mode k.
        nvec        = 0:(M-1);
        state_modes = zeros(3, M, Npsi);
        for k = 1:Npsi
            if stable_mask(k)
                time_factor          = a(k) * (lambda(k) .^ nvec);   % 1 x M
                state_modes(:, :, k) = Vkoop(:, k) * time_factor;    % 3 x M
            end
        end

        % Unit-circle projected modes
        % ----------------------------
        % On any periodic attractor the true Koopman eigenvalues for sustained
        % oscillation lie exactly on the unit circle (|lambda| = 1).  EDMD
        % approximates the global operator with a finite dictionary, so
        % eigenvalues are pulled slightly inside: |lambda| < 1.  Over thousands
        % of samples this causes a_k * lambda_k^n -> 0 even though the physical
        % signal never decays.
        %
        % Fix: for oscillatory modes (Im(lambda) != 0) replace lambda_k with
        % lambda_k / |lambda_k| in the time propagation only.  This preserves
        % the frequency and phase that EDMD identified while removing the
        % spurious amplitude decay.  Real/DC modes (Im ≈ 0) are left unchanged.
        state_modes_unit = zeros(3, M, Npsi);
        for k = 1:Npsi
            if stable_mask(k)
                lk = lambda(k);
                if abs(imag(lk)) > 1e-8
                    lk = lk / abs(lk);   % project to unit circle
                end
                time_factor              = a(k) * (lk .^ nvec);
                state_modes_unit(:,:,k)  = Vkoop(:, k) * time_factor;
            end
        end

        % Measured signals for all three states
        x_true = X_all;   % 3 x M  (rows: x, y, z)

        % --- All bounded modes ---
        x_rec_all      = real(sum(state_modes(     :, :, stable_mask), 3));
        x_rec_unit_all = real(sum(state_modes_unit(:, :, stable_mask), 3));

        % Lifted-space modal consistency (diagnostic)
        stable_idx = find(stable_mask);
        Psi_modal  = zeros(size(Psi_X));
        for n = 1:M
            tmp = zeros(Npsi, 1);
            for kk = stable_idx.'
                tmp = tmp + Phi(:, kk) * (a(kk) * lambda(kk)^(n-1));
            end
            Psi_modal(:, n) = tmp;
        end
        Psi_err = norm(Psi_X - Psi_modal, 'fro') / norm(Psi_X, 'fro');
        fprintf('  Bounded lifted-space modal consistency error: %.4e\n', Psi_err);

        % --- Rank modes by x-signal energy, preserving conjugate pairs ---
        % Score using the unit-projected modes so decaying amplitudes do not
        % bias the ranking (matters most for the period-doubled regime).
        mode_score = zeros(Npsi, 1);
        for k = 1:Npsi
            if stable_mask(k)
                mode_score(k) = norm(real(state_modes_unit(1, :, k)));
            end
        end
        [~, ord_modes] = sort(mode_score, 'descend');

        selected_3  = pick_conjugate_safe_modes(lambda, ord_modes, 3,  stable_mask);
        selected_10 = pick_conjugate_safe_modes(lambda, ord_modes, 10, stable_mask);

        % Reconstruct all 3 states from selected modes (both raw and projected)
        x_rec_3       = real(sum(state_modes(     :, :, selected_3),  3));
        x_rec_10      = real(sum(state_modes(     :, :, selected_10), 3));
        x_rec_unit_3  = real(sum(state_modes_unit(:, :, selected_3),  3));
        x_rec_unit_10 = real(sum(state_modes_unit(:, :, selected_10), 3));

        % Reconstruction errors per state (unit-projected, as the honest metric
        % for periodic attractors)
        state_sym = {'x', 'y', 'z'};
        fprintf('  Modal reconstruction errors per state (unit-circle projected):\n');
        for j = 1:3
            ea  = norm(x_true(j,:) - x_rec_unit_all(j,:)) / norm(x_true(j,:));
            e3  = norm(x_true(j,:) - x_rec_unit_3(j,:))   / norm(x_true(j,:));
            e10 = norm(x_true(j,:) - x_rec_unit_10(j,:))  / norm(x_true(j,:));
            fprintf('    %s:  all-bounded=%.3e  top-3=%.3e  top-10=%.3e\n', ...
                state_sym{j}, ea, e3, e10);
        end

        % Selected mode diagnostics
        fprintf('  Selected eigenvalues (top ~3 modes):\n');
        disp(lambda(selected_3));
        fprintf('  Selected eigenvalues (top ~10 modes):\n');
        disp(lambda(selected_10));

        % Frequencies of dominant oscillatory modes
        % For period-doubled (r=3) the lowest frequency should be ~f0/2
        % where f0 is the limit-cycle fundamental.  If all listed frequencies
        % are multiples of f0, the half-period mode is not being captured by
        % the dictionary and the reconstruction will not show two loops.
        fprintf('  Selected modal frequencies (top-10 bounded modes):\n');
        sel_freqs = [];
        for idx = selected_10
            if abs(imag(lambda(idx))) > 1e-8
                fk = angle(lambda(idx)) / (2*pi*dt);
                sel_freqs(end+1) = abs(fk); %#ok<AGROW>
                fprintf('    mode %d: |lambda|=%.4f, f=%.3f Hz\n', ...
                    idx, abs(lambda(idx)), fk);
            end
        end
        if r == 3 && ~isempty(sel_freqs)
            f_min = min(sel_freqs);
            f_max = max(sel_freqs);
            fprintf('  Period-doubled check: lowest f=%.3f Hz, highest f=%.3f Hz\n', ...
                f_min, f_max);
            fprintf('  Ratio highest/lowest = %.2f  (expect ~%d for %d harmonics)\n', ...
                f_max/f_min, round(f_max/f_min), round(f_max/f_min));
        end

        % ============================================================
        % Window length for time-series plots: show ~6 complete periods
        % so individual cycles are legible rather than a solid block.
        % Estimate period from the lowest-frequency oscillatory mode selected.
        osc_freqs = arrayfun(@(k) abs(angle(lambda(k)))/(2*pi*dt), selected_10);
        osc_freqs(osc_freqs < 1e-4) = inf;   % ignore DC
        f_dom  = min(osc_freqs);             % dominant (lowest) frequency
        n_plot = min(M, max(200, round(6 / (f_dom * dt))));
        t_win  = t(1:n_plot);
        fprintf('  Plot window: %d samples (%.1f periods of %.2f Hz)\n', ...
            n_plot, n_plot*dt*f_dom, f_dom);

        % ============================================================
        % Three separate figures (no uitab — easier to save):
        %   Figure 1 — Time Series   : 3 states × 3 mode sets (3×3 grid)
        %   Figure 2 — Modal Contrib : per-mode x/y/z contributions
        %   Figure 3 — Phase Plane   : reconstructed attractor shape
        rec_data   = {x_rec_unit_all, x_rec_unit_3,  x_rec_unit_10};
        rec_labels = {'All bounded (proj.)', 'Top ~3 (proj.)', 'Top ~10 (proj.)'};
        rec_colors = {'g--', 'r--', 'b--'};

        % ------ Figure 1: Time Series ------
        fig_ts = figure('Name', sprintf('%s — Time Series Reconstruction', regime_name), ...
                        'Position', [80 60 1300 900]);
        tl_ts  = tiledlayout(fig_ts, 3, 3, ...
                             'Padding', 'compact', 'TileSpacing', 'compact');
        title(tl_ts, ...
            sprintf('%s: Koopman reconstruction of x, y, z (unit-circle projected)', ...
                regime_name), 'FontSize', 12, 'FontWeight', 'bold');

        for j = 1:3
            for ci = 1:3
                ax = nexttile(tl_ts);
                plot(ax, t_win, x_true(j, 1:n_plot), 'k-', 'LineWidth', 1.0); hold(ax,'on');
                plot(ax, t_win, rec_data{ci}(j, 1:n_plot), rec_colors{ci}, 'LineWidth', 1.0);
                err_val = norm(x_true(j,:) - rec_data{ci}(j,:)) / norm(x_true(j,:));
                xlabel(ax, 't (s)');
                ylabel(ax, state_names{j});
                title(ax, sprintf('%s  |  %s  (err=%.2e)', ...
                    state_sym{j}, rec_labels{ci}, err_val));
                if j == 1 && ci == 1
                    legend(ax, 'Measured', 'Reconstruction', 'Location', 'best');
                end
                grid(ax, 'on');
            end
        end

        % ------ Figure 2: Modal Contributions ------
        nshow  = min(6, numel(selected_10));
        fig_mc = figure('Name', sprintf('%s — Modal Contributions', regime_name), ...
                        'Position', [80 60 1300 900]);
        tl_mc  = tiledlayout(fig_mc, nshow, 3, ...
                             'Padding', 'compact', 'TileSpacing', 'compact');
        title(tl_mc, ...
            sprintf('%s: per-mode contributions to x, y, z (unit-circle projected)', ...
                regime_name), 'FontSize', 12, 'FontWeight', 'bold');

        for ii = 1:nshow
            k = selected_10(ii);
            for j = 1:3
                ax = nexttile(tl_mc);
                plot(ax, t_win, real(state_modes_unit(j, 1:n_plot, k)), 'LineWidth', 1.0);
                ylabel(ax, state_names{j});
                if j == 1
                    title(ax, sprintf('Mode %d  \\lambda=%.3f%+.3fi  |\\lambda|=%.4f', ...
                        k, real(lambda(k)), imag(lambda(k)), abs(lambda(k))));
                end
                if ii == nshow
                    xlabel(ax, 't (s)');
                end
                grid(ax, 'on');
            end
        end

        % ------ Figure 3: Phase Plane ------
        % 3 projections (x-y, x-z, y-z) × 2 rows (measured | top-10 recon).
        %
        % For the period-doubled regime we use the f0/2 Koopman eigenfunction
        % as a loop label.  That eigenfunction is the unique observable that
        % changes sign exactly once per doubled period — positive on loop A,
        % negative on loop B.  Coloring measured and reconstructed samples by
        % that sign (blue = A, red = B) separates the two loops even when they
        % are geometrically close and the trajectory is undersampled.
        proj_pairs = [1 2; 1 3; 2 3];
        proj_xlbls = {'x (V)', 'x (V)', 'y (V)'};
        proj_ylbls = {'y (V)', 'z (V)', 'z (V)'};
        proj_names = {'x-y', 'x-z', 'y-z'};

        % Compute loop label for period-doubled only.  The f0/2 eigenfunction's
        % complex time factor a_k * lambda_proj^n has phase that increases
        % monotonically by angle(lambda_proj) each step, completing 2pi over
        % one doubled period.  Its sign labels the two loops.
        if r == 3
            [~, ii_half] = min(osc_freqs);
            k_half       = selected_10(ii_half);
            loop_sign = real(state_modes_unit(1, :, k_half));   % 1×M
            loop_A    = loop_sign >= 0;

            fprintf('  Loop label: %d loop-A, %d loop-B samples\n', ...
                nnz(loop_A), nnz(~loop_A));
        end

        fig_pp = figure('Name', sprintf('%s — Phase Plane Reconstruction', regime_name), ...
                        'Position', [80 60 1300 800]);
        tl_pp  = tiledlayout(fig_pp, 2, 3, ...
                             'Padding', 'compact', 'TileSpacing', 'compact');
        title(tl_pp, ...
            sprintf('%s: phase portraits — measured (top) vs top-10 reconstruction (bottom)', ...
                regime_name), 'FontSize', 11, 'FontWeight', 'bold');

        x_rec_plot = x_rec_unit_10;

        for row = 1:2
            for col = 1:3
                ax  = nexttile(tl_pp);
                pi_ = proj_pairs(col, :);
                dat = x_true;
                lbl = 'Measured';
                if row == 2
                    dat = x_rec_plot;
                    lbl = 'Top-10 recon.';
                end

                if r == 3
                    plot(ax, dat(pi_(1), loop_A),  dat(pi_(2), loop_A),  'b.', 'MarkerSize', 3);
                    hold(ax, 'on');
                    plot(ax, dat(pi_(1), ~loop_A), dat(pi_(2), ~loop_A), 'r.', 'MarkerSize', 3);
                    if row == 1 && col == 1
                        legend(ax, 'Loop A', 'Loop B', ...
                               'Location', 'northwest', 'FontSize', 8);
                    end
                else
                    plot(ax, dat(pi_(1),:), dat(pi_(2),:), 'k-', 'LineWidth', 0.4);
                end

                xlabel(ax, proj_xlbls{col});
                ylabel(ax, proj_ylbls{col});
                title(ax, sprintf('%s — %s', proj_names{col}, lbl));
                grid(ax, 'on');
            end
        end

        if r == 3
            exportgraphics(fig_pp, 'writeup/periodDoubledPP.png', 'Resolution', 150);
        end
    end

    %% Save data for dictionary sweep
    X_datasets{r} = X_data;
end

%% Export per-regime composite figures (after all regimes have filled them in)
exportgraphics(fig_eig,   'writeup/eigenvalues.png',  'Resolution', 150);
exportgraphics(fig_phase, 'writeup/MeasuredChua.png', 'Resolution', 150);

% =========================================================================
%% Dictionary Sweep — compare types and sizes across all 4 regimes

sweep_cfgs = {
    struct('type','linear'), ...
    struct('type','pwl'), ...
    struct('type','poly',  'degree',1), ...
    struct('type','poly',  'degree',2), ...
    struct('type','rbf',   'n_centers',5,  'sigma',1.0), ...
    struct('type','poly',  'degree',3), ...
    struct('type','rbf',   'n_centers',10, 'sigma',1.0), ...
    struct('type','poly',  'degree',4), ...
    struct('type','rbf',   'n_centers',20, 'sigma',1.0), ...
    struct('type','poly',  'degree',5), ...
    struct('type','rbf',   'n_centers',50, 'sigma',1.0), ...
    struct('type','full',  'degree',2, 'n_centers',20, 'sigma',1.0, 'use_rbf',false), ...
    struct('type','full',  'degree',3, 'n_centers',20, 'sigma',1.0, 'use_rbf',false), ...
    struct('type','full',  'degree',4, 'n_centers',20, 'sigma',1.0, 'use_rbf',false), ...
    struct('type','full',  'degree',2, 'n_centers',10, 'sigma',1.0, 'use_pwl',false), ...
    struct('type','full',  'degree',3, 'n_centers',10, 'sigma',1.0, 'use_pwl',false), ...
    struct('type','full',  'degree',4, 'n_centers',10, 'sigma',1.0, 'use_pwl',false), ...
    struct('type','full',  'degree',2, 'n_centers',5,  'sigma',1.0), ...
    struct('type','full',  'degree',3, 'n_centers',10, 'sigma',1.0), ...
    struct('type','full',  'degree',3, 'n_centers',20, 'sigma',1.0), ...
    struct('type','full',  'degree',4, 'n_centers',20, 'sigma',1.0), ...
};
N_cfgs = numel(sweep_cfgs);

sweep_Npsi = zeros(n_regimes, N_cfgs);
sweep_err  = zeros(n_regimes, N_cfgs);
sweep_time = zeros(n_regimes, N_cfgs);   % milliseconds

fprintf('\n=== Dictionary Sweep ===\n');
for r = 1:n_regimes
    Xd = X_datasets{r};

    % Normalize each state to [-1, 1] so high-degree monomials stay O(1).
    % Use full-dataset min/max so Xs and Ys share the same scaling.
    xmin  = min(Xd, [], 2);
    xmax  = max(Xd, [], 2);
    xrng  = xmax - xmin;
    xrng(xrng < eps) = 1;
    Xd_n  = 2*(Xd - xmin) ./ xrng - 1;

    Xs = Xd_n(:, 1:end-1);
    Ys = Xd_n(:, 2:end);

    for ci = 1:N_cfgs
        tic;
        [Psi_X_sw, ~] = build_dictionary(Xs, sweep_cfgs{ci});
        [Psi_Y_sw, ~] = build_dictionary(Ys, sweep_cfgs{ci});
        [K_sw, ~, ~]  = edmd(Psi_X_sw, Psi_Y_sw, edmd_opts);
        elapsed = toc * 1e3;

        Psi_Y_pred_sw = K_sw * Psi_X_sw;
        sweep_Npsi(r, ci) = size(Psi_X_sw, 1);
        sweep_err(r,  ci) = norm(Psi_Y_sw - Psi_Y_pred_sw, 'fro') / norm(Psi_Y_sw, 'fro');
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

    for ci = 1:N_cfgs
        col = type_colors(sweep_cfgs{ci}.type);
        np  = sweep_Npsi(r, ci);
        semilogy(ax_err, np, sweep_err(r, ci), [col 'o'], ...
            'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(ax_time, np, sweep_time(r, ci), [col 's'], ...
            'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end

    for fam = {'poly', 'rbf', 'full'}
        idx = find(cellfun(@(c) strcmp(c.type, fam{1}), sweep_cfgs));
        col = type_colors(fam{1});
        [npsi_sorted, so] = sort(sweep_Npsi(r, idx));
        semilogy(ax_err, npsi_sorted, sweep_err(r, idx(so)), [col '-'], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
        plot(ax_time, npsi_sorted, sweep_time(r, idx(so)), [col '--'], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
    end

    ylabel(ax_err,  'Relative error (log scale)');
    xlabel(ax_err,  'N_{\psi}  (observables)');
    title(ax_err,   sprintf('%s: error', regimes{r, 2}));
    grid(ax_err, 'on');

    ylabel(ax_time, 'Computation time (ms)');
    xlabel(ax_time, 'N_{\psi}  (observables)');
    title(ax_time,  sprintf('%s: timing', regimes{r, 2}));
    grid(ax_time, 'on');
end

ax1 = subplot(4, 2, 1, 'Parent', fig_sw);
h = gobjects(5, 1);
fam_names = {'linear', 'pwl', 'poly', 'rbf', 'full'};
for fi = 1:5
    col   = type_colors(fam_names{fi});
    h(fi) = plot(ax1, NaN, NaN, [col 'o-'], 'MarkerSize', 6, 'LineWidth', 1.2);
end
legend(ax1, h, fam_names, 'Location', 'best', 'FontSize', 7);

exportgraphics(fig_sw, 'writeup/dictionarySweep.png', 'Resolution', 150);

% =========================================================================
%% Double Scroll Extended Sweep — richer dictionaries for regime 4 only

ds_idx = 4;   % Double Scroll is regime 4
Xds   = X_datasets{ds_idx};
xmin_ds = min(Xds, [], 2);
xmax_ds = max(Xds, [], 2);
xrng_ds = xmax_ds - xmin_ds;
xrng_ds(xrng_ds < eps) = 1;
Xds_n = 2*(Xds - xmin_ds) ./ xrng_ds - 1;
Xs_ds = Xds_n(:, 1:end-1);
Ys_ds = Xds_n(:, 2:end);

sweep_cfgs_ds = {
    struct('type','rbf',  'n_centers',100, 'sigma',1.0), ...
    struct('type','rbf',  'n_centers',200, 'sigma',1.0), ...
    struct('type','rbf',  'n_centers',20,  'sigma',0.5), ...
    struct('type','rbf',  'n_centers',20,  'sigma',2.0), ...
    struct('type','poly', 'degree',6), ...
    struct('type','poly', 'degree',7), ...
    struct('type','full', 'degree',4, 'n_centers',50,  'sigma',1.0), ...
    struct('type','full', 'degree',3, 'n_centers',50,  'sigma',0.5), ...
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
    [K_ds, ~, ~]  = edmd(Psi_X_ds, Psi_Y_ds, edmd_opts);
    elapsed = toc * 1e3;

    Psi_Y_pred_ds = K_ds * Psi_X_ds;
    ds_Npsi(ci) = size(Psi_X_ds, 1);
    ds_err(ci)  = norm(Psi_Y_ds - Psi_Y_pred_ds, 'fro') / norm(Psi_Y_ds, 'fro');
    ds_time(ci) = elapsed;
    fprintf('  Config %d (%s): N_psi=%d, err=%.4e, t=%.1f ms\n', ...
        ci, sweep_cfgs_ds{ci}.type, ds_Npsi(ci), ds_err(ci), ds_time(ci));
end

% =========================================================================
%% Local helper functions

function idx_sel = pick_conjugate_safe_modes(lambda, ord_modes, target_count, stable_mask)
idx_sel = [];
tol_imag = 1e-8;
tol_pair = 1e-6;

for ii = 1:numel(ord_modes)
    k = ord_modes(ii);

    if ~stable_mask(k)
        continue;
    end
    if ismember(k, idx_sel)
        continue;
    end

    idx_sel(end+1) = k; %#ok<AGROW>

    if abs(imag(lambda(k))) > tol_imag
        diffs = abs(lambda - conj(lambda(k)));
        diffs(~stable_mask) = inf;
        diffs(k) = inf;
        [bestdiff, j] = min(diffs);
        if bestdiff < tol_pair && ~ismember(j, idx_sel)
            idx_sel(end+1) = j; %#ok<AGROW>
        end
    end

    if numel(idx_sel) >= target_count
        break;
    end
end
end

function k_best = pick_visual_eigenfunction(lambda, Xi)
% Pick a nontrivial bounded oscillatory eigenfunction that varies
% substantially along the trajectory, for visualization purposes.
tol_imag = 1e-8;
tol_one  = 1e-3;

osc_mask = abs(imag(lambda)) > tol_imag & abs(lambda) <= 1.05;
osc_mask = osc_mask & abs(lambda - 1) > tol_one;

score = -inf(size(lambda));
for k = 1:numel(lambda)
    if osc_mask(k)
        score(k) = abs(lambda(k)) * std(real(Xi(k, :)));
    end
end

if any(isfinite(score))
    [~, k_best] = max(score);
    return;
end

bounded_mask = abs(lambda) <= 1.05 & abs(lambda - 1) > tol_one;
for k = 1:numel(lambda)
    if bounded_mask(k)
        score(k) = abs(lambda(k)) * std(real(Xi(k, :)));
    end
end

[~, k_best] = max(score);
end

function data = load_chua_csv(filepath)
fid = fopen(filepath, 'r');
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

lines = lines(2:end);

data = cellfun(@(l) sscanf(l, '%f,%f,%f,%f')', lines, 'UniformOutput', false);
data = vertcat(data{:});
end
