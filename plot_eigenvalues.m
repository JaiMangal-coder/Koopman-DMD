function plot_eigenvalues(varargin)
% PLOT_EIGENVALUES  Plot Koopman eigenvalues in the complex plane.
%
%   plot_eigenvalues(lambda)
%   plot_eigenvalues(lambda, title_str)
%   plot_eigenvalues(ax, lambda, title_str)
%
%   Conjugate pairs are identified and each pair is drawn in a distinct
%   color.  Real eigenvalues are shown in black.

% --- parse arguments ---
if isa(varargin{1}, 'matlab.graphics.axis.Axes')
    ax        = varargin{1};
    lambda    = varargin{2};
    title_str = '';
    if nargin >= 3, title_str = varargin{3}; end
else
    lambda    = varargin{1};
    title_str = '';
    if nargin >= 2, title_str = varargin{2}; end
    figure;
    ax = axes;
end

% --- unit circle ---
theta = linspace(0, 2*pi, 361);
hold(ax, 'on');
plot(ax, cos(theta), sin(theta), 'k-', 'LineWidth', 1.0);

% --- identify conjugate pairs and real eigenvalues ---
tol_imag = 1e-6;   % threshold to call an eigenvalue real
tol_pair = 1e-4;   % tolerance for matching conjugates

used = false(size(lambda));
pairs  = {};   % each cell: indices of the two conjugate members
reals  = [];   % indices of real eigenvalues

for k = 1:numel(lambda)
    if used(k), continue; end
    if abs(imag(lambda(k))) < tol_imag
        reals(end+1) = k; %#ok<AGROW>
        used(k) = true;
    else
        % find conjugate partner
        diffs = abs(lambda - conj(lambda(k)));
        diffs(used) = inf;
        diffs(k)    = inf;
        [best, j] = min(diffs);
        if best < tol_pair
            pairs{end+1} = [k, j]; %#ok<AGROW>
            used(k) = true;
            used(j) = true;
        else
            % unpaired complex eigenvalue — treat as its own pair
            pairs{end+1} = [k]; %#ok<AGROW>
            used(k) = true;
        end
    end
end

% --- colors: one per pair, cycling through a qualitative palette ---
palette = lines(max(numel(pairs), 1));

for p = 1:numel(pairs)
    idx = pairs{p};
    col = palette(p, :);
    scatter(ax, real(lambda(idx)), imag(lambda(idx)), 40, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', col * 0.6, ...
        'MarkerFaceAlpha', 0.85);
end

% --- real eigenvalues in black ---
if ~isempty(reals)
    scatter(ax, real(lambda(reals)), imag(lambda(reals)), 40, ...
        'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.85);
end

% --- formatting ---
axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'Re(\lambda)');
ylabel(ax, 'Im(\lambda)');
if ~isempty(title_str)
    title(ax, title_str);
end
hold(ax, 'off');
end
