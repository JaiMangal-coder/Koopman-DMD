function plot_eigenvalues(varargin)
% PLOT_EIGENVALUES  Plot Koopman eigenvalues in the complex plane with unit circle.
%
%   plot_eigenvalues(lambda)
%   plot_eigenvalues(lambda, title_str)
%   plot_eigenvalues(ax, lambda, title_str)
%
%   Inputs:
%     ax         axes handle to draw into (optional; creates a new figure if omitted)
%     lambda     (N_psi x 1) complex eigenvalues from edmd()
%     title_str  string title (optional)
%
%   Eigenvalues inside the unit circle are plotted as filled blue circles;
%   eigenvalues outside as open red circles.

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

% --- eigenvalues: split inside / outside unit circle ---
inside = abs(lambda) <= 1;
if any(inside)
    scatter(ax, real(lambda(inside)),  imag(lambda(inside)),  30, ...
        'b', 'filled', 'DisplayName', 'Inside unit circle');
end
if any(~inside)
    scatter(ax, real(lambda(~inside)), imag(lambda(~inside)), 30, ...
        'r', 'LineWidth', 1.2, 'DisplayName', 'Outside unit circle');
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
