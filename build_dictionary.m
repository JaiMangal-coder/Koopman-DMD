function [Psi, labels] = build_dictionary(X, opts)
% BUILD_DICTIONARY  Lift state snapshots using a dictionary of observables.
%
%   [Psi, labels] = build_dictionary(X, opts)
%
%   Inputs:
%     X     (n x M)   state snapshot matrix
%     opts  struct    dictionary options:
%             .type       'linear' | 'poly' | 'trig' | 'hybrid' |
%                         'rbf' | 'pwl' | 'full'
%             .degree     max polynomial degree     (poly / hybrid / full; default 3)
%             .n_freq     max trig frequency        (trig / hybrid; default 3)
%             .n_centers  number of RBF centers     (rbf  / full; default 20)
%             .sigma      RBF width                 (rbf  / full; default 1.0)
%
%   Outputs:
%     Psi    (N_psi x M)   lifted snapshot matrix
%     labels (N_psi x 1)   cell array of observable names
%
%   Dictionary types
%   ----------------
%   'linear'  — state observables [x1; …; xn]  (no constant; DMD baseline)
%   'poly'    — monomials up to total degree opts.degree  (includes constant)
%   'trig'    — [1; sin(k*xj); cos(k*xj)] for j=1..n, k=1..opts.n_freq
%   'hybrid'  — poly ∪ trig(no constant)
%   'vdp'     — physics-tailored for Van der Pol (requires n == 2):
%               {1, x1, x2, x1^2, x1*x2, x2^2, x1^2*x2, r^2, r^2*x1,
%                sin(theta), cos(theta)}  where r=norm, theta=atan2(x2,x1)
%   'rbf'     — [1; Gaussian RBFs] with centers from k-means on X
%   'pwl'     — [x1;x2;x3; |x1+1|; |x1-1|; x1|x1+1|; x1|x1-1|]
%               (Chua-matched; requires n == 3)
%   'full'    — poly ∪ rbf(no constant) ∪ pwl kink terms only
%               (Chua-matched; requires n == 3)
%
%   Note: 'rbf' and 'full' require the Statistics & Machine Learning Toolbox
%   for kmeans().

[n, M] = size(X);
opts   = parse_opts(opts);

switch lower(opts.type)
    case 'linear'
        [Psi, labels] = linear_dict(X, n, M);

    case 'poly'
        [Psi, labels] = poly_dict(X, n, M, opts.degree);

    case 'trig'
        [Psi, labels] = trig_dict(X, n, M, opts.n_freq);

    case 'vdp'
        if n ~= 2
            error('vdp dictionary requires n == 2.');
        end
        [Psi, labels] = vdp_dict(X, M);

    case 'hybrid'
        [Psi_p, labels_p] = poly_dict(X, n, M, opts.degree);
        [Psi_t, labels_t] = trig_dict(X, n, M, opts.n_freq);
        % poly has constant; drop constant row from trig
        Psi    = [Psi_p;          Psi_t(2:end, :)];
        labels = [labels_p(:); labels_t(2:end)];

    case 'rbf'
        [Psi, labels] = rbf_dict(X, n, M, opts.n_centers, opts.sigma);

    case 'pwl'
        if n ~= 3
            error('pwl dictionary requires n == 3 (Chua state dimension).');
        end
        [Psi, labels] = pwl_dict(X, M);

    case 'full'
        if n ~= 3
            error('full dictionary requires n == 3 (Chua state dimension).');
        end
        [Psi_p, labels_p] = poly_dict(X, n, M, opts.degree);
        Psi    = Psi_p;
        labels = labels_p(:);
        if opts.use_rbf
            [Psi_r, labels_r] = rbf_dict(X, n, M, opts.n_centers, opts.sigma);
            % poly already has constant; drop constant row from rbf
            Psi    = [Psi;    Psi_r(2:end, :)];
            labels = [labels; labels_r(2:end)];
        end
        if opts.use_pwl
            [Psi_k, labels_k] = kink_dict(X, M);
            Psi    = [Psi;    Psi_k];
            labels = [labels; labels_k(:)];
        end

    otherwise
        error('Unknown dictionary type "%s". Use ''linear'', ''poly'', ''trig'', ''hybrid'', ''vdp'', ''rbf'', ''pwl'', or ''full''.', opts.type);
end
end

% =========================================================================
% Dictionary builders

function [Psi, labels] = linear_dict(X, n, M)
% Observables: x1, x2, ..., xn  (no constant term)
Psi    = X;   % already (n x M)
labels = arrayfun(@(j) sprintf('x%d', j), 1:n, 'UniformOutput', false)';
end

function [Psi, labels] = poly_dict(X, n, M, degree)
% Monomials x1^a1 * ... * xn^an with sum(a) <= degree.
exponents = monomial_exponents(n, degree);
N_psi     = size(exponents, 1);

Psi = ones(N_psi, M);
for k = 1:N_psi
    for j = 1:n
        if exponents(k, j) > 0
            Psi(k, :) = Psi(k, :) .* X(j, :) .^ exponents(k, j);
        end
    end
end

labels = cell(N_psi, 1);
for k = 1:N_psi
    labels{k} = monomial_label(exponents(k, :));
end
end

function [Psi, labels] = trig_dict(X, n, M, n_freq)
% [1; x1;..;xn; sin(k*xj); cos(k*xj)] for j=1..n, k=1..n_freq.
% State variables are included so x1/x2 labels exist for prediction extraction.
N_psi  = 1 + n + 2*n*n_freq;
Psi    = zeros(N_psi, M);
labels = cell(N_psi, 1);

Psi(1, :) = 1;
labels{1} = '1';
for j = 1:n
    Psi(1+j, :)  = X(j, :);
    labels{1+j}  = sprintf('x%d', j);
end
idx = 1 + n + 1;
for k = 1:n_freq
    for j = 1:n
        if k == 1
            Psi(idx, :)   = sin(X(j, :));
            labels{idx}   = sprintf('sin(x%d)', j);
            Psi(idx+1, :) = cos(X(j, :));
            labels{idx+1} = sprintf('cos(x%d)', j);
        else
            Psi(idx, :)   = sin(k * X(j, :));
            labels{idx}   = sprintf('sin(%dx%d)', k, j);
            Psi(idx+1, :) = cos(k * X(j, :));
            labels{idx+1} = sprintf('cos(%dx%d)', k, j);
        end
        idx = idx + 2;
    end
end
end

function [Psi, labels] = vdp_dict(X, M)
% Physics-tailored dictionary for the Van der Pol oscillator (n=2).
%
%   {1, x1, x2, x1^2, x1*x2, x2^2, x1^2*x2, r^2, r^2*x1, sin(theta), cos(theta)}
%
%   x1^2*x2  — the VdP nonlinearity (appears explicitly in dx2/dt)
%   r^2      — amplitude squared; converges to a constant on the limit cycle
%   r^2*x1   — amplitude-weighted state for radial modulation
%   sin/cos(theta) — phase around the limit cycle via atan2(x2,x1),
%                    not sinusoids of state values
x1 = X(1, :);
x2 = X(2, :);
r2 = x1.^2 + x2.^2;
r  = sqrt(r2);
r_safe = max(r, eps);   % avoid /0 at the unstable fixed point

Psi = [ones(1, M);
       x1;
       x2;
       x1.^2;
       x1 .* x2;
       x2.^2;
       x1.^2 .* x2;
       r2;
       r2 .* x1;
       x2 ./ r_safe;    % sin(theta)
       x1 ./ r_safe];   % cos(theta)

labels = {'1'; 'x1'; 'x2'; 'x1^2'; 'x1*x2'; 'x2^2'; ...
          'x1^2*x2'; 'r2'; 'r2*x1'; 'sin_theta'; 'cos_theta'};
end

function [Psi, labels] = rbf_dict(X, n, M, n_centers, sigma)
% [1; exp(-||x - c_l||^2 / sigma^2)] for l = 1..n_centers.
% Centers chosen by k-means on the training columns of X.
[~, C] = kmeans(X', n_centers);   % C is (n_centers x n)
C = C';                            % (n x n_centers)

N_psi  = 1 + n_centers;
Psi    = zeros(N_psi, M);
labels = cell(N_psi, 1);

Psi(1, :)  = 1;
labels{1}  = '1';

for l = 1:n_centers
    diff_sq    = sum((X - C(:, l)) .^ 2, 1);   % (1 x M)
    Psi(l+1, :) = exp(-diff_sq / sigma^2);
    labels{l+1} = sprintf('rbf_%d', l);
end
end

function [Psi, labels] = pwl_dict(X, M)
% Piecewise-linear dictionary for the Chua circuit (n=3, state = [x1;x2;x3]).
% Observables: x1, x2, x3, |x1+1|, |x1-1|, x1*|x1+1|, x1*|x1-1|
x1 = X(1, :);
x2 = X(2, :);
x3 = X(3, :);

Psi = [x1;
       x2;
       x3;
       abs(x1 + 1);
       abs(x1 - 1);
       x1 .* abs(x1 + 1);
       x1 .* abs(x1 - 1)];

labels = {'x1'; 'x2'; 'x3'; '|x1+1|'; '|x1-1|'; 'x1|x1+1|'; 'x1|x1-1|'};
end

function [Psi, labels] = kink_dict(X, M)
% The four Chua kink terms only (no linear base), for use in 'full'.
x1 = X(1, :);

Psi = [abs(x1 + 1);
       abs(x1 - 1);
       x1 .* abs(x1 + 1);
       x1 .* abs(x1 - 1)];

labels = {'|x1+1|'; '|x1-1|'; 'x1|x1+1|'; 'x1|x1-1|'};
end

% =========================================================================
% Helpers

function opts = parse_opts(opts)
if ~isstruct(opts)
    % Scalar shorthand: treat as polynomial degree
    opts = struct('type', 'poly', 'degree', opts);
end
if ~isfield(opts, 'type'),      opts.type      = 'poly'; end
if ~isfield(opts, 'degree'),    opts.degree    = 3;      end
if ~isfield(opts, 'n_freq'),    opts.n_freq    = 3;      end
if ~isfield(opts, 'n_centers'), opts.n_centers = 20;     end
if ~isfield(opts, 'sigma'),     opts.sigma     = 1.0;    end
if ~isfield(opts, 'use_rbf'),   opts.use_rbf   = true;   end
if ~isfield(opts, 'use_pwl'),   opts.use_pwl   = true;   end
end

function exps = monomial_exponents(n, max_deg)
exps = [];
for d = 0:max_deg
    exps = [exps; fixed_sum_exponents(n, d)]; %#ok<AGROW>
end
end

function result = fixed_sum_exponents(n, total)
if n == 1
    result = total;
    return
end
result = [];
for k = 0:total
    sub    = fixed_sum_exponents(n - 1, total - k);
    result = [result; k * ones(size(sub, 1), 1), sub]; %#ok<AGROW>
end
end

function s = monomial_label(exp_vec)
if all(exp_vec == 0)
    s = '1';
    return
end
parts = {};
for j = 1:length(exp_vec)
    e = exp_vec(j);
    if e == 1
        parts{end+1} = sprintf('x%d', j); %#ok<AGROW>
    elseif e > 1
        parts{end+1} = sprintf('x%d^%d', j, e); %#ok<AGROW>
    end
end
s = strjoin(parts, '*');
end
