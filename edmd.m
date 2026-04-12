function [K, lambda, Phi] = edmd(Psi_X, Psi_Y, opts)
% EDMD  Least-squares approximation of the Koopman operator
%
%   Inputs:
%     Psi_X  (N_psi x M)  lifted snapshots at current time
%     Psi_Y  (N_psi x M)  lifted snapshots one step forward
%     opts   (struct)      optional fields:
%                            .method  'tikhonov' (default) | 'tsvd'
%                            .alpha   Tikhonov parameter (default 1e-10)
%                            .tol     SVD truncation tolerance (default 1e-10)
%
%   Outputs:
%     K      (N_psi x N_psi)  EDMD Koopman matrix
%     lambda (N_psi x 1)      Koopman eigenvalues
%     Phi    (N_psi x N_psi)  Koopman eigenvectors (columns)

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'method'), opts.method = 'tikhonov'; end
if ~isfield(opts, 'alpha'),  opts.alpha  = 1e-10;      end
if ~isfield(opts, 'tol'),    opts.tol    = 1e-10;      end

M = size(Psi_X, 2);
G = (Psi_X * Psi_X') / M;
A = (Psi_Y * Psi_X') / M;

switch opts.method

    case 'tikhonov'
        % K = A * (G + alpha*I)^{-1}
        N = size(G, 1);
        K = A / (G + opts.alpha * eye(N));

    case 'tsvd'
        % Truncate small singular values of G before inverting
        [U, S, V] = svd(G);
        sv        = diag(S);
        r         = sum(sv > opts.tol * sv(1));   % rank estimate
        K         = A * V(:,1:r) * diag(1./sv(1:r)) * U(:,1:r)';

end

% Diagnostics
fprintf('G condition number: %.2e  |  rank estimate: %d / %d\n', ...
        cond(G), rank(G, opts.tol * max(diag(G))), size(G,1));

[Phi, D] = eig(K);
lambda   = diag(D);
end