function results = compute(obj, xs, us, xt, ur, c)
% COMPUTE Computes the results.

% Algorithm implementation.

M = size(xs, 2);
n = size(ur, 2); %#ok<NASGU>

t_start = tic;

% Compute the kernel matrix.
Gx = rbf_kernel(xs, xs, obj.Sigma);
Gu = rbf_kernel(us, us, obj.Sigma);
G = Gx.*Gu;
clear Gx Gu

W = inv(G + obj.Lambda*M*eye(size(G)));

% Compute the "kernel" matrix of input samples vs. input points.
Ups = rbf_kernel(us, ur, obj.Sigma);

% Variable to hold results corresponding to each test point.
U = zeros(size(ur, 1), size(xt, 2));

Z = c*W; %#ok<MINV>
w = zeros(size(xt, 2), size(ur, 2));
beta = repmat(rbf_kernel(xs, xt, obj.Sigma), 1, 1, size(ur, 2));
beta = permute(beta, [1 3 2]).*Ups;
for k = 1:size(ur, 2)
    w(:, k) = Z*squeeze(beta(:, k, :));
end

[~, Idx] = min(w, [], 2);
U = ur(:, Idx);

t_elapsed = toc(t_start);

results.u_opt = U;
results.comp_time = t_elapsed;

end
