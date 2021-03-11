function results = compute(obj, xs, us, x0, ur, N, cost, dynamics)
% COMPUTE Computes the results.

% Algorithm implementation.

M = size(xs, 2);

t_start = tic;

% Compute the kernel matrix.
Gx = rbf_kernel(xs, xs, obj.Sigma);
Gu = rbf_kernel(us, us, obj.Sigma);
G = Gx.*Gu;
clear Gx Gu

W = inv(G + obj.Lambda*M*eye(size(G)));

% Compute the "kernel" matrix of input samples vs. input points.
Ups = rbf_kernel(us, ur, obj.Sigma);

X = [];
U = [];

X = [X, x0]; 

for t = 1:N
    
    c = cost(t);
    
    w = zeros(size(X(:, end), 2), size(ur, 2));
    beta = repmat(rbf_kernel(xs, X(:, end), obj.Sigma), 1, 1, size(ur, 2));
    beta = permute(beta, [1 3 2]).*Ups;
    beta = W*beta; %#ok<MINV>
    for p = 1:size(ur, 2)
        w(:, p) = c*squeeze(beta(:, p, :));
    end

    [~, Idx] = min(w, [], 2);
    U = [U, ur(:, Idx)]; %#ok<AGROW>
    X = [X, dynamics(X(:, end), U(:, end))]; %#ok<AGROW>
    
end

t_elapsed = toc(t_start);

results.x_traj = X;
results.u_traj = U;
results.comp_time = t_elapsed;

end
