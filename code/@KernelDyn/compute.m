function results = compute(obj, xs, us, ys, x0, ur, N, g, f)
% COMPUTE Computes the results.

% Algorithm implementation.

M = size(xs, 2);
n = size(ur, 2); %#ok<NASGU>

t_start = tic;

% Compute the kernel matrix.
Gx = rbf_kernel(xs, xs, obj.sigma_);
Gu = rbf_kernel(us, us, obj.sigma_);

G = Gx.*Gu;

W = inv(G + obj.lambda_*M*eye(size(G))); %#ok<*MINV>

Vk = zeros(N, M);

cxy = rbf_kernel(xs, ys, obj.sigma_);
cut = rbf_kernel(us, ur, obj.sigma_);

beta = repmat(W*cxy, 1, 1, n);
beta = permute(beta, [1 3 2]).*cut;

Vk(N, :) = g(N, :);

for k = N-1:-1:2
    
    fprintf('Computing V(%d)...\n', k);
    
    c = Vk(k+1, :); 
    w = zeros(M, n);
    for p = 1:n
        w(:, p) = c*squeeze(beta(:, p, :));
    end
    [~, Idx] = max(w, [], 2);
    Uopt = ur(:, Idx);

    cup = rbf_kernel(us, Uopt, obj.sigma_);

    gamma = cxy.*cup;
    gamma = W*gamma;

    Vk(k, :) = g(k, :) + (Vk(k+1, :))*gamma;
end

Xtraj = x0;

for k = N-1:-1:1
    
    fprintf('Computing X(%d)...\n', k);
    
    cxt = rbf_kernel(xs, Xtraj(:, end), obj.sigma_);
    beta = repmat(W*cxt, 1, 1, n);
    beta = permute(beta, [1 3 2]).*cut;

    c = Vk(k+1, :); 
    w = zeros(1, n);
    for p = 1:n
        w(:, p) = c*squeeze(beta(:, p, :));
    end
    [~, Idx] = max(w, [], 2);
    Uopt = ur(:, Idx);
    
    Xtraj = [Xtraj, f(Xtraj(:, end), Uopt)];

%     cup = rbf_kernel(us, Uopt, algorithm.sigma_);
% 
%     gamma = cxt.*cup;
%     gamma = W*gamma;

%     Pr(k, :) = (Vk(k+1, :)*gamma);
    
end

% % Compute the "kernel" matrix of input samples vs. input points.
% Ups = rbf_kernel(us, ur, obj.sigma_);
% 
% % Variable to hold results corresponding to each test point.
% U = zeros(size(ur, 1), size(xt, 2));
% 
% Z = c*W; %#ok<MINV>
% w = zeros(size(xt, 2), size(ur, 2));
% beta = repmat(rbf_kernel(xs, xt, obj.sigma_), 1, 1, size(ur, 2));
% beta = permute(beta, [1 3 2]).*Ups;
% for k = 1:size(ur, 2)
%     w(:, k) = Z*squeeze(beta(:, k, :));
% end
% 
% [~, Idx] = min(w, [], 2);
% U = ur(:, Idx);

t_elapsed = toc(t_start);

results.traj = Xtraj;
% results.u_opt = U;
results.comp_time = t_elapsed;

end
