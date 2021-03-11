function results = compute(obj, xs, us, ys, x0, ur, N, cost, dynamics)
% COMPUTE Computes the results.

% Algorithm implementation.

M = size(xs, 2);
n = size(ur, 2); 

t_start = tic;

% Compute the kernel matrix.
Gx = rbf_kernel(xs, xs, obj.Sigma);
Gu = rbf_kernel(us, us, obj.Sigma);
G = Gx.*Gu;
clear Gx Gu

W = inv(G + obj.Lambda*M*eye(size(G)));

Vk = zeros(N, M);

% cxy = rbf_kernel(xs, ys, obj.Sigma);
cut = rbf_kernel(us, ur, obj.Sigma);

Vk(N, :) = cost(N);

for t = N-1:-1:1
    
    fprintf('Computing V(%d)...\n', t);
    
    c = Vk(t+1, :);
    Z = c*W; %#ok<MINV>
    
    w = zeros(M, n);

    for p = 1:M
        w(p, :) = Z*(rbf_kernel(xs, ys(:, p), obj.Sigma).*cut);
    end
    
    [V, ~] = min(w, [], 2);
    Vk(t, :) = cost(t) + V.';
    
end

X = [];
U = [];

X = [X, x0]; 

for t = 1:N
    
    cxt = rbf_kernel(xs, X(:, end), obj.Sigma);
    beta = W*(cxt.*cut); %#ok<MINV>

    c = Vk(t, :); 
    
    w = c*beta;
    [~, Idx] = min(w, [], 2);
    Uopt = ur(:, Idx);
    
    X = [X, dynamics(X(:, end), Uopt)]; %#ok<AGROW>
    U = [U, Uopt]; %#ok<AGROW>
    
end

t_elapsed = toc(t_start);

results.x_traj = X;
results.u_traj = U;
results.comp_time = t_elapsed;

end
