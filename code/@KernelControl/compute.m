function results = compute(obj, xs, us, xt, ur, c)
% COMPUTE Computes the results.

% Algorithm implementation.

M = size(xs, 2);
n = size(ur, 2); %#ok<NASGU>

t_start = tic;

% Compute the kernel matrix.
Gx = rbf_kernel(xs, xs, obj.sigma_);
Gu = rbf_kernel(us, us, obj.sigma_);
G = Gx.*Gu;
clear Gx Gu

K = (G + obj.lambda_*M*eye(size(G)));
W = inv(K);

% Compute the "kernel" matrix of input samples vs. input points.
Ups = rbf_kernel(us, ur, obj.sigma_);

% Variable to hold results corresponding to each test point.
U = zeros(size(ur, 1), size(xt, 2));

% Version 1

Z = c*W; %#ok<MINV>
w = zeros(size(xt, 2), size(ur, 2));
beta = repmat(rbf_kernel(xs, xt, obj.sigma_), 1, 1, size(ur, 2));
beta = permute(beta, [1 3 2]).*Ups;
for k = 1:size(ur, 2)
    w(:, k) = Z*squeeze(beta(:, k, :));
end

[~, Idx] = min(w, [], 2);
U = ur(Idx);

% Version 2

% Z = W*c.'; %#ok<MINV>
% w = zeros(size(xt, 2), size(ur, 2));
% beta = repmat(rbf_kernel(xs, xt, obj.sigma_), 1, 1, size(ur, 2));
% beta = permute(beta, [1 3 2]).*Ups;
% for k = 1:size(ur, 2)
%     w(:, k) = Z*squeeze(beta(:, k, :));
% end
% 
% [~, Idx] = min(w, [], 2);
% U = ur(Idx);

% Version 3

% for k = 1:size(xt, 2)
% 
%     % Compute the "kernel" matrix of samples vs. the test point.
%     Phi = rbf_kernel(xs, xt(:, k), obj.sigma_);
%     beta = W*(Phi.*Ups);        %#ok<MINV>
% 
%     % Evaluate the cost function at the test point under each control
%     % input so that we get an approximation of the cost c under each
%     % input.
%     w = c*beta;
% 
% %     % Optimization.
% %     cvx_begin quiet
% %         variable a(n);
% %         minimize( dot(w, a) );
% %         subject to
% %             sum(a) == 1;        %#ok<*EQEFF>
% %             0 <= a;             %#ok<*VUNUS>
% %     cvx_end
% % 
% %     % Store the results.
% %     [~, Idx] = max(a);
% %     U(:, k) = ur(Idx);
% 
%     [~, Idx] = min(w);
%     U(:, k) = ur(Idx);
% 
% end

t_elapsed = toc(t_start);

results.u_opt = U;
results.comp_time = t_elapsed;

end
