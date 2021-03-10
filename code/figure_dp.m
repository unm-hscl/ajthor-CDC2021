%% Figure (3)
% Figure showing the phase portrait of a nonholonomic vehicle.

rng(0);

% ur = linspace(0, 1, 100);
r1 = 0:0.1:1;
r2 = -10:1:10;
[U1, U2] = meshgrid(r1, r2);
ur = [
    reshape(U1, 1, []);
    reshape(U2, 1, []);
    ];

%% Compute kernel solution.
disp('Computing kernel based solution...');

% Kernel parameters.
lambda = 1/(M^2);
sigma = 3;

discount_factor = 0.95;

x0 = [
    -0.8;
    0;
    pi;
    ];

tic

T_end = 25;
r_tracking = linspace(-1, 1, T_end);
X_tracking = [r_tracking; r_tracking; pi/4*ones(size(r_tracking))];

g = zeros(T_end, size(ys, 2));
for k = T_end:-1:1
    g(k, :) = vecnorm(ys([1 2], :) - X_tracking([1 2], k), 2).^2;
    
    % Sanity constraint to ensure non-trivial solution. 
    % Since the approximation tends to zero outside the region for which we
    % have samples, the dynamic program tends to favor leaving the sampled
    % region instead of staying within the desired area. Thus, we add a
    % high 'penalty' for leaving the region of interest. 
%     g(k, :) = g(k, :) + 10*double(any(abs(ys([1, 2], :)) > 1));
end

% dynamics = @(x, u) nh_dynamics(x, u, Ts);
% 
% algorithm = KernelDyn('Lambda', 1/M^2, 'Sigma', 1);
% result = algorithm.compute(xs, us, ys, x0, ur, T_end, g, dynamics);
% 
% X = result.traj;

%% START 

M = size(xs, 2);
n = size(ur, 2); %#ok<NASGU>
N = T_end;

%%

% Compute the kernel matrix.
Gx = rbf_kernel(xs, xs, sigma);
Gu = rbf_kernel(us, us, sigma);

G = Gx.*Gu;

W = inv(G + lambda*M*eye(size(G))); %#ok<*MINV>

Vk = zeros(N, M);

cxy = rbf_kernel(xs, ys, sigma);
cut = rbf_kernel(us, ur, sigma);

% a = W*cxy;
% beta = repmat(W*cxy, 1, 1, n);
% beta = permute(beta, [1 3 2]).*cut;

%%

Vk(N, :) = g(N, :);

for k = N-1:-1:1
    
    fprintf('Computing V(%d)...\n', k);
    
    c = Vk(k+1, :);
    Z = c*W;
    
    w = zeros(M, n);
%     for p = 1:n
%         w(:, p) = c*squeeze(beta(:, p, :));
% %         w(:, p) = c*(a*cut);
%     end

    for p = 1:M
% %         w(p, :) = c*(a(:, p).*cut);
%         sth = (a(:, p).*cut);
% %         sth = sth./sum(abs(sth), 1);
%         w(p, :) = c*sth;
% %         w(p, :) = c*a;
        w(p, :) = Z*(rbf_kernel(xs, ys(:, p), sigma).*cut);
    end
    
    [V, Idx] = min(w, [], 2);
    Uopt = ur(:, Idx);

%     cup = rbf_kernel(us, Uopt, sigma);
% 
%     gamma = cxy.*cup;
%     gamma = W*gamma;
% 
%     Vk(k, :) = g(k, :) + (Vk(k+1, :))*gamma;
% %     Vk(k, :) = c*gamma;

    Vk(k, :) = g(k, :) + V.';
    
%     figure, surf(reshape(Vk(k, 1:144), 12, 12)), view([0 90]), title(num2str(k))

end

%%

Xtraj = x0;
Utraj = [];

% for k = N-1:-1:1
for k = 1:N
    
%     fprintf('Computing X(%d)...\n', k);
    
    cxt = rbf_kernel(xs, Xtraj(:, end), sigma);
%     beta = repmat(W*cxt, 1, 1, n);
%     beta = permute(beta, [1 3 2]).*cut;
    beta = W*(cxt.*cut);
%     beta = beta./sum(abs(beta), 1);

    c = Vk(k, :); 
%     c = g(k, :);

    w = zeros(1, n);
    
%     for p = 1:n
%         w(:, p) = c*squeeze(beta(:, p, :));
%         w(:, p) = c*beta(:, p);
%     end
    w = c*beta;
    [~, Idx] = min(w, [], 2);
    Uopt = ur(:, Idx);
    
    Xtraj = [Xtraj, nh_dynamics(Xtraj(:, end), Uopt, Ts)]; %#ok<AGROW>
    Utraj = [Utraj, Uopt]; %#ok<AGROW>

%     cup = rbf_kernel(us, Uopt, algorithm.sigma_);
% 
%     gamma = cxt.*cup;
%     gamma = W*gamma;

%     Pr(k, :) = (Vk(k+1, :)*gamma);
    
end

toc

X = Xtraj;

%% Plot the results.

% figure
figure('Units', 'points', ...
       'Position', [0, 0, 240, 180]);
ax = axes;
ax.NextPlot = 'add';
ax.Units = 'points';
grid on

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = '$$x_1$$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$$x_2$$';
set(ax, 'FontSize', 8);

xlim([-1, 1])
ylim([-1, 1])

plot(X_tracking(1, :), X_tracking(2, :), 'kx:');
plot(X(1, :), X(2, :), '^--');

%%
% pause(2); close
% for k = 1:N, figure, surf(reshape(Vk(k, 1:144), 12, 12)), view([0 90]), title(num2str(k)), end

function G = rbf_kernel(X, Y, sigma)
% RBF_KERNEL Compute kernel matrix.

    M = size(X, 2);
    T = size(Y, 2);

    G = zeros(M, T);

    for k = 1:size(X, 1)
        G = G + (repmat(Y(k, :), [M, 1]) - repmat(X(k, :)', [1, T])).^2;
    end

    G = exp(-G/(2*sigma^2));

end