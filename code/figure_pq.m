%% Figure (3)
% Figure showing the phase portrait of a nonholonomic vehicle.

rng(0);

% Number of test points.
Mt = 11;

% Range of the test points.
r1 = linspace(-1, 1, Mt);
r2 = linspace(0, 1, Mt);
[XX, YY] = meshgrid(r1, r2);

xt = [
    reshape(XX, 1, []);
    zeros(1, Mt^2);
    reshape(YY, 1, []);
    zeros(1, Mt^2);
    zeros(1, Mt^2);
    zeros(1, Mt^2);
    ];

ru = linspace(0, 10, 10);
[U1, U2] = meshgrid(ru, ru);

ur = [
    reshape(U1, 1, []);
    reshape(U2, 1, []);
    ];

%% Compute kernel solution.
disp('Computing kernel based solution...');

% Kernel parameters.
lambda = 1/M^2;
sigma = 0.1;

tic 

% Compute the cost function for our samples. 
c = vecnorm(ys);

% Compute the results of the kernel based algorithm.
alg = KernelControl('Lambda', lambda, 'Sigma', sigma);
results = alg.compute(xs, us, xt, ur, c);
U = results.u_opt;

toc

yt = zeros(size(xt));

for k = 1:size(xt, 2)
    yt(:, k) = pq_dynamics(xt(:, k), U(:, k), Ts);
end

%% Plot the results.

figure
% figure('Units', 'points', ...
%        'Position', [0, 0, 180, 180]);
ax = axes; 
ax.NextPlot = 'add';
ax.Units = 'points';

% ax.Position = [30, 30, 200, 200];

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = '$$x_1$$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$$x_2$$';
set(ax, 'FontSize', 8);

xlim([-1, 1])
ylim([0, 1])

D = yt - xt;
quiver(xt(1, :), xt(3, :), D(1, :), D(3, :));

%% Compute kernel solution.
disp('Computing kernel based solution...');

% Kernel parameters.
lambda = 1/(M^2);
sigma = 0.25;

tic

% Compute the cost function.

c = vecnorm(ys);
% c = vecnorm(ys - [0; 0; 1; 0; 0; 0]) + (1 - exp(-abs(ys(5)).^2/2));
% c = vecnorm(ys([1, 3], :) - [0; 1]);
% c = vecnorm(ys - [1; 0; 1; 0; 0; 0]);

% Compute the kernel matrix.
Gx = rbf_kernel(xs, xs, sigma);
Gu = rbf_kernel(us, us, sigma);
G = Gx.*Gu;
clear Gx Gu

K = (G + lambda*M*eye(size(G)));
W = inv(K);

Z = c*W; %#ok<MINV>

% Compute the "kernel" matrix of input samples vs. input points.
Ups = rbf_kernel(us, ur, sigma);

X = [
    -0.5;
    0;
    0.5;
    0;
    0;
    0;
    ];

U = [];

for t = 1:50
    
    x0 = X(:, end);

    w = zeros(size(x0, 2), size(ur, 2));
    beta = repmat(rbf_kernel(xs, x0, sigma), 1, 1, size(ur, 2));
    beta = permute(beta, [1 3 2]).*Ups;
    for k = 1:size(ur, 2)
        w(:, k) = Z*squeeze(beta(:, k, :));
    end

    [~, Idx] = min(w, [], 2);
    U = [U, ur(:, Idx)];
    
    X = [X, pq_dynamics(x0, U(:, end), Ts)];
    
end

toc

% %% Plot the results.

figure 
ax = axes;
ax.NextPlot = 'add';
ax.Units = 'points';

plot(X(1, :), X(3, :))

% figure('Units', 'points', ...
%        'Position', [0, 0, 500, 500]);
% ax = axes; 
% ax.NextPlot = 'add';
% ax.Units = 'points';
% 
% % ax.Position = [30, 30, 200, 200];
% 
% ax.XLabel.Interpreter = 'latex';
% ax.XLabel.String = '$$x_1$$';
% ax.YLabel.Interpreter = 'latex';
% ax.YLabel.String = '$$x_2$$';
% set(ax, 'FontSize', 8);
% 
% xlim([-0.75, 0.75])
% ylim([-0.75, 0.75])
% 
% % scatter(xt(1, :), xt(2, :))
% % scatter(yt(1, :), yt(2, :))
% % for k = 1:size(xt, 2)
% %     plot([xt(1, k), yt(1, k)], [xt(2, k), yt(2, k)], 'k');
% % end
% 
% D = yt - xt;
% quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));

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