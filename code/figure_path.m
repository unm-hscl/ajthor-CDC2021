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
sigma = 1;

tic

% Compute the cost function.
% c = abs(ys(3, :));
% c = vecnorm(ys([1, 2], :) - [0; 0]);
c = vecnorm(ys);

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
    -0.25;
    -0.25;
    0;
    ];

U = [];

T_end = 500;
for t = 1:T_end
    
    x0 = X(:, end);

    w = zeros(size(x0, 2), size(ur, 2));
    beta = repmat(rbf_kernel(xs, x0, sigma), 1, 1, size(ur, 2));
    beta = permute(beta, [1 3 2]).*Ups;
    for k = 1:size(ur, 2)
        w(:, k) = Z*squeeze(beta(:, k, :));
    end

    [~, Idx] = min(w, [], 2);
    U = [U, ur(:, Idx)];
    
    X = [X, nh_dynamics(x0, U(:, end), Ts)];
    
end

toc

%% Plot the results.

% figure 
% ax = axes;
% ax.NextPlot = 'add';
% ax.Units = 'points';

plot(X(1, :), X(2, :));
% plot3(X(1, :), X(2, :), 0:1:T_end)
% xlabel('x1')
% ylabel('x3')
% zlabel('t')

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