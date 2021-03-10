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

tic

% Compute the cost function.
% c = abs(ys(3, :));
% c = vecnorm(ys([1, 2], :) - [0; 0]);
% c = vecnorm(ys);

% Compute the kernel matrix.
Gx = rbf_kernel(xs, xs, sigma);
Gu = rbf_kernel(us, us, sigma);
G = Gx.*Gu;
clear Gx Gu

K = (G + lambda*M*eye(size(G)));
W = inv(K);

% Compute the "kernel" matrix of input samples vs. input points.
Ups = rbf_kernel(us, ur, sigma);

X = [
    -0.8;
    0;
    pi;
    ];

U = [];

T_end = 25;
r_tracking = linspace(-1, 1, T_end);
X_tracking = [r_tracking; r_tracking; pi/4*ones(size(r_tracking))];

for t = 1:T_end
    
    x0 = X(:, end);
    
    c = vecnorm(ys([1 2], :) - X_tracking([1 2], t), 2).^2;

    w = zeros(size(x0, 2), size(ur, 2));
    beta = repmat(rbf_kernel(xs, x0, sigma), 1, 1, size(ur, 2));
    beta = permute(beta, [1 3 2]).*Ups;
    beta = W*beta; %#ok<MINV>
    for k = 1:size(ur, 2)
        w(:, k) = c*squeeze(beta(:, k, :));
    end

    [~, Idx] = min(w, [], 2);
    U = [U, ur(:, Idx)]; %#ok<AGROW>
    
    X = [X, nh_dynamics(x0, U(:, end), Ts)]; %#ok<AGROW>
    
end

toc

%% Plot the results.

figure('Units', 'points', ...
       'Position', [0, 0, 240, 120]);
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