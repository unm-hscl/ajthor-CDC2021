%% Figure (3)
% Figure showing the phase portrait of a nonholonomic vehicle.

rng(0);

% Number of test points.
Mt = 21;

% Range of the test points.
lb = -0.5;
ub = 0.5;
r1 = linspace(lb, ub, Mt);
r2 = linspace(lb, ub, Mt);
[XX, YY] = meshgrid(r1, r2);

xt = [
    reshape(XX, 1, []);
    zeros(1, Mt^2);
    reshape(YY, 1, []);
    zeros(1, Mt^2);
    ];

ru = linspace(-10, 10, 101);

ur = ru;

%% Compute kernel solution.
disp('Computing kernel based solution...');

% Kernel parameters.
lambda = 1/M^2;
sigma = 1;

tic 

% Compute the cost function for our samples. 
c = vecnorm(ys);
% c = abs(ys(3, :));
% c = vecnorm(ys([1 3], :)).^2;

% Compute the results of the kernel based algorithm.
alg = KernelControl('Lambda', lambda, 'Sigma', sigma);
results = alg.compute(xs, us, xt, ur, c);
U = results.u_opt;

toc

% U = zeros(1, size(xt, 2));

% yt = zeros(size(xt));

for k = 1:size(xt, 2)
    yt(:, k) = cp_dynamics(xt(:, k), U(:, k), Ts);
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

% xlim([-1, 1])
% ylim([-1, 1])

D = yt - xt;
% quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));
% quiver(xt(3, :), xt(4, :), D(3, :), D(4, :));
quiver(xt(1, :), xt(3, :), D(1, :), D(3, :));

% X = [-0.25; 0.25; 0; 0];
% for k = 1:20
%     X = [X, cp_dynamics(X(:, end), -1, Ts)];
% end
% plot(X(1, :), X(2, :))
