%% Figure (1)
% Figure showing the vector field of the optimal solution vs. the kernel
% based solution.

rng(0);

% Load the sample.
load('../data/sample_di.mat');

% Create a grid of test points.
r = linspace(-1, 1, 11);
[XX, YY] = meshgrid(r, r);

xt = [
    reshape(XX, 1, []);
    reshape(YY, 1, []);
    ];

% Cleanup.
clear XX YY

%% Compute optimal solution.
disp('Computing optimal solution...');

A = [1, Ts; 0 1];
B = [(Ts^2)/2; Ts];
dynamics = @(x, u) A*x + B*u;

tic

% Compute the results using the optimization based algorithm.
alg_opt = OptimalControl();
results_opt = alg_opt.compute(dynamics, xt, [-1, 1]);
U_opt = results_opt.u_opt;

toc

yt_opt = zeros(size(xt));
for p = 1:size(xt, 2)
    yt_opt(:, p) = dynamics_di(xt(:, p), U_opt(:, p), Ts, false);
end

%% Compute kernel solution.
disp('Computing kernel based solution...');

% Set of admissible control actions.
ur = linspace(-1, 1, 100);

tic 

% Compute the cost function for our samples. 
c = vecnorm(ys);

% Compute the results of the kernel based algorithm.
alg = KernelControl('Sigma', 1, 'Lambda', 1/(M^2));
results = alg.compute(xs, us, xt, ur, c);
U_alg = results.u_opt;

toc

yt_alg = zeros(size(xt));
for p = 1:size(xt, 2)
    yt_alg(:, p) = dynamics_di(xt(:, p), U_alg(:, p), Ts, false);
end

%% Plot the results.

figure;
ax = axes; 
ax.NextPlot = 'add';
ax.Units = 'points';

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = '$$x_1$$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$$x_2$$';
set(ax, 'FontSize', 8);

xlim([-0.75, 0.75])
ylim([-0.75, 0.75])

D = yt_opt - xt;
quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));

D = yt_alg - xt;
quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));

% Save the figure as 'figure1'.
saveas(gcf, '../results/figure1.png')
savefig('../results/figure1.fig')

