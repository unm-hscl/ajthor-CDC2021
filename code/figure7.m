%% Figure (7)
% Figure showing the backward in time trajectory of a nonholonomic vehicle.

rng(0);

% Load the sample.
load('../data/sample_nh.mat');

% Time horizon.
N = 25;

% Set of admissible control actions.
r1 = 0:0.1:1;
r2 = -10:1:10;
[U1, U2] = meshgrid(r1, r2);
ur = [
    reshape(U1, 1, []);
    reshape(U2, 1, []);
    ];

%% Define the target trajectory.

r = linspace(-1, 1, N);
R = [r; r; pi/4*ones(size(r))];

%% Compute kernel solution.
disp('Computing backward in time target tracking trajectory...');

% Specify the initial condition.
x0 = [-0.8; 0; pi];

alg = KernelDyn('Sigma', 3, 'Lambda', 1/(M^2));

cost        = @(t) vecnorm(ys([1 2], :) - R([1 2], t), 2).^2;
dynamics    = @(x, u) dynamics_nh(x, u, Ts);

tic 

results = alg.compute(xs, us, ys, x0, ur, N, cost, dynamics);

toc

X = results.x_traj;

%% Plot the results.

figure;
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

plot(R(1, :), R(2, :), 'kx:');
plot(X(1, :), X(2, :), '^--');

% Save the figure as 'figure7'.
saveas(gcf, '../results/figure7.png')
savefig('../results/figure7.fig')
