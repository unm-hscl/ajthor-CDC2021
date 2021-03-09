%% Figure (1)
% Figure showing the phase portrait of the optimal solution vs. the kernel
% based solution.

rng(0);

% Sampling time.
T = 0.25;

% System matrices (double integrator).
A = [1 T; 0 1];
B = [(T^2)/2; T];

% Number of test points.
Mt = 11;

% Range of the test points.
r = linspace(-1, 1, Mt);
[XX, YY] = meshgrid(r, r);

xt = [
    reshape(XX, 1, []);
    reshape(YY, 1, []);
    ];

% Cleanup.
clear XX YY

%% Compute optimal solution.
disp('Computing optimal solution...');
tic

dynamics = @(x, u) A*x + B*u;

tic

% Compute the results using the optimization based algorithm.
alg_opt = OptimalControl();
results_opt = alg_opt.compute(dynamics, xt, [-1, 1]);
U = results_opt.u_opt;

toc

yt = A*xt + B*U;

%% Plot the results.

figure('Units', 'points', ...
       'Position', [0, 0, 240, 240]);
ax = axes; 
ax.NextPlot = 'add';
ax.Units = 'points';

% ax.Position = [30, 30, 200, 200];

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = '$$x_1$$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$$x_2$$';
set(ax, 'FontSize', 8);

xlim([-0.75, 0.75])
ylim([-0.75, 0.75])

D = yt - xt;
quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));

%% Compute kernel solution.
disp('Computing kernel based solution...');

% Number of samples.
M = 5625;

% Kernel parameters.
lambda = 1/(M^2);
sigma = 0.25;

% Generate the sample used by the kernel regression. 
xs = [
    -1 + 2*rand(1, M);
    -1 + 2*rand(1, M)
    ];

us = -1.1 + 2.2*rand(1, size(xs, 2));

% Compute the disturbance.
W = 0.01*randn(size(xs));

% Propagate the dynamics forward.
ys = A*xs + B*us + W;

% Cleanup.
clear W XX YY UU

% Number of control inputs to evaluate at. 
n = 100;

% Set of admissible control actions.
ur = linspace(-1, 1, n);

tic 

% Compute the cost function for our samples. 
c = vecnorm(ys);

% Compute the results of the kernel based algorithm.
alg = KernelControl('Lambda', lambda, 'Sigma', sigma);
results = alg.compute(xs, us, xt, ur, c);
U = results.u_opt;

toc

yt = A*xt + B*U;

%% Plot the results.

D = yt - xt;
quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));

% Save the figure as 'figure1'.
% saveas(gcf, '../results/figure1.png')
