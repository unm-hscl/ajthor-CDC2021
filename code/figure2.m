%% Figure (2)
% Figure showing the phase portrait of the double integrator system under a
% beta distribution disturbance and exponential distribution disturbance.

rng(0);

% Create a grid of test points.
r = linspace(-1, 1, 11);
[XX, YY] = meshgrid(r, r);

xt = [
    reshape(XX, 1, []);
    reshape(YY, 1, []);
    ];

% Cleanup.
clear XX YY

%% Compute the kernel solution (beta disturbance).
disp('Computing kernel based solution (beta)...');

% Set of admissible control actions.
ur = linspace(-1, 1, 100);

% Load the sample.
load('sample_di_beta.mat');

tic 

% Compute the cost function for our samples. 
c = vecnorm(ys);

% Compute the results of the kernel based algorithm.
alg = KernelControl('Sigma', 1, 'Lambda', 1/(M^2));
results = alg.compute(xs, us, xt, ur, c);
U_beta = results.u_opt;

toc

yt_beta = zeros(size(xt));
for p = 1:size(xt, 2)
    yt_beta(:, p) = dynamics_di(xt(:, p), U_beta(:, p), Ts, false);
end

%% Compute the kernel solution (exponential disturbance).
disp('Computing kernel based solution (exponential)...');

% Set of admissible control actions.
ur = linspace(-1, 1, 100);

% Load the sample.
load('sample_di_exp.mat');

tic 

% Compute the cost function for our samples. 
c = vecnorm(ys);

% Compute the results of the kernel based algorithm.
alg = KernelControl('Sigma', 1, 'Lambda', 1/(M^2));
results = alg.compute(xs, us, xt, ur, c);
U_exp = results.u_opt;

toc

yt_exp = zeros(size(xt));
for p = 1:size(xt, 2)
    yt_exp(:, p) = dynamics_di(xt(:, p), U_exp(:, p), Ts, false);
end

%% Plot the results.

figure;
ax1 = axes; 
ax1.NextPlot = 'add';
ax1.Units = 'points';

ax1.XLabel.Interpreter = 'latex';
ax1.XLabel.String = '$$x_1$$';
ax1.YLabel.Interpreter = 'latex';
ax1.YLabel.String = '$$x_2$$';
set(ax1, 'FontSize', 8);

xlim([-0.75, 0.75])
ylim([-0.75, 0.75])

D = yt_beta - xt;
quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));

% Save the figure as 'figure2a'.
saveas(gcf, '../results/figure2a.png')
savefig('../results/figure2a.fig')

%% Plot the results.

figure;
ax2 = axes; 
ax2.NextPlot = 'add';
ax2.Units = 'points';

ax2.XLabel.Interpreter = 'latex';
ax2.XLabel.String = '$$x_1$$';
ax2.YLabel.Interpreter = 'latex';
ax2.YLabel.String = '$$x_2$$';
set(ax2, 'FontSize', 8);

xlim([-0.75, 0.75])
ylim([-0.75, 0.75])

D = yt_exp - xt;
quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));

% Save the figure as 'figure2b'.
saveas(gcf, '../results/figure2b.png')
savefig('../results/figure2b.fig')
