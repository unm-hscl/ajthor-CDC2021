%% Generate Nonholonomic vehicle sample.
% In this file, we generate a sample of observations taken from a system
% with double integrator dynamics. 
% 
% Note that this file is not run during a reproducible run. In order to
% re-generate the sample, this file must be run manually.
% 

rng(0);

%% Define constants.

% Sampling time.
Ts = 0.25;

% Sample size.
M = 1600;

%% Generate the sample. 

xs = [
    -1 + 2*rand(1, M);
    -1 + 2*rand(1, M)
    ];

us = -1.1 + 2.2*rand(1, size(xs, 2));

% % System matrices (double integrator).
% A = [1 Ts; 0 1];
% B = [(Ts^2)/2; Ts];
% 
% % Compute the disturbance.
% W = 0.01*randn(size(xs));
% 
% % Propagate the dynamics forward.
% ys = A*xs + B*us + W;

% Propagate the dynamics forward.
ys = zeros(size(xs));

for p = 1:size(xs, 2)
    ys(:, p) = dynamics_di(xs(:, p), us(:, p), Ts);
end

%% Plot the sample. For debugging.
% Plot the samples under each control input. "Hair" plot.
figure
ax = axes; ax.NextPlot = 'add';
scatter(xs(1, :), xs(2, :))
scatter(ys(1, :), ys(2, :))
for p = 1:size(xs, 2)
    plot([xs(1, p), ys(1, p)], [xs(2, p), ys(2, p)], 'k');
end

%% Save the sample.
% We save the sample to a file in order to reuse the sample across multiple
% experiments. 

save('../data/sample_di.mat', 'Ts', 'M', 'xs', 'us', 'ys');

%% Generate a sample with a Beta disturbance.

% Propagate the dynamics forward.
ys = zeros(size(xs));

for p = 1:size(xs, 2)
    ys(:, p) = dynamics_di_beta(xs(:, p), us(:, p), Ts);
end

%% Save the sample.

save('../data/sample_di_beta.mat', 'Ts', 'M', 'xs', 'us', 'ys');

%% Generate a sample with an Exponential disturbance.

% Propagate the dynamics forward.
ys = zeros(size(xs));

for p = 1:size(xs, 2)
    ys(:, p) = dynamics_di_exp(xs(:, p), us(:, p), Ts);
end

%% Save the sample.

save('../data/sample_di_exp.mat', 'Ts', 'M', 'xs', 'us', 'ys');
