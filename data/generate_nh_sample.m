%% Generate Nonholonomic vehicle sample.
% In this file, we generate a sample of observations taken from a system
% with nonholonomic vehicle dynamics. We choose the states uniformly across
% the range [-1.1, 1.1] x [-1.1, 1.1] x [-2*pi, 2*pi] to fully capture the
% area of interest, and then apply random control actions at each state. 
% 
% Note that this file is not run during a reproducible run. In order to
% re-generate the sample, this file must be run manually.
% 

rng(0);

%% Define constants.

% Sampling time.
Ts = 0.1;

% Sample size.
M = 1600;

%% Generate the sample. 

r = linspace(-1.1, 1.1, 12);
t = linspace(-2*pi, 2*pi, 12);
[XX, YY, TT] = ndgrid(r, r, t);

xs = [
    reshape(XX, 1, []);
    reshape(YY, 1, []);
    reshape(TT, 1, []);
    ];

% xs = [
%     -1.1 + 2.2*rand(1, M);
%     -1.1 + 2.2*rand(1, M);
%     -2*pi + 2*pi*rand(1, M)
%     ];

% Cleanup.
clear XX YY TT

us = [
    -0.1 + 1.2*rand(1, size(xs, 2));
    -10.1 + 20.2*rand(1, size(xs, 2));
    ];

% Propagate the dynamics forward.
ys = zeros(size(xs));

for p = 1:size(xs, 2)
    ys(:, p) = dynamics_nh(xs(:, p), us(:, p), Ts);
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

save('../data/sample_nh.mat', 'Ts', 'M', 'xs', 'us', 'ys');
