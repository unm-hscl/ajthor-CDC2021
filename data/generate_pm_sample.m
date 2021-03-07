
% Number of test points.
M = 5000;

% Sampling time.
Ts = 0.25;

xs = [
    -5 + 10*rand(1, M);
    -5 + 10*rand(1, M);
%     -5 + 10*rand(1, Mt);
%     -5 + 10*rand(1, Mt);
    ];

us = [
    -2 + 4*rand(1, M);
%     -2 + 4*rand(1, Mt);
    ];

% Propagate the dynamics forward.
% ys = pq_dynamics(xs, us, Ts);
ys = zeros(size(xs));

for k = 1:size(xs, 2)
    ys(:, k) = pm_dynamics(xs(:, k), us(:, k), Ts);
end

% Plot the samples under each control input. "Hair" plot.
figure
ax = axes; ax.NextPlot = 'add';
scatter(xs(1, :), xs(2, :))
scatter(ys(1, :), ys(2, :))
for k = 1:size(xs, 2)
    plot([xs(1, k), ys(1, k)], [xs(2, k), ys(2, k)], 'k');
end
