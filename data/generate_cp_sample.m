Ts = 0.1;

M = 5625;

rng(0);

lb = -0.5 - 0.1;
ub =  0.5 + 0.1;

xs = [
    lb + (ub - lb)*rand(1, M);
    lb + (ub - lb)*rand(1, M);
    lb + (ub - lb)*rand(1, M);
    lb + (ub - lb)*rand(1, M);
    ];

us = -10.1 + 20.2*rand(1, M);
% us = lb + (ub - lb)*rand(1, M);

% Propagate the dynamics forward.
ys = zeros(size(xs));

for k = 1:size(xs, 2)
    ys(:, k) = cp_dynamics(xs(:, k), us(:, k), Ts);
end

% Plot the samples under each control input. "Hair" plot.
figure
ax = axes; ax.NextPlot = 'add';
scatter(xs(1, :), xs(3, :))
scatter(ys(1, :), ys(3, :))
for k = 1:size(xs, 2)
    plot([xs(1, k), ys(1, k)], [xs(3, k), ys(3, k)], 'k');
end
