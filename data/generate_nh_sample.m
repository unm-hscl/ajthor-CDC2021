Ts = 0.1;

M = 1600;

rng(0);

% Generate the sample. 
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

% us = rand(1, size(xs, 2));
us = [
    -0.1 + 1.2*rand(1, size(xs, 2));
    -10.1 + 20.2*rand(1, size(xs, 2));
    ];

% Propagate the dynamics forward.
ys = zeros(size(xs));

for k = 1:size(xs, 2)
    ys(:, k) = nh_dynamics(xs(:, k), us(:, k), Ts);
end

% Plot the samples under each control input. "Hair" plot.
figure
ax = axes; ax.NextPlot = 'add';
scatter(xs(1, :), xs(2, :))
scatter(ys(1, :), ys(2, :))
for k = 1:size(xs, 2)
    plot([xs(1, k), ys(1, k)], [xs(2, k), ys(2, k)], 'k');
end
