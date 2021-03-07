Ts = 0.1;

M = 5625;

xs = [
    -1.1 + 2.2*rand(1, M);
    -1.1 + 2.2*rand(1, M);
    -1.1 + 2.2*rand(1, M);
    -1.1 + 2.2*rand(1, M);
    -1.1 + 2.2*rand(1, M);
    -1.1 + 2.2*rand(1, M);
    ];
    
% xs = [
%     -10 + 20*rand(1, M);
%     -10 + 20*rand(1, M);
%     -10 + 20*rand(1, M);
%     -10 + 20*rand(1, M);
%     -2*pi + 4*pi*rand(1, M);
%     -2*pi + 4*pi*rand(1, M);
%     ];

us = [
    -0.1 + 10.2*rand(1, M);
    -0.1 + 10.2*rand(1, M);
    ];

% r1 = linspace(-1.1, 1.1, 10);
% r2 = linspace(-0.5, 0.5, 5);
% r3 = linspace(-1.1, 1.1, 10);
% r4 = linspace(-0.5, 0.5, 5);
% r5 = linspace(-pi, pi, 10);
% r6 = linspace(-pi, pi, 5);
% 
% [XX, DX, YY, DY, TT, DT] = ndgrid(r1, r2, r3, r4, r5, r6);
% 
% xs = [
%     reshape(XX, 1, []);
%     reshape(DX, 1, []);
%     reshape(YY, 1, []);
%     reshape(DY, 1, []);
%     reshape(TT, 1, []);
%     reshape(DT, 1, []);
%     ];

% Propagate the dynamics forward.
% ys = pq_dynamics(xs, us, Ts);
ys = zeros(size(xs));

for k = 1:size(xs, 2)
    ys(:, k) = pq_dynamics(xs(:, k), us(:, k), Ts);
end

% Plot the samples under each control input. "Hair" plot.
figure
ax = axes; ax.NextPlot = 'add';
scatter(xs(1, :), xs(3, :))
scatter(ys(1, :), ys(3, :))
for k = 1:size(xs, 2)
    plot([xs(1, k), ys(1, k)], [xs(3, k), ys(3, k)], 'k');
end
