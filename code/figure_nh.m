%% Figure (3)
% Figure showing the phase portrait of a nonholonomic vehicle.

rng(0);

% Number of test points.
Mt = 21;

% Range of the test points.
lb = -0.5;
ub = 0.5;
r1 = linspace(lb, ub, Mt);
r2 = linspace(lb, ub, Mt);
[XX, YY] = meshgrid(r1, r2);

xt = [
    reshape(XX, 1, []);
    reshape(YY, 1, []);
    pi/4*ones(1, Mt^2);
    ];

% ur = linspace(0, 1, 100);
r1 = 0:0.1:1;
r2 = -10:1:10;
[U1, U2] = meshgrid(r1, r2);
ur = [
    reshape(U1, 1, []);
    reshape(U2, 1, []);
    ];

%% Compute kernel solution.
disp('Computing kernel based solution...');

% Kernel parameters.
lambda = 1/M^2;
sigma = 1;

tic 

% Compute the cost function for our samples. 
% c = vecnorm(ys);
% c = abs(ys(3, :));
% c = vecnorm(ys([1 2], :)).^2;
c = abs(ys(1, :) - ys(2, :)).^2;

% Compute the results of the kernel based algorithm.
alg = KernelControl('Lambda', lambda, 'Sigma', sigma);
results = alg.compute(xs, us, xt, ur, c);
U = results.u_opt;

toc

% U = zeros(1, size(xt, 2));

yt = zeros(size(xt));

for k = 1:size(xt, 2)
    yt(:, k) = nh_dynamics(xt(:, k), U(:, k), Ts);
end

%% Plot the results.

figure
% figure('Units', 'points', ...
%        'Position', [0, 0, 180, 180]);
ax = axes; 
ax.NextPlot = 'add';
ax.Units = 'points';

% ax.Position = [30, 30, 200, 200];

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = '$$x_1$$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$$x_2$$';
set(ax, 'FontSize', 8);

% xlim([-1, 1])
% ylim([-1, 1])

D = yt - xt;
quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));

% X = [-0.25; 0.25; 0; 0];
% for k = 1:20
%     X = [X, cp_dynamics(X(:, end), -1, Ts)];
% end
% plot(X(1, :), X(2, :))
