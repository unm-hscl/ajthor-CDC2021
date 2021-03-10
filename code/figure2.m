%% Figure (2)
% Figure showing the phase portrait of the double integrator system under a
% beta distribution disturbance and exponential distribution disturbance.

rng(0);

% Sampling time.
T = 0.25;
% System matrices (double integrator).
A = [1 T; 0 1];
B = [(T^2)/2; T];

% Number of samples.
m = 11;
M = m^2;

% Specify the points where we wish to evaluate at. 

% Number of test points.
Mt = 11;
% Number of control inputs to evaluate at. 
n = 10;

% Range of the test points.
r = linspace(-1, 1, Mt);
[XX, YY] = meshgrid(r, r);

xt = [
    reshape(XX, 1, []);
    reshape(YY, 1, []);
    ];

% Cleanup.
clear XX YY

ut = linspace(-1, 1, n);

% Kernel parameters.
lambda = 1/(M^2);
sigma = 1;

%% Compute the kernel solution (beta disturbance).

% Generate the sample used by the kernel regression. 
% [XX, YY] = meshgrid(linspace(-1, 1, 75), linspace(-1, 1, 75));
% 
% xs = [
%     reshape(XX, 1, []);
%     reshape(YY, 1, []);
%     ];
xs = [
    -1 + 2*rand(1, 5625);
    -1 + 2*rand(1, 5625)
    ];

us = -1.1 + 2.2*rand(1, size(xs, 2));

% Compute the disturbance.
W = 0.1.*betarnd(2, 0.5, size(X));

% Propagate the dynamics forward.
ys = A*xs + B*us + W;

% Cleanup.
clear XX YY UU

tic

% Compute the cost function for our samples. 

c = vecnorm(ys);

% Compute the kernel (Gram) matrix.
Gx = rbf_kernel(xs, xs, sigma);
Gu = rbf_kernel(us, us, sigma);
G = Gx.*Gu;
clear Gx Gu

% Compute the "weight" matrix. (G + lambda*M*I)^-1
K = (G + lambda*M*eye(size(G)));
W = inv(K);

% Compute the "kernel" matrix of input samples vs. input points.
Ups = rbf_kernel(us, ut, sigma);

% Variable to hold results corresponding to each test point.
U = [];

% **BAD** Loop.
% Choose one test point at a time.
for X = xt

    % Compute the "kernel" matrix of samples vs. the test point.
    Phi = rbf_kernel(xs, X, sigma);
    beta = W*(Phi.*Ups); %#ok<MINV>

    % Evaluate the cost function at the test point under each control
    % input so that we get an approximation of the cost f under each
    % input.
    a = c*beta;

    % Optimization. Minimize the inner product by choosing weights.
    cvx_begin quiet
        variable w(n);
        minimize( dot(a, w) );
        subject to
            sum(w) == 1;
            0 <= w;
    cvx_end

    % Store the results.
    [~, Idx] = max(w);
    U = [U, ut(Idx)]; %#ok<AGROW>

end

toc

yt = A*xt + B*U;

%% Plot the results.

figure('Units', 'points', ...
       'Position', [0, 0, 120, 120]);
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

% scatter(xt(1, :), xt(2, :))
% scatter(yt(1, :), yt(2, :))
% for k = 1:size(xt, 2)
%     plot([xt(1, k), yt(1, k)], [xt(2, k), yt(2, k)], 'k');
% end

D = yt - xt;
quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));

%% Compute the kernel solution (exponential disturbance).

% Generate the sample used by the kernel regression. 
% [XX, YY] = meshgrid(linspace(-1, 1, 75), linspace(-1, 1, 75));
% 
% xs = [
%     reshape(XX, 1, []);
%     reshape(YY, 1, []);
%     ];
xs = [
    -1 + 2*rand(1, 5625);
    -1 + 2*rand(1, 5625)
    ];

us = -1.1 + 2.2*rand(1, size(xs, 2));

% Compute the disturbance.
W = 0.01.*exprnd(3, size(X));

% Propagate the dynamics forward.
ys = A*xs + B*us + W;

% Cleanup.
clear XX YY UU

tic

% Compute the cost function for our samples. 

c = vecnorm(ys);

% Compute the kernel (Gram) matrix.
Gx = rbf_kernel(xs, xs, sigma);
Gu = rbf_kernel(us, us, sigma);
G = Gx.*Gu;
clear Gx Gu

% Compute the "weight" matrix. (G + lambda*M*I)^-1
K = (G + lambda*M*eye(size(G)));
W = inv(K);

% Compute the "kernel" matrix of input samples vs. input points.
Ups = rbf_kernel(us, ut, sigma);

% Variable to hold results corresponding to each test point.
U = [];

% **BAD** Loop.
% Choose one test point at a time.
for X = xt

    % Compute the "kernel" matrix of samples vs. the test point.
    Phi = rbf_kernel(xs, X, sigma);
    beta = W*(Phi.*Ups); %#ok<MINV>

    % Evaluate the cost function at the test point under each control
    % input so that we get an approximation of the cost f under each
    % input.
    a = c*beta;

    % Optimization. Minimize the inner product by choosing weights.
    cvx_begin quiet
        variable w(n);
        minimize( dot(a, w) );
        subject to
            sum(w) == 1;
            0 <= w;
    cvx_end

    % Store the results.
    [~, Idx] = max(w);
    U = [U, ut(Idx)]; %#ok<AGROW>

end

toc

yt = A*xt + B*U;

%% Plot the results.

figure('Units', 'points', ...
       'Position', [0, 0, 120, 120]);
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

% scatter(xt(1, :), xt(2, :))
% scatter(yt(1, :), yt(2, :))
% for k = 1:size(xt, 2)
%     plot([xt(1, k), yt(1, k)], [xt(2, k), yt(2, k)], 'k');
% end

D = yt - xt;
quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));
