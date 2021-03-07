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

y_opt = A*xt + B*U;

%% Compute kernel solution.
disp('Computing kernel based solution...');

% Number of control inputs to evaluate at. 
n = 100;

% Set of admissible control actions.
ur = linspace(-1, 1, n);

% Number of samples.
M_range = 100:500:8500;

%%
% Variable to hold the error.
E = zeros(size(M_range));

for k = 1:length(M_range)
    
    M = M_range(k);

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
    w = 0.01*randn(size(xs));

    % Propagate the dynamics forward.
    ys = A*xs + B*us + w;

    tic 
    
    % Compute the cost function for our samples. 
    c = vecnorm(ys);

    % Compute the results of the kernel based algorithm.
    alg = KernelControl('Lambda', lambda, 'Sigma', sigma);
    results = alg.compute(xs, us, xt, ur, c);
    U = results.u_opt;

    toc

    yt = A*xt + B*U;
    
    E(k) = mean(vecnorm(yt - y_opt, 2).^2);
    
end

%%
% Variable to hold the error.
E_beta = zeros(size(M_range));

for k = 1:length(M_range)
    
    M = M_range(k);

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
    w = 0.1.*betarnd(2, 0.5, size(xs));

    % Propagate the dynamics forward.
    ys = A*xs + B*us + w;

    tic 

    % Compute the cost function for our samples. 
    c = vecnorm(ys);

    % Compute the results of the kernel based algorithm.
    alg = KernelControl('Lambda', lambda, 'Sigma', sigma);
    results = alg.compute(xs, us, xt, ur, c);
    U = results.u_opt;

    toc

    yt = A*xt + B*U;
    
    E_beta(k) = mean(vecnorm(yt - y_opt, 2).^2);
    
end

%%
% Variable to hold the error.
E_exp = zeros(size(M_range));

for k = 1:length(M_range)
    
    M = M_range(k);

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
    w = 0.01.*exprnd(3, size(xs));

    % Propagate the dynamics forward.
    ys = A*xs + B*us + w;

    tic 

    % Compute the cost function for our samples. 
    c = vecnorm(ys);

    % Compute the results of the kernel based algorithm.
    alg = KernelControl('Lambda', lambda, 'Sigma', sigma);
    results = alg.compute(xs, us, xt, ur, c);
    U = results.u_opt;

    toc

    yt = A*xt + B*U;
    
    E_exp(k) = mean(vecnorm(yt - y_opt, 2).^2);
    
end

%% Plot the results.

figure('Units', 'points', ...
       'Position', [0, 0, 240, 120]);
ax = axes; 
ax.NextPlot = 'add';
ax.Units = 'points';

% ax.Position = [30, 30, 200, 200];

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Sample Size $$M$$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Mean Error';
set(ax, 'FontSize', 8);

grid on

xlim([0, 8000])
ylim([0, 0.11])
yticks(0.01:0.02:0.11)

plot(M_range, E);
plot(M_range, E_beta);
plot(M_range, E_exp);

legend(...
    'Gaussian Disturbance', ...
    'Beta Disturbance', ...
    'Exponential Disturbance' ...
    );


% Save the figure as 'figure_err'.
% saveas(gcf, '../results/figure_err.png')
