%% Figure (3)
% Figure showing the maximum and mean error of the kernel based approach.

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

%% Compute optimal solution.
disp('Computing optimal solution...');

Ts = 0.25;
A = [1, Ts; 0 1];
B = [(Ts^2)/2; Ts];
dynamics = @(x, u) A*x + B*u;

tic

% Compute the results using the optimization based algorithm.
alg_opt = OptimalControl();
results_opt = alg_opt.compute(dynamics, xt, [-1, 1]);
ut = results_opt.u_opt;

toc

yt_opt = zeros(size(xt));
for p = 1:size(xt, 2)
    yt_opt(:, p) = dynamics_di(xt(:, p), ut(:, p), Ts, false);
end

%% Compute kernel solution.
disp('Computing kernel based solution (error)...');

% Set of admissible control actions.
ur = linspace(-1, 1, 100);

% Number of samples.
M_range = 100:500:8500;

% Variables to hold the error.
max_error_gaussian  = zeros(size(M_range));
max_error_beta      = zeros(size(M_range));
max_error_exp       = zeros(size(M_range));

mean_error_gaussian = zeros(size(M_range));
mean_error_beta     = zeros(size(M_range));
mean_error_exp      = zeros(size(M_range));

for k = 1:length(M_range)
    
    M = M_range(k);
    
    fprintf('Computing for M = %d\n', M);
    
    alg = KernelControl('Sigma', 1, 'Lambda', 1/(M^2));

    % Generate the sample used by the kernel regression. 
    xs = [
        -1 + 2*rand(1, M);
        -1 + 2*rand(1, M)
        ];

    us = -1.1 + 2.2*rand(1, size(xs, 2));
    
    %%%%% GAUSSIAN %%%%%
    
    ys = zeros(size(xs));

    for p = 1:size(xs, 2)
        ys(:, p) = dynamics_di(xs(:, p), us(:, p), Ts);
    end

    tic 
    
    % Compute the cost function for our samples. 
    c = vecnorm(ys);

    % Compute the results of the kernel based algorithm.
    results = alg.compute(xs, us, xt, ur, c);
    ut = results.u_opt;

    toc
    
    yt_gauss = zeros(size(xt));

    for p = 1:size(xt, 2)
        yt_gauss(:, p) = dynamics_di(xt(:, p), ut(:, p), Ts);
    end
    
    max_error_gaussian(k)   = max(vecnorm(yt_gauss - yt_opt, 2).^2);
    mean_error_gaussian(k)  = mean(vecnorm(yt_gauss - yt_opt, 2).^2);
    
    %%%%% BETA %%%%%
    
    ys = zeros(size(xs));

    for p = 1:size(xs, 2)
        ys(:, p) = dynamics_di_beta(xs(:, p), us(:, p), Ts);
    end

    tic 
    
    % Compute the cost function for our samples. 
    c = vecnorm(ys);

    % Compute the results of the kernel based algorithm.
    results = alg.compute(xs, us, xt, ur, c);
    ut = results.u_opt;

    toc
    
    yt_beta = zeros(size(xt));

    for p = 1:size(xt, 2)
        yt_beta(:, p) = dynamics_di(xt(:, p), ut(:, p), Ts);
    end
    
    max_error_beta(k)   = max(vecnorm(yt_beta - yt_opt, 2).^2);
    mean_error_beta(k)  = mean(vecnorm(yt_beta - yt_opt, 2).^2);
    
    %%%%% EXPONENTIAL %%%%%
    
    ys = zeros(size(xs));

    for p = 1:size(xs, 2)
        ys(:, p) = dynamics_di_exp(xs(:, p), us(:, p), Ts);
    end

    tic 
    
    % Compute the cost function for our samples. 
    c = vecnorm(ys);

    % Compute the results of the kernel based algorithm.
    results = alg.compute(xs, us, xt, ur, c);
    ut = results.u_opt;

    toc
    
    yt_exp = zeros(size(xt));

    for p = 1:size(xt, 2)
        yt_exp(:, p) = dynamics_di(xt(:, p), ut(:, p), Ts);
    end
    
    max_error_exp(k)    = max(vecnorm(yt_exp - yt_opt, 2).^2);
    mean_error_exp(k)   = mean(vecnorm(yt_exp - yt_opt, 2).^2);
    
end

%% Plot the results.

figure;
ax1 = axes; 
ax1.NextPlot = 'add';
ax1.Units = 'points';

ax1.XLabel.Interpreter = 'latex';
ax1.XLabel.String = 'Sample Size $$M$$';
ax1.YLabel.Interpreter = 'latex';
ax1.YLabel.String = 'Max Error';
set(ax1, 'FontSize', 8);

grid on

xlim([0, 8000])

plot(M_range, max_error_gaussian);
plot(M_range, max_error_beta);
plot(M_range, max_error_exp);

legend(...
    'Gaussian Disturbance', ...
    'Beta Disturbance', ...
    'Exponential Disturbance' ...
    );

% Save the figure as 'figure3a'.
saveas(gcf, '../results/figure3a.png')
savefig('../results/figure3a.fig')

figure;
ax = axes; 
ax.NextPlot = 'add';
ax.Units = 'points';

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Sample Size $$M$$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Mean Error';
set(ax, 'FontSize', 8);

grid on

xlim([0, 8000])

plot(M_range, mean_error_gaussian);
plot(M_range, mean_error_beta);
plot(M_range, mean_error_exp);

legend(...
    'Gaussian Disturbance', ...
    'Beta Disturbance', ...
    'Exponential Disturbance' ...
    );

% Save the figure as 'figure3b'.
saveas(gcf, '../results/figure3b.png')
savefig('../results/figure3b.fig')
