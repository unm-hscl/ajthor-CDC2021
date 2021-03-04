%% Figure (3)
% Figure showing the phase portrait of a nonholonomic vehicle.

rng(0);

m = 20;
M = m^3;

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
    pi*ones(1, Mt^2);
    ];

clear XX YY

r1 = linspace(0, 0.25, 5);
r2 = linspace(-pi/4, pi/4, 20);
[U1, U2] = meshgrid(r1, r2);

ut = [
    reshape(U1, 1, []);
    reshape(U2, 1, [])
    ];

clear U1 U2

%% Compute kernel solution.
disp('Computing kernel based solution...');

% Kernel parameters.
lambda = 1/(M^2);
sigma = 1;

tic

% Compute the cost function.
% f = vecnorm(ys([1, 2], :)).^2;
f = vecnorm(ys([1, 2], :) - [0.5; 0]).^2;

% Idx = any(abs(ys([1, 2], :)) >= 0.75);
% f(Idx) = f(Idx) + 10;

x_obs = [0; 0];
r_obs = [0.25];

% for k = 1:5
%     x_obs = [x_obs, -0.75 + 1.5*rand(2, 1)];
%     r_obs = [r_obs, 0.5*rand(1)];
    
    f = f + 10*double(vecnorm(ys([1, 2], :) - x_obs(:, end)).^2 <= r_obs(end));
% end

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

%%

x0 = [-0.5; 0; pi/2];
xt = [];
yt = [];

for k = 1:100
    
    xt = [xt, x0];

    Phi = rbf_kernel(xs, xt(:, end), sigma);
    beta = W*(Phi.*Ups); %#ok<MINV>

    a = f*beta;

    cvx_begin quiet
        variable w(100);
        minimize( dot(a, w) );
        subject to
            sum(w) == 1;
            0 <= w;
    cvx_end
    
    [~, Idx] = max(w);
    U = ut(:, Idx);
    
    x0 = nh_dynamics(xt(:, end), U, Ts);
    yt = [yt, x0];
end

toc

%% Plot the results.

figure 
ax = axes;
ax.NextPlot = 'add';
ax.Units = 'points';

plot(xt(1, :), xt(2, :))

for k = 1:size(x_obs, 2)
    pos = [
        x_obs(1, k) - r_obs(k)/2, 
        x_obs(2, k) - r_obs(k)/2, 
        r_obs(k), 
        r_obs(k)
        ];
    rectangle('Position', pos, 'Curvature', [1 1], ...
        'EdgeColor', [0.6350, 0.0780, 0.1840], 'LineWidth', 1)
end

% figure('Units', 'points', ...
%        'Position', [0, 0, 500, 500]);
% ax = axes; 
% ax.NextPlot = 'add';
% ax.Units = 'points';
% 
% % ax.Position = [30, 30, 200, 200];
% 
% ax.XLabel.Interpreter = 'latex';
% ax.XLabel.String = '$$x_1$$';
% ax.YLabel.Interpreter = 'latex';
% ax.YLabel.String = '$$x_2$$';
% set(ax, 'FontSize', 8);
% 
% xlim([-0.75, 0.75])
% ylim([-0.75, 0.75])
% 
% % scatter(xt(1, :), xt(2, :))
% % scatter(yt(1, :), yt(2, :))
% % for k = 1:size(xt, 2)
% %     plot([xt(1, k), yt(1, k)], [xt(2, k), yt(2, k)], 'k');
% % end
% 
% D = yt - xt;
% quiver(xt(1, :), xt(2, :), D(1, :), D(2, :));