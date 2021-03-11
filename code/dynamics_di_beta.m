function y = dynamics_di_beta(x0, u, time)
% DYNAMICS_DI_BETA Double integrator dynamics.

function dxdt = nh(~, x)

    dxdt =[
        x(2) + 0.1.*betarnd(2, 0.5, 1);
        u + 0.1.*betarnd(2, 0.5, 1);
    ];

end

[~, y] = ode45(@nh, [0, time], x0);

y = y(size(y, 1), :).';

end

