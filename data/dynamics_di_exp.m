function y = dynamics_di_exp(x0, u, time)
% DYNAMICS_DI_EXP Double integrator dynamics.

function dxdt = nh(~, x)

    dxdt =[
        x(2) + 0.01.*exprnd(3, 1);
        u + 0.01.*exprnd(3, 1);
    ];

end

[~, y] = ode45(@nh, [0, time], x0);

y = y(size(y, 1), :).';

end

