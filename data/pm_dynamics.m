function y = pm_dynamics(x0, u, time)

function dxdt = pq(~, x)

    dxdt =[
        x(2);
        u(1);   % + 0.1*randn(1);
    ];

end

[~, y] = ode45(@pq, [0, time], x0);

y = y(size(y, 1), :).';

end