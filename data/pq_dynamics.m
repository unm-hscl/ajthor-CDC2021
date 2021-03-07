function y = pq_dynamics(x0, u, time)

function dxdt = pq(~, x)
    
    I = 0.25;
    r = 1;
    m = 1;
    g = 9.8;
    
    dxdt =[
        x(2);
        -(1/m)*(u(1) + u(2))*sin(x(5));     % + 0.1*randn(1);
        x(4);
        (1/m)*(u(1) + u(2))*cos(x(5)) - g;  % + 0.1*randn(1);
        x(6);
        (r/I)*(u(1) - u(2));
    ];
end

[~, y] = ode45(@pq, [0, time], x0);

y = y(size(y, 1), :).';

end