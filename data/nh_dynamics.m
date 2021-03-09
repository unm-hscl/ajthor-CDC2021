function y = nh_dynamics(x0, u, time)

function dxdt = nh(~, x)

%     dxdt =[
%         u(1)*sin(x(3));% + 0.1*randn(1);
%         u(1)*cos(x(3));% + 0.1*randn(1);
%         0.1;
%     ];

    dxdt =[
        (u(1) + 0.1)*sin(x(3));% + 0.1*randn(1);
        (u(1) + 0.1)*cos(x(3));% + 0.1*randn(1);
        u(2);% + 0.1*randn(1);
    ];
end

[~, y] = ode45(@nh, [0, time], x0);

y = y(size(y, 1), :).';

end