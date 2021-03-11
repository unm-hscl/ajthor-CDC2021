function y = dynamics_nh(x0, u, time)
% DYNAMICS_NH Nonholonomic vehicle dynamics.

function dxdt = nh(~, x)

    dxdt =[
        (u(1) + 0.1)*sin(x(3)) + 0.1*randn(1);
        (u(1) + 0.1)*cos(x(3)) + 0.1*randn(1);
        u(2) + 0.1*randn(1);
    ];

end

[~, y] = ode45(@nh, [0, time], x0);

y = y(size(y, 1), :).';

% We constrain the angle of the system to be within the range [0, 2*pi].
if abs(y(3)) >= 2*pi
    y(3) = mod(y(3), 2*pi);
end

end