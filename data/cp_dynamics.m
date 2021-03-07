function y = cp_dynamics(x0, u, time)

g = 9.8;
m = 0.1;
l = 0.5;
mt = 1.1;

function dxdt = cp(~, x)
    
    A = (u + m*l*(x(4)^2)*sin(x(3)))/mt;
    B = (l*((4/3) - m*(cos(x(3))^2)/mt)*mt);
    
    dxdt =[
        x(2);
        A - (m*l*(g*sin(x(3)) - cos(x(3)))*A*cos(x(3)))/B;
        x(4);
        (g*sin(x(3)) - cos(x(3))*A*cos(x(3)))/B;
    ];

%     dxdt =[
%         x(2);
%         -x(1) + 0.1*sin(x(3));
%         x(4);
%         u;
%     ];

%     dxdt =[
%         x(2);
%         0.0043*x(4) - 2.75*x(3) + 1.94*u - 10.95*x(2);
%         x(3);
%         28.58*x(3) - 0.044*x(4) - 4.44*u + 24.92*x(2);
%     ];

end

[~, y] = ode45(@cp, [0, time], x0);

y = y(size(y, 1), :).';

end