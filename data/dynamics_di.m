function y = dynamics_di(x0, u, time, varargin)
% DYNAMICS_DI Double integrator dynamics.

function dxdt = nh(~, x)

    if ~isempty(varargin) && varargin{1}
        dxdt =[
            x(2);
            u;
        ];
    else
        dxdt =[
            x(2) + 0.01*randn(1);
            u + 0.01*randn(1);
        ];
    end

end

[~, y] = ode45(@nh, [0, time], x0);

y = y(size(y, 1), :).';

end

