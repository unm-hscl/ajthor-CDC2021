function results = compute(~, dynamics, xt, ur)
%COMPUTE Summary of this function goes here
%   Detailed explanation goes here

lb = ur(1);
ub = ur(2);

assert(lb < ub);

% Variable to hold results corresponding to each test point.
U = zeros(size(ur, 1), size(xt, 2));

t_start = tic;

% Choose one test point at a time.
for k = 1:size(xt, 2)

    % Optimization.
    cvx_begin quiet
        variable v(1);
        minimize( norm(dynamics(xt(:, k), v)) );
        subject to 
            lb <= v <= ub; %#ok<*CHAIN,*VUNUS>
    cvx_end

    % Store the results.
    U(:, k) = v;
end

t_elapsed = toc(t_start);

results.u_opt = U;
results.comp_time = t_elapsed;

end

