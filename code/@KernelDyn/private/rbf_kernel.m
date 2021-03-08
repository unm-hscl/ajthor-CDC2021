function G = rbf_kernel(X, Y, sigma)
% RBF_KERNEL Compute kernel matrix.

    M = size(X, 2);
    T = size(Y, 2);

    G = zeros(M, T);

    for k = 1:size(X, 1)
        G = G + (repmat(Y(k, :), [M, 1]) - repmat(X(k, :)', [1, T])).^2;
    end

    G = exp(-G/(2*sigma^2));

end