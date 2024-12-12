function J = lagrangian(D, R, step_size, param)
    lambda = param * step_size^2;
    J = D + lambda * R;
end

