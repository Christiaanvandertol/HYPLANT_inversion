function J = numericalJacobian(func, x)
    n = length(x); % Number of inputs
    m = length(func(x)); % Number of outputs
    J = zeros(m, n); % Initialize Jacobian matrix

    delta = 1e-6; % Small perturbation

    for i = 1:n
        x1 = x; 
        x1(i) = x1(i) + delta; % Perturb input
        J(:, i) = (func(x1) - func(x)) / delta; % Compute partial derivatives
    end
end