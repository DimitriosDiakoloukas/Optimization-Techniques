clc;
clear;

function run_steepest_descent()
    syms x y
    f = (1/3)*x^2 + 3*y^2;
    grad_f = matlabFunction(gradient(f, [x, y]), 'Vars', [x, y]);
    eval_f = matlabFunction(f, 'Vars', [x, y]);

    epsilon = 0.01;
    constraints = struct('x_min', -10, 'x_max', 5, 'y_min', -8, 'y_max', 12);

    test_cases = {
        struct('start_point', [5, -5], 'gamma', 0.5, 'step_size', 5, 'name', '\gamma_k = 0.5, s_k = 5'),
        struct('start_point', [-5, 10], 'gamma', 0.1, 'step_size', 2.5, 'name', '\gamma_k = 0.1, s_k = 2.5'),
        struct('start_point', [-5, 10], 'gamma', 0.1, 'step_size', 15, 'name', '\gamma_k = 0.1, s_k = 15'),
        struct('start_point', [8, -10], 'gamma', 0.2, 'step_size', 0.1, 'name', '\gamma_k = 0.2, s_k = 0.1')
    };

    for i = 1:length(test_cases)
        case_data = test_cases{i};
        [xk, k, x_vals, y_vals, f_vals] = steepest(eval_f, grad_f, case_data.start_point, epsilon, case_data.gamma, case_data.step_size, constraints);
        create_subplot(eval_f, x_vals, y_vals, f_vals, k, case_data.name, case_data.start_point);
    end
end

function [xk, k, x_vals, y_vals, f_vals] = steepest(f, grad_f, start_point, epsilon, gamma, step_size, constraints)
    xk = start_point;
    k = 0;
    x_vals = [xk(1)];
    y_vals = [xk(2)];
    f_vals = [f(xk(1), xk(2))];

    while true
        grad_eval = grad_f(double(xk(1)), double(xk(2)));
        if norm(grad_eval) <= epsilon || k > 40000
            break;
        end

        dk = -grad_eval;
        a = xk(1) + dk(1) * step_size;
        b = xk(2) + dk(2) * step_size;

        if (a <= constraints.x_max && a >= constraints.x_min)
            x_next(1) = a;
        elseif (a < constraints.x_min)
            x_next(1) = constraints.x_min;
        elseif (a > constraints.x_max)
            x_next(1) = constraints.x_max;
        end

        if (b <= constraints.y_max && b >= constraints.y_min)
            x_next(2) = b;
        elseif (b < constraints.y_min)
            x_next(2) = constraints.y_min;
        elseif (b > constraints.y_max)
            x_next(2) = constraints.y_max;
        end

        xk(1) = xk(1) + gamma * (x_next(1) - xk(1));
        xk(2) = xk(2) + gamma * (x_next(2) - xk(2));

        x_vals = [x_vals, xk(1)];
        y_vals = [y_vals, xk(2)];
        f_vals = [f_vals, f(double(xk(1)), double(xk(2)))];
        k = k + 1;
    end
end

function create_subplot(f, x_vals, y_vals, f_vals, k, case_name, start_point)
    [x_grid, y_grid] = meshgrid(-10:0.1:10, -10:0.1:10);
    z_vals = arrayfun(f, x_grid, y_grid);

    figure;
    sgtitle(sprintf('Results for %s', case_name));

    subplot(1, 3, 1);
    contour(x_grid, y_grid, z_vals, 30);
    hold on;
    plot(x_vals, y_vals, '-o', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');
    scatter(x_vals(1), y_vals(1), 'filled', 'DisplayName', 'Start Point');
    scatter(x_vals(end), y_vals(end), 'filled', 'DisplayName', 'End Point');
    title(sprintf('Contour Plot\nStart = (%.2f, %.2f), Iter = %d', start_point(1), start_point(2), k));
    xlabel('x'); ylabel('y');
    legend;
    grid on;

    subplot(1, 3, 2);
    plot(1:k+1, f_vals, '-o', 'LineWidth', 1.5);
    title('Convergence of f(x)');
    xlabel('Iteration'); ylabel('f(x)');
    grid on;

    subplot(1, 3, 3);
    plot(1:k+1, x_vals, '-o', 'LineWidth', 1.5, 'DisplayName', 'x1');
    hold on;
    plot(1:k+1, y_vals, '-o', 'LineWidth', 1.5, 'DisplayName', 'x2');
    title('Trajectory of x1 and x2');
    xlabel('Iteration'); ylabel('Value');
    legend;
    grid on;
end

run_steepest_descent();