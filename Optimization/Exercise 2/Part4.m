clc;
clear;

function run_levenberg()
    syms x y
    objective_func = x^5 * exp(-x^2 - y^2);
    gradient_func = matlabFunction(gradient(objective_func, [x, y]), 'Vars', [x, y]);
    hessian_func = matlabFunction(hessian(objective_func, [x, y]), 'Vars', [x, y]);
    eval_func = matlabFunction(objective_func, 'Vars', [x, y]);
    e = 1e-6;
    gamma = 0.1;
    b = 10;
    a = 0;
    l = 1e-4;
    start_points = [0, 0; -1, 1; 1, -1]';
    
    methods = {@levenberg_marquardt_fixed, @levenberg_armijo, @levenberg_fminbnd, @levenberg_golden_section};
    method_names = {"Fixed Gamma", "Armijo Rule", "fminbnd", "Golden Section"};

    for i = 1:size(start_points, 2)
        plot_results_for_start(eval_func, gradient_func, hessian_func, methods, method_names, start_points(:, i), e, gamma, a, b, l, i);
    end
end

function plot_results_for_start(f, grad, hess, methods, method_names, start_point, epsilon, gamma, a, b, l, figure_index)
    figure('Name', sprintf('Levenberg Results for Starting Point (%.2f, %.2f)', start_point(1), start_point(2)));
    num_methods = length(methods); 
    rows = 2; 
    cols = num_methods; 
    
    for j = 1:length(methods)
        meth = methods{j};
        if j == 1
            [xk, n, f_values] = meth(f, grad, hess, epsilon, gamma, start_point);
        elseif j == 2
            [xk, n, f_values] = meth(f, grad, hess, epsilon, gamma, start_point);
        elseif j == 3
            [xk, n, f_values] = meth(f, grad, hess, epsilon, b, start_point);
        else
            [xk, n, f_values] = meth(f, grad, hess, epsilon, a, b, l, start_point);
        end

        subplot(rows, cols, j);
        [X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));
        Z = arrayfun(f, X, Y);
        contour(X, Y, Z, 20);
        hold on;
        plot(xk(1, :), xk(2, :), '-o', 'DisplayName', method_names{j});
        scatter(xk(1, end), xk(2, end), '*', 'DisplayName', 'Minimum');
        title(sprintf('%s\n# Iter = %d, Min = (%.2f, %.2f)', method_names{j}, n, xk(1, end), xk(2, end)));
        xlabel('x');
        ylabel('y');
        legend show;
        hold off;

        subplot(rows, cols, j + cols);
        plot(0:n, f_values, '-o');
        hold on;
        scatter(n, f_values(end), '*');
        title(sprintf('f(x_k) vs Iterations\n%s', method_names{j}));
        xlabel('Iteration (k)');
        ylabel('f(x_k, y_k)');
        grid on;
        hold off;
    end
end

function [xk, n, f_values] = levenberg_marquardt_fixed(f, grad, hess, epsilon, gamma, x0)
    xk = x0;
    k = 1;
    f_values = f(xk(1), xk(2));
    mu = 0;
    while norm(grad(xk(1, k), xk(2, k))) > epsilon
        grad_k = grad(xk(1, k), xk(2, k));
        hess_matrix = hess(xk(1, k), xk(2, k));
        eigen_values = eig(hess_matrix);
        if any(eigen_values < 0)
            mu = abs(max(eigen_values)) + 1;
        else
            mu = 0;
        end
        A = hess_matrix + mu * eye(size(hess_matrix));
        d_k = -A \ grad_k;
        xk_next = xk(:, k) + gamma * d_k;
        xk = [xk, xk_next];
        f_values = [f_values, f(xk_next(1), xk_next(2))];
        k = k + 1;
        if k > 40000
          break;
        end
    end
    n = k - 1;
end

function [xk, n, f_values] = levenberg_armijo(f, f_grad, f_hess, epsilon, gamma, x0)
    alpha = 1e-3;
    beta = 0.5;
    xk = x0;
    k = 1;
    mk = 0;
    s = 5;
    mu = 0;
    f_values = f(xk(1), xk(2));
    while norm(f_grad(xk(1, k), xk(2, k))) > epsilon
        grad_k = f_grad(xk(1, k), xk(2, k));
        hess_k = f_hess(xk(1, k), xk(2, k));
        eigen_values = eig(hess_k);
        if any(eigen_values < 0)
            mu = abs(max(eigen_values)) + 1;
        else
            mu = 0;
        end
        hess_k = hess_k + mu * eye(size(hess_k));
        dk = -hess_k \ grad_k;
        fx_k = f(xk(1, k), xk(2, k));
        while fx_k - f(xk(1, k) + gamma * dk(1), xk(2, k) + gamma * dk(2)) < -alpha * gamma * (dk' * grad_k)
            mk = mk + 1;
            gamma = s * beta^mk;
        end
        x_next = xk(:, k) + gamma * dk;
        xk = [xk, x_next];
        f_values = [f_values, f(x_next(1), x_next(2))];
        k = k + 1;
        if k > 40000
            break;
        end
        mk = 0;
        gamma = s * beta^mk;
    end
    n = k - 1;
end

function [xk, n, f_values] = levenberg_fminbnd(f, f_grad, f_hess, epsilon, b, x0)
    xk = x0;
    k = 1;
    f_values = f(xk(1), xk(2));
    mu = 0;
    while norm(f_grad(xk(1, k), xk(2, k))) > epsilon
        grad_k = f_grad(xk(1, k), xk(2, k));
        hess_k = f_hess(xk(1, k), xk(2, k));
        eigen_values = eig(hess_k);
        if any(eigen_values < 0)
            mu = abs(max(eigen_values)) + 1;
        else
            mu = 0;
        end
        hess_k = hess_k + mu * eye(size(hess_k));
        dk = -hess_k \ grad_k;
        step_func = @(gamma) f(xk(1, k) + gamma * dk(1), xk(2, k) + gamma * dk(2));
        optimal_gamma = fminbnd(step_func, 0, b);
        x_next = xk(:, k) + optimal_gamma * dk;
        xk = [xk, x_next];
        f_values = [f_values, f(x_next(1), x_next(2))];
        k = k + 1;
        if k > 40000
          break;
        end
    end
    n = k - 1;
end

function [xk, n, f_values] = levenberg_golden_section(f, f_grad, f_hess, epsilon, a, b, l, x0)
    xk = x0;
    k = 1;
    f_values = f(xk(1), xk(2));
    mu = 0;
    while norm(f_grad(xk(1, k), xk(2, k))) > epsilon
        grad_k = f_grad(xk(1, k), xk(2, k));
        hess_k = f_hess(xk(1, k), xk(2, k));
        eigen_values = eig(hess_k);
        if any(eigen_values < 0)
            mu = abs(max(eigen_values)) + 1;
        else
            mu = 0;
        end
        hess_k = hess_k + mu * eye(size(hess_k));
        dk = -hess_k \ grad_k;
        step_func = @(gamma) f(xk(1, k) + gamma * dk(1), xk(2, k) + gamma * dk(2));
        [~, ~, ~, min_vals] = golden_section_method(step_func, a, b, l);
        optimal_gamma = min_vals(end);
        x_next = xk(:, k) + optimal_gamma * dk;
        xk = [xk, x_next];
        f_values = [f_values, f(x_next(1), x_next(2))];
        k = k + 1;
        if k > 40000
            break;
        end
    end
    n = k - 1;
end

function [k, a_vals, b_vals, min_vals] = golden_section_method(f, a, b, l)
    a_vals = [];
    b_vals = [];
    k = 0;
    gamma = 0.618;
    one_minus_gamma = 1 - gamma;
    x1 = a + one_minus_gamma * (b-a);
    x2 = a + gamma * (b-a);
    f1 = f(x1);
    f2 = f(x2);

    while abs(b-a) > l
        k = k + 1;
        
        a_vals = [a_vals, a];
        b_vals = [b_vals, b];
      
        if f1 < f2
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + one_minus_gamma * (b-a);
            f1 = f(x1);
        else
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + gamma * (b-a);
            f2 = f(x2);
        end
    end
    min_vals = (a_vals+b_vals) / 2;
end


run_levenberg();
