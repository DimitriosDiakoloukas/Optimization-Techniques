clc;
clear;

function [k, a_vals, b_vals, min_vals] = fibonacci_method(f, a, b, epsilon, N)
    fib = fibonacci_sequence(N);
    a_vals = [];
    b_vals = [];
    min_vals = [];
    x1 = a + (fib(N-2) / fib(N)) * (b-a);
    x2 = a + (fib(N-1) / fib(N)) * (b-a);
    f1 = f(x1);
    f2 = f(x2);
    
    k = 0;
    for i = 1:N-1
        a_vals = [a_vals, a];
        b_vals = [b_vals, b];
        min_vals = [min_vals, (a+b) / 2];
        if i ~= N-2
            if f1 < f2
                b = x2;
                x2 = x1;
                x1 = a + (fib(N-i-2) / fib(N-i)) * (b-a);
                f2 = f1;
                f1 = f(x1);
            else
                a = x1;
                x1 = x2;
                x2 = a + (fib(N-i-1) / fib(N-i)) * (b-a);
                f1 = f2;
                f2 = f(x2);
            end
        else
            x2 = x1 + epsilon;
            f2 = f(x2);
            if f1 < f2
                b = x2;
            else
                a = x1;
            end
            a_vals = [a_vals, a];
            b_vals = [b_vals, b];
            min_vals = [min_vals, (a+b) / 2];
            k = k + 1;
            break;
        end
        k = k + 1;
    end
end

function fib = fibonacci_sequence(N)
    fib = zeros(1, N);
    fib(1) = 1; 
    fib(2) = 1;
    for i = 3:N
        fib(i) = fib(i-1) + fib(i-2);
    end
end

function N = N_fibonacci(a, b, l)
    target_ratio = (b - a) / l;
    F_prev = 1;  
    F_curr = 1;  
    N = 2;     
    
    while F_curr < target_ratio
        F_next = F_prev + F_curr;
        F_prev = F_curr;
        F_curr = F_next;
        N = N + 1;
    end
end

 
f1 = @(x) (x-2).^2 + x .* log(x+3);
f2 = @(x) exp(-2*x) + (x-2).^2;
f3 = @(x) exp(x) .* (x.^3-1) + (x-1) .* sin(x);
a = -1;
b = 3;
epsilon = 0.001;
ls = linspace(0.001, 0.01, 100);
computations_f1 = zeros(1, length(ls));
computations_f2 = zeros(1, length(ls));
computations_f3 = zeros(1, length(ls));

for i = 1:length(ls)
    l = ls(i);
    N = N_fibonacci(a, b, l);
    [k1, ~, ~, ~] = fibonacci_method(f1, a, b, epsilon, N);
    [k2, ~, ~, ~] = fibonacci_method(f2, a, b, epsilon, N);
    [k3, ~, ~, ~] = fibonacci_method(f3, a, b, epsilon, N);
    computations_f1(i) = k1;
    computations_f2(i) = k2;
    computations_f3(i) = k3;
end

figure;
subplot(3, 1, 1);
plot(ls, computations_f1, 'Color', 'g');
xlabel('l');
ylabel('Iterations');
title('f_1(x): Variation of Iterations with l');
grid on;

subplot(3, 1, 2);
plot(ls, computations_f2, 'Color', 'r');
xlabel('l');
ylabel('Iterations');
title('f_2(x): Variation of Iterations with l');
grid on;

subplot(3, 1, 3);
plot(ls, computations_f3, 'Color', 'b');
xlabel('l');
ylabel('Iterations');
title('f_3(x): Variation of Iterations with l');
grid on;

l_values = [0.1, 0.05, 0.01, 0.005]; 
epsilon = 0.001;
figure;
for i = 1:3
    if i == 1
        f = f1;
        title_text = 'f_1(x) για διάφορες τιμές του l';
    elseif i == 2
        f = f2;
        title_text = 'f_2(x) για διάφορες τιμές του l';
    else
        f = f3;
        title_text = 'f_3(x) για διάφορες τιμές του l';
    end

    subplot(3, 1, i);
    hold on;
    for j = 1:length(l_values)
        l = l_values(j); 
        N = N_fibonacci(a, b, l);

        [~, a_vals, b_vals, ~] = fibonacci_method(f, a, b, epsilon, N);

        plot(1:length(a_vals), a_vals, 'DisplayName', ['a_k, l = ' num2str(l)]);
        plot(1:length(b_vals), b_vals, 'DisplayName', ['b_k, l = ' num2str(l)]);
    end
    xlabel('k');
    ylabel('Τιμή');
    title(title_text);
    legend('Location', 'best');
    grid on;
    hold off;
end