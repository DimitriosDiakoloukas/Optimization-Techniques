clc;
clear;

function [k, a_vals, b_vals] = dichotomy_derivative(df, a, b, N)
    k = 1;
    a_vals = a;  
    b_vals = b; 
    while k < N
        xk = (a + b) / 2;
        df_xk = df(xk);

        if df_xk == 0
            break;
        end
        
        if df_xk > 0
            b = xk; 
        else
            a = xk; 
        end

        a_vals = [a_vals, a];
        b_vals = [b_vals, b];
        k = k + 1;
    end
end

function N = calculate_N(a, b, l)
    ratio = (b-a) / l;
    N = ceil(log(ratio) / log(2));
end

f1 = @(x) (x-2).^2 + x .* log(x+3);
df1 = @(x) 2 * (x-2) + log(x+3) + x / (x+3);

f2 = @(x) exp(-2*x) + (x-2).^2;
df2 = @(x) -2 * exp(-2*x) + 2 * (x-2);

f3 = @(x) exp(x) .* (x.^3-1) + (x-1) .* sin(x);
df3 = @(x) exp(x) * (3*x.^2) + (x.^3-1) * exp(x) + sin(x) + (x-1) * cos(x);

a = -1;
b = 3;
ls = linspace(0.001, 0.01, 100);
computations_f1 = zeros(1, length(ls));
computations_f2 = zeros(1, length(ls));
computations_f3 = zeros(1, length(ls));

for i = 1:length(ls)
    l = ls(i);
    N = calculate_N(a, b, l);
    disp(N);
    [k1, ~, ~] = dichotomy_derivative(df1, a, b, N);
    [k2, ~, ~] = dichotomy_derivative(df2, a, b, N);
    [k3, ~, ~] = dichotomy_derivative(df3, a, b, N);
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
figure;
for i = 1:3
    if i == 1
        f = f1;
        df = df1;
        title_text = 'f_1(x) για διάφορες τιμές του l';
    elseif i == 2
        f = f2;
        df = df2;
        title_text = 'f_2(x) για διάφορες τιμές του l';
    else
        f = f3;
        df = df3;
        title_text = 'f_3(x) για διάφορες τιμές του l';
    end

    subplot(3, 1, i);
    hold on;
    for j = 1:length(l_values)
        l = l_values(j); 
        N = calculate_N(a, b, l);
        [~, a_vals, b_vals] = dichotomy_derivative(df, a, b, N);

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

