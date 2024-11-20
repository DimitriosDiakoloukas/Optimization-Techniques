clc;
clear;

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

f1 = @(x) (x-2).^2 + x .* log(x+3);
f2 = @(x) exp(-2*x) + (x-2).^2;
f3 = @(x) exp(x) .* (x.^3-1) + (x-1) .* sin(x);
a = -1;
b = 3;

% --- Μεταβολή του l με τη Μέθοδο Χρυσής Τομής ---
ls = linspace(0.001, 0.01, 100);
computations_f1 = zeros(1, length(ls));
computations_f2 = zeros(1, length(ls));
computations_f3 = zeros(1, length(ls));

for i = 1:length(ls)
    l = ls(i);
    computations_f1(i) = golden_section_method(f1, a, b, l);
    computations_f2(i) = golden_section_method(f2, a, b, l);
    computations_f3(i) = golden_section_method(f3, a, b, l);
end

figure;
subplot(3, 1, 1);
plot(ls, computations_f1, 'Color', 'g');
xlabel('l');
ylabel('Υπολογισμοί');
title('f_1(x): Μεταβολή Αριθμού Υπολογισμών ως προς l');
grid on;

subplot(3, 1, 2);
plot(ls, computations_f2, 'Color', 'r');
xlabel('l');
ylabel('Υπολογισμοί');
title('f_2(x): Μεταβολή Αριθμού Υπολογισμών ως προς l');
grid on;

subplot(3, 1, 3);
plot(ls, computations_f3, 'Color', 'b');
xlabel('l');
ylabel('Υπολογισμοί');
title('f_3(x): Μεταβολή Αριθμού Υπολογισμών ως προς l');
grid on;

l_values = [0.1, 0.05, 0.01, 0.005]; 
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

        [~, a_vals, b_vals, ~] = golden_section_method(f, a, b, l);

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

