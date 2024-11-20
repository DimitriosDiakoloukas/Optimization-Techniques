clc;
clear;

function [k, a_vals, b_vals, min_val] = dichotomy_method(f, a, b, l, epsilon)
    a_vals = []; 
    b_vals = [];  
    k = 0;
 
    while (abs(b-a) >= l)
        k = k + 1;
        x1 = (a+b) / 2 - epsilon;
        x2 = (a+b) / 2 + epsilon;
       
        if f(x1) < f(x2)
            b = x2;
        else
            a = x1;
        end
        a_vals = [a_vals, a];
        b_vals = [b_vals, b];
    end
    min_val = (a_vals+b_vals) / 2;
end

f1 = @(x) (x-2).^2 + x .* log(x+3);
f2 = @(x) exp(-2*x) + (x-2).^2;
f3 = @(x) exp(x) .* (x.^3-1) + (x-1) .* sin(x);
a = -1;
b = 3;  
l = 0.01; 
epsilons = linspace(0.0001, 0.007, 100); 

% Uncomment for testing while keeping commenting epsilons as linspace
% epsilons = 0.0001
% [~, ~, ~, min_val1] = dichotomy_method(f1, -1, 3, 0.01, epsilons);
% [~, ~, ~, min_val2] = dichotomy_method(f2, -1, 3, 0.01, epsilons);
% [~, ~, ~, min_val3] = dichotomy_method(f3, -1, 3, 0.01, epsilons);
% disp(min_val1);
% disp(min_val2);
% disp(min_val3);

computations_f1 = zeros(1, length(epsilons));
computations_f2 = zeros(1, length(epsilons));
computations_f3 = zeros(1, length(epsilons));

for i = 1:length(epsilons)
    epsilon = epsilons(i);
    if 2 * epsilon < l
        [k1, ~, ~, ~] = dichotomy_method(f1, a, b, l, epsilon);
        [k2, ~, ~, ~] = dichotomy_method(f2, a, b, l, epsilon);
        [k3, ~, ~, ~] = dichotomy_method(f3, a, b, l, epsilon);
    
        computations_f1(i) = k1;
        computations_f2(i) = k2;
        computations_f3(i) = k3;     
    end
end


figure;

subplot(3, 1, 1);
plot(epsilons, computations_f1, 'Color', 'g', 'DisplayName', 'f_1(x)');
xlabel('\epsilon');
ylabel('Υπολογισμοί');
title('f_1(x): Μεταβολή Αριθμού Υπολογισμών ως προς \epsilon');
legend('Location', 'best');
grid on;

subplot(3, 1, 2);
plot(epsilons, computations_f2, 'Color', 'r','DisplayName', 'f_2(x)');
xlabel('\epsilon');
ylabel('Υπολογισμοί');
title('f_2(x): Μεταβολή Αριθμού Υπολογισμών ως προς \epsilon');
legend('Location', 'best');
grid on;

subplot(3, 1, 3);
plot(epsilons, computations_f3, 'Color', 'b', 'DisplayName', 'f_3(x)');
xlabel('\epsilon');
ylabel('Υπολογισμοί');
title('f_3(x): Μεταβολή Αριθμού Υπολογισμών ως προς \epsilon');
legend('Location', 'best');
grid on;

epsilon = 0.001;
ls = linspace(0.001, 0.01, 100);  
computations_f1_l = zeros(1, length(ls));
computations_f2_l = zeros(1, length(ls));
computations_f3_l = zeros(1, length(ls));

for i = 1:length(ls)
    l = ls(i);
    if 2 * epsilon < l
        computations_f1_l(i) = dichotomy_method(f1, a, b, l, epsilon);
        computations_f2_l(i) = dichotomy_method(f2, a, b, l, epsilon);
        computations_f3_l(i) = dichotomy_method(f3, a, b, l, epsilon);
    end
end

figure;
subplot(3, 1, 1);
plot(ls, computations_f1_l, 'Color', 'g');
xlabel('l');
ylabel('Υπολογισμοί');
title('f_1(x): Μεταβολή Αριθμού Υπολογισμών ως προς l');
grid on;

subplot(3, 1, 2);
plot(ls, computations_f2_l, 'Color', 'r');
xlabel('l');
ylabel('Υπολογισμοί');
title('f_2(x): Μεταβολή Αριθμού Υπολογισμών ως προς l');
grid on;

subplot(3, 1, 3);
plot(ls, computations_f3_l, 'Color', 'b');
xlabel('l');
ylabel('Υπολογισμοί');
title('f_3(x): Μεταβολή Αριθμού Υπολογισμών ως προς l');
grid on;


epsilon = 0.001;
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

        [~, a_vals, b_vals, ~] = dichotomy_method(f, a, b, l, epsilon);

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
