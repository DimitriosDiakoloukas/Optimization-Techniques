clc;
clear;

syms x y
f = (1/3)*x^2 + 3*y^2;

eval_f = matlabFunction(f, 'Vars', [x, y]);

[x_vals, y_vals] = meshgrid(-10:0.1:10, -10:0.1:10);

z_vals = eval_f(x_vals, y_vals);

figure;
surf(x_vals, y_vals, z_vals);
title('3D Surface Plot of f(x)');
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
shading interp; 
colorbar; 
grid on;
