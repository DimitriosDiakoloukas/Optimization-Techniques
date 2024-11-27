clc; 
clear;

x = linspace(-4, 4, 1000);
y = linspace(-4, 4, 1000);
[X, Y] = meshgrid(x, y);
func  = X.^5.*exp(-X.^2-Y.^2);

figure;
surf(X, Y, func);
shading interp;
colormap('viridis'); 
colorbar;
title('3D Αναπαράσταση της f(x,y) = x^5 e^{-x^2-y^2}');
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
