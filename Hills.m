clc; clear; close all;
m1 = 5.972*10^24;
m2 = 7.34767309*10^22;

mu = m2/(m1+m2); % Mass 
%mu = 0.1; % basic

% Positions of masses
z1 = [-mu; 0];    % Earth
z2 = [1-mu; 0];   % Moon

% Grid
[X, Y] = meshgrid(linspace(-1.8, 1.8, 2500), linspace(-1.8, 1.8, 2500));

% Potential function
Phi = @(x,y) 0.5*(x.^2 + y.^2) + ...               
            (1-mu)./sqrt((x-z1(1)).^2 + (y-z1(2)).^2) + ... 
            mu./sqrt((x-z2(1)).^2 + (y-z2(2)).^2) + ...   
            0.5*mu*(1-mu);                      

% Potential over grid
Z = Phi(X, Y);

% Jacobi constants to plot (C values)
C_values = [3.1, 3.01, 3.5, 4.0]; 

% Figure
figure;
hold on;
axis equal;

% Primaries
plot(z1(1), z1(2), 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r'); 
plot(z2(1), z2(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', '#6433ff'); 

%Libration points
options = optimset('TolX',1e-12);

% Euler
L1 = fzero(@(x) x - (1-mu)*(x+mu)/abs(x+mu)^3 - mu*(x-1+mu)/abs(x-1+mu)^3, [-mu+0.01 1-mu-0.01], options);
L2 = fzero(@(x) x - (1-mu)*(x+mu)/abs(x+mu)^3 - mu*(x-1+mu)/abs(x-1+mu)^3, [1-mu+0.01 2], options);
L3 = fzero(@(x) x - (1-mu)*(x+mu)/abs(x+mu)^3 - mu*(x-1+mu)/abs(x-1+mu)^3, [-2 -mu-0.01], options);

% Lagrange
L4x = 0.5 - mu;
L4y = sqrt(3)/2;
L5x = L4x;
L5y = -L4y;

% Plot Libration points 
plot(L1, 0, 'm*', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'm'); 
plot(L2, 0, 'c*', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'c'); 
plot(L3, 0, 'g*', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'g'); 
plot(L4x, L4y, '*', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'b');
plot(L5x, L5y, 'k*', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'k');

% Colors
colors = lines(length(C_values));

% Hill's regions for different C values
for i = 1:length(C_values)
    C = C_values(i);
    contour(X, Y, 2*Z, [C C], 'LineWidth', 2, 'LineColor', colors(i,:), ...
           'DisplayName', sprintf('C = %.2f', C));
end

legend('m1', 'm2', 'L1', 'L2', 'L3', 'L4', 'L5', 'Location', 'best');
grid on;

hold off;