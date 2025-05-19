clc; clear; close all;
m1 = 5.972*10^24;
m2 = 7.34767309*10^22;

mu = m2/(m1+m2); % Mass 
r = 1;     % distance between primaries 
%mu=0.1; %basic case

% Positions of primaries
x1 = -mu*r;    % Earth
x2 = (1-mu)*r; % Moon

% Distances from origin 
r1 = abs(x1);
r2 = abs(x2);

% Libration points
options = optimset('TolX',1e-12);

% Euler
L1 = fzero(@(x) x - (1-mu)*(x+mu)/abs(x+mu)^3 - mu*(x-1+mu)/abs(x-1+mu)^3, [x1+0.01 x2-0.01], options);
L2 = fzero(@(x) x - (1-mu)*(x+mu)/abs(x+mu)^3 - mu*(x-1+mu)/abs(x-1+mu)^3, [x2+0.01 x2+0.5], options);
L3 = fzero(@(x) x - (1-mu)*(x+mu)/abs(x+mu)^3 - mu*(x-1+mu)/abs(x-1+mu)^3, [-1.5 x1-0.01], options);

% Lagrange
L4x = (1-2*mu)*r/2;
L4y = sqrt(3)*r/2;
L5x = L4x;
L5y = -L4y;


figure;
hold on;
axis equal;
grid on;
xlim([-1.25 1.25]);
ylim([-1.25 1.25]);

% Circles rad is distance to primary
theta = linspace(0, 2*pi, 100);
plot(r1*cos(theta), r1*sin(theta), 'r-', 'LineWidth', 1); 
plot(r2*cos(theta), r2*sin(theta), 'b-', 'LineWidth', 1); 

% Plot Primaries
plot(x1, 0, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r'); % m1
plot(x2, 0, 'bo', 'MarkerSize', 10*(1-mu), 'MarkerFaceColor', 'b'); % m2

% Plot Lagrange points
plot(L1, 0, 'k*', 'MarkerSize', 10); % L1
plot(L2, 0, 'k*', 'MarkerSize', 10); % L2
plot(L3, 0, 'k*', 'MarkerSize', 10); % L3
plot(L4x, L4y, 'k*', 'MarkerSize', 10); % L4
plot(L5x, L5y, 'k*', 'MarkerSize', 10); % L5

% Plot triangle Lagrange pts
plot([x1, L4x, x2, x1], [0, L4y, 0, 0], 'm--', 'LineWidth', 1); % Triangle to L4
plot([x1, L5x, x2, x1], [0, L5y, 0, 0], 'm--', 'LineWidth', 1); % Triangle to L5

% Plot axes
plot([-1.5 1.5], [0 0], 'k-', 'LineWidth', 0.5); % x-axis
plot([0 0], [-1.5 1.5], 'k-', 'LineWidth', 0.5); % y-axis

% Labels
text(x1, 0.1, 'm_1', 'HorizontalAlignment', 'center');
text(x2, 0.1, 'm_2', 'HorizontalAlignment', 'center');
text(L1, 0.1, 'L_1', 'HorizontalAlignment', 'center');
text(L2, 0.1, 'L_2', 'HorizontalAlignment', 'center');
text(L3, 0.1, 'L_3', 'HorizontalAlignment', 'center');
text(L4x, L4y+0.1, 'L_4', 'HorizontalAlignment', 'center');
text(L5x, L5y-0.1, 'L_5', 'HorizontalAlignment', 'center');