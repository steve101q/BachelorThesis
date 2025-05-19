clc; clear; close all;

% Mass
m1 = 5.972e24;
m2 = 7.34767309e22;
mu2 = m2 / (m1 + m2);

% mu
mu_list = [0.001, mu2, 0.03852, 0.08]; 
% mu1= sun-jupiter, mu2=earth-moon, mu3 = threshold, mu4 =arbitrary high

% Jacobi C
C_values = 3.1;

% Grid
[X, Y] = meshgrid(linspace(-1.8, 1.8, 2500), linspace(-1.8, 1.8, 2500));

% Colors
colors = lines(length(C_values));

figure;

for k = 1:length(mu_list)
    mu = mu_list(k);

    % Positions primariess
    z1 = [-mu; 0];
    z2 = [1 - mu; 0];

    % Potential function 
    Phi = @(x, y) 0.5*(x.^2 + y.^2) + ...
                  (1 - mu)./sqrt((x - z1(1)).^2 + (y - z1(2)).^2) + ...
                  mu./sqrt((x - z2(1)).^2 + (y - z2(2)).^2) + ...
                  0.5 * mu * (1 - mu);

    % Potential on grid
    Z = Phi(X, Y);

    % Euler 
    options = optimset('TolX',1e-12);
    f = @(x) x - (1 - mu)*(x + mu)./abs(x + mu).^3 - mu*(x - 1 + mu)./abs(x - 1 + mu).^3;
    L1 = fzero(f, [-mu + 0.01, 1 - mu - 0.01], options);
    L2 = fzero(f, [1 - mu + 0.01, 2], options);
    L3 = fzero(f, [-2, -mu - 0.01], options);

    % Lagrange
    L4x = 0.5 - mu;
    L4y = sqrt(3)/2;
    L5x = L4x;
    L5y = -L4y;

    % Plotting
    subplot(2, 2, k);
    hold on;
    axis equal;
    title(sprintf('\\mu = %.5f', mu), 'FontWeight', 'bold');

    % Plot primaries
    plot(z1(1), z1(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(z2(1), z2(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', '#6433ff');

    % Plot libration points
    plot(L1, 0, 'm*', 'MarkerSize', 10, 'LineWidth', 1.5);
    plot(L2, 0, 'c*', 'MarkerSize', 10, 'LineWidth', 1.5);
    plot(L3, 0, 'g*', 'MarkerSize', 10, 'LineWidth', 1.5);
    plot(L4x, L4y, 'b*', 'MarkerSize', 10, 'LineWidth', 1.5);
    plot(L5x, L5y, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5);

    % Plot zero-velocity contours
    for i = 1:length(C_values)
        C = C_values(i);
        contour(X, Y, 2*Z, [C C], 'LineWidth', 1.5, 'LineColor', colors(i,:));
    end

    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    grid on;
end

% legend
legend({'m1', 'm2', 'L1', 'L2', 'L3', 'L4', 'L5'}, ...
    'Position', [0.9 0.25 0.1 0.5], 'FontSize', 8);
