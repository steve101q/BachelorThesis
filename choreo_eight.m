function dydt = three_body_rhs(~, y, G, m)
    r1 = y(1:2); r2 = y(3:4); r3 = y(5:6);
    v1 = y(7:8); v2 = y(9:10); v3 = y(11:12);

    a1 = G * m * ((r2 - r1)/norm(r2 - r1)^3 + (r3 - r1)/norm(r3 - r1)^3);
    a2 = G * m * ((r1 - r2)/norm(r1 - r2)^3 + (r3 - r2)/norm(r3 - r2)^3);
    a3 = G * m * ((r1 - r3)/norm(r1 - r3)^3 + (r2 - r3)/norm(r2 - r3)^3);

    dydt = [v1; v2; v3; a1; a2; a3];
end

% Constants
G = 1;             
m = 1; % same mass

% Initial conditions, given by Simó-Moore-Chenciner" solution
r1 = [-0.97000436,  0.24308753];
r2 = [ 0.97000436, -0.24308753];
r3 = [0, 0];

v1 = [0.4662036850, 0.4323657300];
v2 = [0.4662036850, 0.4323657300];
v3 = [-0.93240737, -0.86473146];

y0 = [r1, r2, r3, v1, v2, v3];

T = 6.3259;
tspan = linspace(0, T, 1000);

% Solve ODE
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[t, Y] = ode45(@(t, y) three_body_rhs(t, y, G, m), tspan, y0, opts);

r1 = Y(:, 1:2);
r2 = Y(:, 3:4);
r3 = Y(:, 5:6);

figure;
hold on;
axis equal;
axis([-1.5 1.5 -0.5 0.5]);
grid on;

% Plot orbits
plot(r1(:,1), r1(:,2), 'r--', 'LineWidth', 0.5);
plot(r2(:,1), r2(:,2), 'g--', 'LineWidth', 0.5);
plot(r3(:,1), r3(:,2), 'b--', 'LineWidth', 0.5);

% Animated masses
h1 = plot(r1(1,1), r1(1,2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
h2 = plot(r2(1,1), r2(1,2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
h3 = plot(r3(1,1), r3(1,2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

% Add text labels
label1 = text(r1(1,1)+0.05, r1(1,2)+0.05, 'm1', 'Color', 'r', 'FontWeight', 'bold');
label2 = text(r2(1,1)+0.05, r2(1,2)+0.05, 'm2', 'Color', 'g', 'FontWeight', 'bold');
label3 = text(r3(1,1)+0.05, r3(1,2)+0.05, 'm3', 'Color', 'b', 'FontWeight', 'bold');

% Animation 
while true
    for i = 1:length(t)
        % Masses positions update
        set(h1, 'XData', r1(i,1), 'YData', r1(i,2));
        set(h2, 'XData', r2(i,1), 'YData', r2(i,2));
        set(h3, 'XData', r3(i,1), 'YData', r3(i,2));

        % Label positions update
        set(label1, 'Position', r1(i,:) + [0.05 0.05]);
        set(label2, 'Position', r2(i,:) + [0.05 0.05]);
        set(label3, 'Position', r3(i,:) + [0.05 0.05]);

        pause(0.01);
     end
end
