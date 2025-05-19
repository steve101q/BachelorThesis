function dydt = CR3BP_equations(~, y, mu)
    % Distances to primary bodies (2D only)
    r1 = sqrt((y(1) + mu)^2 + y(2)^2);
    r2 = sqrt((y(1) - 1 + mu)^2 + y(2)^2);

    % CR3BP eqts
    dydt = zeros(4,1);
    dydt(1) = y(3); 
    dydt(2) = y(4);  
    dydt(3) = 2*y(4) + y(1) - (1-mu)*(y(1)+mu)/r1^3 - mu*(y(1)-1+mu)/r2^3;
    dydt(4) = -2*y(3) + y(2) - (1-mu)*y(2)/r1^3 - mu*y(2)/r2^3; 
end

function simulate_Lagrange_point(mu, tspan, title_str, color, sign_y)
    % Position of L4 or L5, diff sign only 
    L_point = [0.5 - mu; sign_y * sqrt(3)/2];  

    % Initial cond
    y0 = [L_point(1) + 0.01; L_point(2) + sign_y * 0.01; 0; 0];

    % Solve
    options = odeset('RelTol', 1e-6);
    [~, y] = ode45(@(t,y) CR3BP_equations(t, y, mu), tspan, y0, options);

    % Plot trajectory
    plot(y(:,1), y(:,2), [color '-'], 'LineWidth', 1.5); hold on;
    plot(L_point(1), L_point(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    
    title(title_str);
    grid on; axis equal;
    
end


mu_stable =  0.01;    
mu_unstable =  0.03853;  

% Time 
tspan_stable = [0 200];  
tspan_unstable = [0 150];  

figure('Position', [100, 100, 900, 700]);

subplot(2, 2, 1);
simulate_Lagrange_point(mu_stable, tspan_stable, 'μ = 0.001', 'b', 1);

subplot(2, 2, 2);
simulate_Lagrange_point(mu_unstable, tspan_unstable, 'μ = 0.04', 'r', 1);

subplot(2, 2, 3);
simulate_Lagrange_point(mu_stable, tspan_stable, 'μ = 0.001', 'b', -1);

subplot(2, 2, 4);
simulate_Lagrange_point(mu_unstable, tspan_unstable, 'μ = 0.04', 'r', -1);
