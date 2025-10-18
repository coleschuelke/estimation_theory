clear variables;
close all;
clc;

h = @(x) atan(x);
H = @(x) 1 / (1 + x^2);
J = @(z, x) norm(z - h(x))^2;
z = 0;
sigma_squared = 0.5;


function [x_star, P_xx_star, steps] = gauss_newton(J, H, h, R, z, x_g, threshold, adjustStep)
    steps = zeros([1 100]);
    dx = 1;

    % Perform normalization
    Ra = chol(R);
    Rait = inv(Ra.');

    step = 1;
    while norm(dx) > threshold
        Hx = Rait*H(x_g);
        hx = Rait*h(x_g);
        dx = (Hx.'*Hx)\Hx.'*(Rait*z-hx);

        % Calculate step size
        alpha = 1;
        if adjustStep == true
            current_cost = J(z, x_g);
            new_cost = J(z, x_g + alpha * dx);
            while current_cost < new_cost
                alpha = alpha / 2;
                new_cost = J(z, x_g + alpha * dx);
            end
        end
        steps(step) = x_g;
        x_g = x_g + alpha * dx; % Update the estimate
        step = step + 1;
    end
    % Pass out final results
    x_star = x_g;
    P_xx_star = inv(Hx.'*Hx);
    
end

[x_adjust, P_xx_adjust, steps_adjust] = gauss_newton(J, H, h, sigma_squared, z, 1.5, 0.0005, true);
x_adjust
P_xx_adjust
steps_adjust(1:4)

[x, P_xx, steps] = gauss_newton(J, H, h, sigma_squared, z, 1.5, 0.0005, false);
x
P_xx
steps(1:4)

% Plot J
% Generate a range of x values for plotting
x_values = linspace(-6, 6, 100);
J_values = arrayfun(@(x) J(z, x), x_values);
J_steps = arrayfun(@(x) J(z, x), steps(1:4));
J_steps_adjust = arrayfun(@(x) J(z, x), steps_adjust(1:4));

% Create the plot
figure;
plot(x_values, J_values);
xlabel('x');
ylabel('J(z, x)');
title('Cost Function J vs. x');
hold on;
plot(steps(1:4), J_steps, 'r-');
plot(steps_adjust(1:4), J_steps_adjust, 'g-');
grid on;