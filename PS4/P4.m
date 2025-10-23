clear variables;
close all;
clc;

%% Problem
% Load data and initiate constants
data = load("radarmeasdata_missile.mat");
t = data.thist;
rho_a = data.rhoahist;
rho_b = data.rhobhist;
z = zeros([2*length(t), 1]);
z(1:2:end) = rho_a;
z(2:2:end) = rho_b;
Ri = [100, 0; 0 900];
R = kron(eye(length(t)), Ri);

x_g = [8.1547e4; 1.2842e3; 1.2856e3; 1.2842e3]; % Initial guess

% Set up GN params
threshold = 1e-5;
adjustStep = true;

% Set up function handles
Hprime = @(x) Hprime_func(x, t);
hprime = @(x) hprime_func(x, t);
J = @(z, x, h) J_func(z, x, h);


% Run GN
[x_sol, P_xx_sol] = gauss_newton(J, Hprime, hprime, R, z, x_g, threshold, adjustStep);

x_sol
P_xx_sol


%% Functions
function [H] = Hprime_func(x, t)
    % Constants
    d = 3e4; % m
    l = 4.4e5; % m
    g = 9.81; % m/s^2

    H = zeros([2*length(t), 4]); % Initialize the H matrix
    
    % Break out states
    x0 = x(1);
    y0 = x(2);
    vx0 = x(3);
    vy0 = x(4);

    for j=1:length(t)
        tj = t(j);
        % Compute elements for each time step
        h1x1 = (d-l+vx0*tj+x0) / sqrt((d-l+vx0*tj+x0)^2 + (y0+vy0*tj-0.5*g*tj^2)^2);
        h1x2 = (y0+vy0*tj-0.5*g*tj^2) / sqrt((d-l+vx0*tj+x0)^2 + (y0+vy0*tj-0.5*g*tj^2)^2);
        h1x3 = tj*h1x1;
        h1x4 = tj*h1x2;
    
        h2x1 = (x0+vx0*tj-l) / sqrt((x0+vx0*tj-l)^2 + (y0+vy0*tj-0.5*g*tj^2)^2);
        h2x2 = (y0+vy0*tj-0.5*g*tj^2) / sqrt((x0+vx0*tj-l)^2 + (y0+vy0*tj-0.5*g*tj^2)^2);
        h2x3 = tj*h2x1;
        h2x4 = tj*h2x2;

        % Insert computed values into H
        H((j-1)*2 + 1, :) = [h1x1, h1x2, h1x3, h1x4];
        H((j-1)*2 + 2, :) = [h2x1, h2x2, h2x3, h2x4];
    end
end

function [h] = hprime_func(x, t)
    % Constants
    d = 3e4; % m
    l = 4.4e5; % m
    g = 9.81; % m/s^2

    h = zeros([2*length(t), 1]); % Initialize the h vector
    
    % Break out states
    x0 = x(1);
    y0 = x(2);
    vx0 = x(3);
    vy0 = x(4);
    
    % Compute h for each
    for j=1:length(t)
        tj = t(j);
        h1 = sqrt((d-l+vx0*tj+x0)^2 + (y0+vy0*tj-0.5*g*tj^2)^2);
        h2 = sqrt((x0+vx0*tj-l)^2 + (y0+vy0*tj-0.5*g*tj^2)^2);

        h((j-1)*2 + 1) = h1;
        h((j-1)*2 + 2) = h2;
    end

end

function [x_star, P_xx_star] = gauss_newton(J, Hprime, hprime, Rprime, zprime, x_g, threshold, adjustStep)
    dx = 1;

    % Perform normalization
    Ra = chol(Rprime);
    Rait = inv(Ra.');

    % Define normalized quantities
    H = @(x) Rait*Hprime(x);
    h = @(x) Rait*hprime(x);
    z = Rait*zprime;

    % cost = J(z, x_g, h);

    while norm(dx) > threshold

        Hx = H(x_g);
        hx = h(x_g);
        dx = (Hx.'*Hx)\Hx.'*(z-hx);

        % Calculate step size
        alpha = 1;
        if adjustStep == true
            current_cost = J(z, x_g, h);
            new_cost = J(z, x_g + alpha * dx, h);
            while current_cost < new_cost
                alpha = alpha / 2;
                new_cost = J(z, x_g + alpha * dx, h);
            end
        end
        x_g = x_g + alpha * dx; % Update the estimate
        cost = J(z, x_g, h); % Calculate the cost at each step
    end
    % Pass out final results
    x_star = x_g;
    P_xx_star = inv(Hx.'*Hx);
    cost  = J(z, x_g, h);
    
end

function [cost] = J_func(z, x, h)
    cost = norm(z - h(x))^2;
end

