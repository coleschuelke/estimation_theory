clear variables;
close all;
clc;

% Load measurements
thist =[0; 0.1000; 0.2000; 0.3000; 0.4000; 0.5000;
0.6000; 0.7000; 0.8000; 0.9000; 1.0000];

zprime = [7.7969; 1.4177; -3.0970; -7.6810; -9.8749; -6.1828;
-0.8212; 4.5074; 8.2259; 9.5369; 6.2827];

% Set up the problem
off_diag = 0.5 * ones([1, length(zprime)-1]);
R = eye(size(zprime, 1)) + diag(off_diag, 1) + diag(off_diag, -1);
Ra = chol(R);
Rait = inv(Ra.');
function h = hprime(x, t)
    h = zeros(length(t), 1);
    for j=1:length(t)
        h(j) = x(1)*cos(x(2)*t(j) + x(3));
    end
end
function H = Hprime(x, t)    
    H = zeros(length(t), length(x));
    for j=1:length(t)
        H(j, :) = [cos(x(2)*t(j) + x(3)), -t(j)*x(1)*sin(x(2)*t(j) + x(3)), -x(1)*sin(x(2)*t(j) + x(3))];
    end
end

% Gauss Newton Method
x_g = [4 6 0].';
convergence_thresh = 0.0002;
dx = 1;
while dx > convergence_thresh
    % Calculate dx
    Hx = Rait*Hprime(x_g, thist);
    hx = Rait*hprime(x_g, thist);
    dx = (Hx.'*Hx)\Hx.'*(Rait*zprime-hx);
    % Update guess
    x_g = x_g + 0.5 * dx;
end

% Final Results
x_g

Hfinal = Rait*Hprime(x_g, thist);
Pxx = inv(Hfinal.'*Hfinal)

J = @(x1, x2, x3) norm(Rait*zprime - Rait*hprime([x1, x2, x3], thist))^2;

plot3Dfunction(J, x_g, 10, 100)


% Plotting function
function plot3Dfunction(fun, center_point, range_val, num_points)
% plot3Dfunction Visualizes a function of three variables using 2D surface plots.
%
% This function creates three separate 2D surface plots, each representing a
% slice of the 3D function through a specified center point. The plots
% correspond to the XY, XZ, and YZ planes. It also plots a red marker on
% each surface to indicate the exact center point.
%
% This version is robust and uses 'arrayfun' to handle non-vectorized
% function handles.
%
% SYNTAX:
%   plot3Dfunction(fun, center_point, range_val)
%   plot3Dfunction(fun, center_point, range_val, num_points)
%
% INPUTS:
%   fun          - Function handle of the form @(x, y, z).
%   center_point - A 1x3 vector [x0, y0, z0] specifying the center of the
%                  plotting region. The planes will pass through this point.
%   range_val    - The range to plot around the center point.
%                  If scalar: The range for x, y, and z.
%                  If 1x3 vector: The respective ranges for [x, y, z].
%   num_points   - (Optional) The number of points for the grid. Default is 60.

% --- 1. Input Validation and Default Values ---
if nargin < 3
    error('This function requires at least 3 inputs: a function handle, a center point, and a range.');
end

if nargin < 4
    num_points = 60;
end

% Unpack center point and range
x0 = center_point(1);
y0 = center_point(2);
z0 = center_point(3);

if isscalar(range_val)
    rx = range_val; ry = range_val; rz = range_val;
elseif numel(range_val) == 3
    rx = range_val(1); ry = range_val(2); rz = range_val(3);
else
    error('range_val must be a scalar or a 3-element vector.');
end

% --- 2. Create the Grid Vectors ---
x_vec = linspace(x0 - rx, x0 + rx, num_points);
y_vec = linspace(y0 - ry, y0 + ry, num_points);
z_vec = linspace(z0 - rz, z0 + rz, num_points);

% --- 3. Evaluate the Function on Each Plane ---
fprintf('Evaluating function on three %d x %d planes...\n', num_points, num_points);

% XY Plane (z = z0)
[X_xy, Y_xy] = meshgrid(x_vec, y_vec);
V_xy = arrayfun(fun, X_xy, Y_xy, z0 * ones(size(X_xy)));

% XZ Plane (y = y0)
[X_xz, Z_xz] = meshgrid(x_vec, z_vec);
V_xz = arrayfun(fun, X_xz, y0 * ones(size(X_xz)), Z_xz);

% YZ Plane (x = x0)
[Y_yz, Z_yz] = meshgrid(y_vec, z_vec);
V_yz = arrayfun(fun, x0 * ones(size(Y_yz)), Y_yz, Z_yz);

% ADDED: Calculate the function value at the exact center point
v0 = fun(x0, y0, z0);

fprintf('Evaluation complete.\n');

% --- 4. Generate the Plots ---
figure('Name', 'Function Surface Slices');

% Plot 1: XY Plane
subplot(2, 2, 1);
surf(X_xy, Y_xy, V_xy);
hold on;
% ADDED: Plot the center point as a red dot
plot3(x0, y0, v0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Center Point');
hold off;
shading interp;
colorbar;
axis tight;
xlabel('X-axis', 'FontSize', 12);
ylabel('Y-axis', 'FontSize', 12);
zlabel('Function Value', 'FontSize', 12);
title(sprintf('XY Plane (z = %.2f)', z0), 'FontSize', 14);
view(3);

% Plot 2: XZ Plane
subplot(2, 2, 2);
surf(X_xz, Z_xz, V_xz);
hold on;
% ADDED: Plot the center point as a red dot
plot3(x0, z0, v0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Center Point');
hold off;
shading interp;
colorbar;
axis tight;
xlabel('X-axis', 'FontSize', 12);
ylabel('Z-axis', 'FontSize', 12);
zlabel('Function Value', 'FontSize', 12);
title(sprintf('XZ Plane (y = %.2f)', y0), 'FontSize', 14);
view(3);

% Plot 3: YZ Plane
subplot(2, 2, 3);
surf(Y_yz, Z_yz, V_yz);
hold on;
% ADDED: Plot the center point as a red dot
plot3(y0, z0, v0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Center Point');
hold off;
shading interp;
colorbar;
axis tight;
xlabel('Y-axis', 'FontSize', 12);
ylabel('Z-axis', 'FontSize', 12);
zlabel('Function Value', 'FontSize', 12);
title(sprintf('YZ Plane (x = %.2f)', x0), 'FontSize', 14);
view(3);

% Add a main title for the entire figure
sgtitle(sprintf('Function surfaces around [%.2f, %.2f, %.2f]', x0, y0, z0), 'FontSize', 16, 'FontWeight', 'bold');

disp('Plots generated successfully.');

end