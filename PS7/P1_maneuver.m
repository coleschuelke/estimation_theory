clear variables;
close all;
clc;

addpath('..\Functions\')

prop_orbit_handle = @(t, x, u, v) prop_orbit(t, x, u, v);
mu = 3.986005e14;

rk = [780000 + 6378000; 0; 0];
vk = [0; sqrt(mu/rk(1)); 0];

xk = [rk; vk];

T = round(2*pi*sqrt((rk(1)^3)/mu)*10); % 0.1s

x = zeros(6, 3*T+1);
x(:, 1) = [rk; vk];

u = zeros(3, 33*T);
u(:, 25001:30000) = repmat([.5; .1; 0], 1, 5000);
u(:, 95001:100000) = repmat([0; .7; .5], 1, 5000);

for i=1:3*T
    [xk, ~, ~] = discrete_prop(xk, u(:, i), zeros(3, 1), (i-1)/10, 1/10, prop_orbit_handle);
    x(:, i+1) = xk;
end

total_error = norm(x(1:3, 1)-x(1:3, end))

%% Plotting

figure;
plot(x(1, :), x(2, :));
axis equal;

figure;
hold on;
plot3(x(1, :), x(2, :), x(3, :))
plot3(x(1, 1), x(2, 1), x(3, 1), 'bo')
plot3(x(1, end), x(2, end), x(3, end), 'rx')
axis equal;