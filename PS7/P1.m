clear variables;
close all;
clc;

addpath('..\Functions\')

prop_orbit_handle = @(t, x, u, v) prop_orbit(t, x, u, v);
mu = 3.986005e14;

rk = [780000 + 6378000; 0; 0];
vk = [0; sqrt(mu/rk(1))*cosd(30); sqrt(mu/rk(1))*sind(30)];

xk = [rk; vk];

second_denom = 10; % 0.1s

T = 2*pi*sqrt((rk(1)^3)/mu);

num_steps = 10000; 
dt  = T/num_steps; % s
t = 0:dt:T;

x = zeros(6, length(t));
F = zeros(6, 6, length(t));
Gamma = zeros(6, 3, length(t));
x(:, 1) = [rk; vk];

for i=1:length(t)-1
    [xkp1, Fk, Gammak] = discrete_prop(xk, zeros(3, 1), zeros(3, 1), t(i), dt, prop_orbit_handle);
    x(:, i+1) = xkp1;
    F(:, :, i) = Fk;
    Gamma(:, :, i) = Gammak;

    xk = xkp1;
end

total_error = norm(x(1:3, 1)-x(1:3, end))

%% Testing Linearization
xk = [20000, 0, 0, 0, 0, 4000].';
xbark = [20001, 0, 0, 0, 0, 3999].';
[xkp1, Fk, Gammak] = discrete_prop(xk, zeros(3, 1), zeros(3, 1), 0, 0.1, prop_orbit_handle);
[xbarkp1, ~, ~] = discrete_prop(xbark, zeros(3, 1), zeros(3, 1), 0, 0.1, prop_orbit_handle);
d1 = xkp1-xbarkp1;
d2 = Fk*(xk-xbark);
der_error = norm(d1-d2)

%% Plotting

figure;
plot(x(1, :), x(2, :));
axis equal;

figure;
hold on;
plot3(x(1, :), x(2, :), x(3, :))
plot3(x(1, 1), x(2, 1), x(3, 1), 'bo')
plot3(x(1, end), x(2, end), x(3, end), 'rx')
plot3(x(1, end-1), x(2, end-1), x(3, end-1), 'go')
axis equal;