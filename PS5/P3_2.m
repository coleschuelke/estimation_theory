clear variables;
close all;
clc;

addpath('..\Functions\')

% Load problem data into the workspace
run("kf_example02a.m");

% The filter

[t_out, x, P_out] = kalman_filter(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, zhist);

P_final = P_out(:, :, end)
P_mid = P_out(:, :, ceil(end/2))

% Extract the std from cov
x1_cov = squeeze(sqrt(P_out(1, 1, :)));
x2_cov = squeeze(sqrt(P_out(2, 2, :)));

%% Plotting
% CONFORM TO PLOTTING GUIDELINES IF TURNING IN
% x1
figure;
hold on;
plot([0 repelem(thist.', 2)], x(:, 1), 'b');
plot([0 repelem(thist.', 2)], x(:, 1)+x1_cov, 'r');
plot([0 repelem(thist.', 2)], x(:, 1)-x1_cov, 'r');

% x2
figure;
hold on;
plot([0 repelem(thist.', 2)], x(:, 2), 'b');
plot([0 repelem(thist.', 2)], x(:, 2)+x2_cov, 'r');
plot([0 repelem(thist.', 2)], x(:, 2)-x2_cov, 'r');

% Looks a whole lot better than the last one
figure;
plot([0 repelem(thist.', 2)], x1_cov);