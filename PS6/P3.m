clear variables;
close all;
clc;

addpath(['..\Functions\']);

rng(43, "twister");

% Load the matrices into the workspace
run('kf_example03a.m');

% Kalman Filter
[t_kf, x_kf, P_kf, nis_kf] = kalman_filter(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, zhist);

% Square-root Information Filter
[t_srif, x_srif, P_srif, nis_srif] = srif(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, zhist);

% Square-root Information Smoother

[t_sris, x_sris, P_sris, ~, ~] = sris(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, zhist);


% Plotting
plot_kf(t_kf, x_kf, P_kf, [], [], 3);
plot_kf(t_srif, x_srif, P_srif, [], [], 3);
plot_kf(t_sris, x_sris, P_sris, [], [], 3);