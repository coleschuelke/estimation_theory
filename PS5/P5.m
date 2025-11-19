clear variables;
close all;
clc;

% Load problem data into the workspace
run("kf_example02b.m");

%% Running the filter

[t_kf1, x_kf1, P_kf1, nis1] = kalman_filter(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, zhist);
[t_kf2, x_kf2, P_kf2, nis2] = kalman_filter(Fk, Gammak, Hk, Qk/100, Rk, xhat0, P0, zhist);
[t_kf3, x_kf3, P_kf3, nis3] = kalman_filter(Fk, Gammak, Hk, Qk/10000, Rk, xhat0, P0, zhist);

% Calculating the NIS bounds

% Calculating RMS between best estimate and inaccurate models

%% Plotting
plot_kf(t_kf1, x_kf1, P_kf1, [], [], 3);
plot_kf(t_kf2, x_kf2, P_kf2, [], [], 3);
plot_kf(t_kf3, x_kf3, P_kf3, [], [], 3);

