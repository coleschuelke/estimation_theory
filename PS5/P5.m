clear variables;
close all;
clc;

% Set the seed
rng(37, "twister");

% Load problem data into the workspace
run("kf_example02b.m");

%% Running the filter

[t_kf1, x_kf1, P_kf1, nis1] = kalman_filter(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, zhist);
[t_kf2, x_kf2, P_kf2, nis2] = kalman_filter(Fk, Gammak, Hk, Qk/100, Rk, xhat0, P0, zhist);
[t_kf3, x_kf3, P_kf3, nis3] = kalman_filter(Fk, Gammak, Hk, Qk/10000, Rk, xhat0, P0, zhist);

% Calculating the NIS bounds
nis1_sum = sum(nis1);
nis2_sum = sum(nis2);
nis3_sum = sum(nis3);

stat1 = nis1_sum / 50
stat2 = nis2_sum / 50 % This one wins
stat3 = nis3_sum / 50

% Calculating RMS between best estimate and inaccurate models
rms21_1 = sqrt(mean((x_kf1(end-40:end, 1) - x_kf2(end-40:end, 1)).^2))
rms21_2 = sqrt(mean((x_kf1(end-40:end, 2) - x_kf2(end-40:end, 2)).^2))

rms23_1 = sqrt(mean((x_kf3(end-40:end, 1) - x_kf2(end-40:end, 1)).^2))
rms23_2 = sqrt(mean((x_kf3(end-40:end, 2) - x_kf2(end-40:end, 2)).^2))


%% Plotting
plot_kf(t_kf1, x_kf1, P_kf1, [], [], 3);
plot_kf(t_kf2, x_kf2, P_kf2, [], [], 3);
plot_kf(t_kf3, x_kf3, P_kf3, [], [], 3);

