clear variables;
close all;
clc;
time = 1:51;
steps = repelem(time, 2);
steps = steps(2:end);

% Load problem matrices
run('kf_example02b.m');

% Generate truth
[truth, meas] = mcltisim(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, 50);

% Run Kalman Filter
[x_hat, P, nis] = kalman_filter(Fk, Gammak, Hk, Qk/10, Rk, xhat0, P0, meas);

% Statistical evaluation
eps_bar = nis/50 % H0 says this is distributed as chi-squared with 

% Final covariance
P(:, :, 101)

%% Plotting
figure;
hold on;
plot(time, truth(:, 1), 'g');
plot(steps, x_hat(:, 1), 'magenta');
plot(steps, x_hat(:, 1) + squeeze(P(1, 1, :)), 'b');
plot(steps, x_hat(:, 1) - squeeze(P(1, 1, :)), 'b');

figure;
hold on;
plot(time, truth(:, 2), 'g');
plot(steps, x_hat(:, 2), 'magenta');
plot(steps, x_hat(:, 2) + squeeze(P(2, 2, :)), 'b');
plot(steps, x_hat(:, 2) - squeeze(P(2, 2, :)), 'b');

figure;
hold on;
plot(steps, squeeze(P(1, 1, :)), 'y');
plot(steps, squeeze(P(2, 2, :)), 'b');