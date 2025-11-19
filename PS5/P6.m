clear variables;
close all;
clc;

% This script is just testing for the sim 

addpath('..\Functions\');
rng(45, 'twister');

Fk     = [  0.81671934103521,  0.08791146849849;...
          -3.47061412053765,  0.70624978972000];     % for all k
Gammak = [  0.00464254201630;...
           0.08791146849849];                        % for all k
Hk     = [  2.00000000000000,  0.30000000000000];     % for all k

Qk     =    4.00000000000000;                         % for all k
Rk     =    0.01000000000000;                         % for all k

xhat0   = [  0.20000000000000;...
           -2.50000000000000];
P0      = [  0.25000000000000,  0.08000000000000;...
            0.08000000000000,  0.50000000000000];

[ta, xa, za] = mcltisim(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, 500);

[t_kf, x_kf, P_kf, nis_kf] = kalman_filter(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, za);

[t_srif, x_srif, P_srif, nis_srif] = srif(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, za);

% sys = ss(Fk, [eye(size(Fk)), Gammak], Hk, 0, -1);
% [km, L, P_bar_ss, W_ss] = kalman(sys, Qk, Rk);

% P_bar_ss

nis_srif
nis_kf

% Extract the std from cov
x1_cov = squeeze(sqrt(P_kf(1, 1, :)));
x2_cov = squeeze(sqrt(P_kf(2, 2, :)));

%% Plotting
% x1
figure;
hold on;
plot(t_kf, x_kf(:, 1), 'b');
plot(t_kf, x_kf(:, 1)+x1_cov, 'r');
plot(t_kf, x_kf(:, 1)-x1_cov, 'r');
hold off;

% Test the new plotting function
plot_kf(t_kf, x_kf, P_kf, ta, xa, 3);
plot_kf(t_srif, x_srif, P_srif, ta, xa, 3);



