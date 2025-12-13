clear variables;
close all;
clc;

% Import Problem Files
addpath('..\Functions\');
addpath('codeAndDataForStudents\');
load('problem3data.mat');
load('problem3truth.mat');


% Set up problem matrices
dt = 0.1; % s
nsteps = length(lidar);
Q = diag([0.25, 0.25, 3, 40*pi/180, 0.1, 0.1].^2)/dt;
R = diag([2*pi/180, 2*pi/180, 0.1].^2);
nx = size(Q, 1);
nz = size(R, 1);

% Initialize outputs
x_ukf = zeros(nsteps+1, nx);
P_ukf = zeros(nx, nx, nsteps+1);

% Set up filter inputs
dyn_handle = @(t, x) dyn_car(x);
h_handle = @(x, w) h_car([x; w]);
xhat = [90, 4.25, 13, pi, 5, 2].';
P = diag([2, 5, 1, pi/4, 4, 2].^2);
alpha = 1e-3;
beta = 2;
kappa = 0;

% Run the filter
for i=1:nsteps
    raw_meas = lidar(i).z;
    processed_meas = [min(raw_meas(:, 1)), max(raw_meas(:,1)), min(raw_meas(:, 2))];

    % Run the UKF step
    [xhat, P] = ukf_step(dyn_handle, h_handle, Q, R, alpha, beta, kappa, processed_meas, xhat, P);
    
    % Save the estimates
    x_ukf(i+1, :) = xhat.';
    P_ukf(:, :, i+1) = P;

end

%% Plotting
