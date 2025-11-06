clear variables;
close all;
clc;

% Load problem data into the workspace
run("kf_example02b.m");

%% Running the filter

% Initialize output
x = zeros(2, length(thist)*2 + 1);
x(:, 1) = xhat0;
P = repmat(zeros(size(P0)), 1, length(thist)*2 + 1);
P(:, 1:2) = P0;

% The filter
% Inititialize priors
x_post = xhat0;
P_post = P0;
for k=1:length(zhist)
    
    % Prediction step
    x_prior = Fk*x_post;
    P_prior = Fk*P_post*Fk.' + Gammak*Qk*Gammak.';
    
    % Store the prediction 
    x(:, 2*k) = x_prior;
    P(:, 4*k-1:4*k) = P_prior;

    % Update step
    innov = zhist(k) - Hk*x_prior;
    innov_cov = Hk*P_prior*Hk.' + Rk;
    K = P_prior*Hk.'/inv(innov_cov);
    x_post = x_prior + K*innov;
    P_post = P_prior - K*innov_cov*K.';

    % Store the correction
    x(:, 2*k + 1) = x_post;
    P(:, 4*k+1:4*k+2) = P_post;
end

% Extract the std from cov
x1_cov = sqrt(P(1, 1:2:end));
x2_cov = sqrt(P(2, 2:2:end));

%% Plotting
% CONFORM TO PLOTTING GUIDELINES IF TURNING IN
% x1
figure;
hold on;
plot([0 repelem(thist.', 2)], x(1, :), 'b');
plot([0 repelem(thist.', 2)], x(1, :)+x1_cov, 'r');
plot([0 repelem(thist.', 2)], x(1, :)-x1_cov, 'r');

figure;
hold on;
plot([0 repelem(thist.', 2)], x(2, :), 'b');
plot([0 repelem(thist.', 2)], x(2, :)+x2_cov, 'r');
plot([0 repelem(thist.', 2)], x(2, :)-x2_cov, 'r');