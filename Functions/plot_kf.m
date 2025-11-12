function [] = plot_kf(time_kf, x_kf, P_kf, time_truth, x_truth)
%PLOT_KF Summary of this function goes here
%   Detailed explanation goes here

% Meta stuff
nx = size(x_kf, 2);
% Make a truth to align with the KF steps
n = size(x_truth, 1);
idx = repelem(1:n, 2); 
truth_kf = x_truth(idx, :);
truth_kf = truth_kf(2:end, :);

% Plot truth and each element of the state with its standard deviation (should it be variance?)
figure;
for i=1:nx
    subplot(nx, 1, i);
    hold on;
    plot(time_truth, x_truth(:, i), 'g');
    plot(time_kf, x_kf(:, i), 'b');
    plot(time_kf, x_kf(:, i) + sqrt(squeeze(P_kf(i, i, :))), 'r');
    plot(time_kf, x_kf(:, i) - sqrt(squeeze(P_kf(i, i, :))), 'r');
    hold off;
end


% Plot estimation error
figure;
for i=1:nx
    subplot(nx, 1, i);
    hold on;
    plot(time_kf, x_kf(:, i) - truth_kf(:, i), 'k'); % Estimation error
    plot(time_kf, sqrt(squeeze(P_kf(i, i, :))), 'r');
    plot(time_kf, -sqrt(squeeze(P_kf(i, i, :))), 'r');
    hold off;
end