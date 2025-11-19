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

%% Plotting Configuration
% Define aesthetic preferences
lw_signal = 1.5;       % Thicker lines for main signals
lw_bound = 1.0;        % Slightly thinner for bounds
font_size = 12;        % Readable font size

% Colors defined using MATLAB's standard "lines" palette 
% These specific shades work well on both Dark and Light backgrounds
col_truth = [0.4660 0.6740 0.1880]; % Standard Green (Visible on dark)
col_est   = [0.0000 0.4470 0.7410]; % Standard Blue
col_bound = [0.8500 0.3250 0.0980]; % Standard Orange/Red
col_err   = [0.4940 0.1840 0.5560]; % Standard Purple (Replaces Black for Dark Mode visibility)

%% Figure 1: State Estimate Summary
% 'Name' labels the window, but we leave 'Color' default to respect your theme
fig1 = figure('Name', 'State Estimates'); 
sgtitle('State Estimate Summary by State', 'FontSize', font_size+2, 'FontWeight', 'bold');

ax1 = gobjects(nx, 1); % Pre-allocate for axis linking

for i = 1:nx
    ax1(i) = subplot(nx, 1, i);
    hold on; grid on; box on;
    
    % Calculate sigma for this state
    sigma_vec = sqrt(squeeze(P_kf(i, i, :)));
    
    % Plot Truth
    p_truth = plot(time_truth, x_truth(:, i), ...
        'Color', col_truth, 'LineWidth', lw_signal, 'DisplayName', 'Truth');
    
    % Plot Estimate
    p_est = plot(time_kf, x_kf(:, i), ...
        'Color', col_est, 'LineWidth', lw_signal, 'DisplayName', 'Estimate');
    
    % Plot Bounds (Dashed)
    plot(time_kf, x_kf(:, i) + sigma_vec, '--', ...
        'Color', col_bound, 'LineWidth', lw_bound, 'DisplayName', 'Uncertainty (1\sigma)');
    plot(time_kf, x_kf(:, i) - sigma_vec, '--', ...
        'Color', col_bound, 'LineWidth', lw_bound, 'HandleVisibility', 'off');
    
    ylabel(['State x_', num2str(i)], 'FontSize', font_size);
    set(gca, 'FontSize', font_size);
    
    % Legend only on the first subplot
    if i == 1
        legend([p_truth, p_est], 'Location', 'best');
    end
    
    hold off;
end
xlabel('Time [s]', 'FontSize', font_size);
linkaxes(ax1, 'x'); % Zooming one axis zooms them all

%% Figure 2: Estimation Error Summary
fig2 = figure('Name', 'Estimation Errors');
sgtitle('Estimation Error Summary by State', 'FontSize', font_size+2, 'FontWeight', 'bold');

ax2 = gobjects(nx, 1);

for i = 1:nx
    ax2(i) = subplot(nx, 1, i);
    hold on; grid on; box on;
    
    % Calculate sigma and error
    sigma_vec = sqrt(squeeze(P_kf(i, i, :)));
    error_vec = x_kf(:, i) - truth_kf(:, i);
    
    % Plot Error (Purple instead of Black for Dark Mode visibility)
    p_err = plot(time_kf, error_vec, ...
        'Color', col_err, 'LineWidth', lw_signal, 'DisplayName', 'Error');
        
    % Plot Bounds
    p_bnd = plot(time_kf, sigma_vec, '--', ...
        'Color', col_bound, 'LineWidth', lw_bound, 'DisplayName', '\pm 1\sigma Bound');
    plot(time_kf, -sigma_vec, '--', ...
        'Color', col_bound, 'LineWidth', lw_bound, 'HandleVisibility', 'off');
        
    ylabel(['Error x_', num2str(i)], 'FontSize', font_size);
    set(gca, 'FontSize', font_size);
    
    if i == 1
        legend([p_err, p_bnd], 'Location', 'best');
    end
    
    hold off;
end
xlabel('Time [s]', 'FontSize', font_size);
linkaxes(ax2, 'x');