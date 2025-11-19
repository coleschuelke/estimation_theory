function [] = plot_kf(time_kf, x_kf, P_kf, time_truth, x_truth, num_sigma)
%PLOT_KF Plots State Estimates and Errors (if truth is provided)
%   Inputs:
%       time_kf, x_kf, P_kf: Filter output (Required)
%       time_truth, x_truth: Truth data (Optional, pass [] to skip)
%       num_sigma: Number of standard deviations for bounds (Optional, default 1)

%% 1. Input Handling & Pre-processing
if nargin < 6 || isempty(num_sigma)
    num_sigma = 1; % Default to 1 sigma if not provided
end

% Determine if truth data is available
has_truth = (nargin >= 5) && ~isempty(time_truth) && ~isempty(x_truth);

nx = size(x_kf, 2);

% Align truth only if it exists
truth_kf = [];
if has_truth
    n_truth = size(x_truth, 1);
    % Alignment logic from original request
    idx = repelem(1:n_truth, 2); 
    truth_aligned = x_truth(idx, :);
    truth_kf = truth_aligned(2:end, :);
    
    % Safety check to ensure alignment didn't break dimensions
    if size(truth_kf, 1) ~= size(x_kf, 1)
        warning('Truth and KF dimensions mismatch after alignment. Error plotting might fail.');
    end
end

%% 2. Plotting Configuration
% Define aesthetic preferences
lw_signal = 1.5;       
lw_bound = 1.0;        
font_size = 12;        

% Colors (Dark/Light mode compatible)
col_truth = [0.4660 0.6740 0.1880]; % Green
col_est   = [0.0000 0.4470 0.7410]; % Blue
col_bound = [0.8500 0.3250 0.0980]; % Orange/Red
col_err   = [0.4940 0.1840 0.5560]; % Purple

%% 3. Figure 1: State Estimate Summary
fig1 = figure('Name', 'State Estimates'); 
sgtitle('State Estimate Summary by State', 'FontSize', font_size+2, 'FontWeight', 'bold');

ax1 = gobjects(nx, 1); 

for i = 1:nx
    ax1(i) = subplot(nx, 1, i);
    hold on; grid on; box on;
    
    % Calculate sigma bounds scaled by num_sigma
    std_dev = sqrt(squeeze(P_kf(i, i, :)));
    sigma_vec = num_sigma * std_dev;
    
    % Plot Truth (If available)
    if has_truth
        plot(time_truth, x_truth(:, i), ...
            'Color', col_truth, 'LineWidth', lw_signal, 'DisplayName', 'Truth');
    end
    
    % Plot Estimate
    plot(time_kf, x_kf(:, i), ...
        'Color', col_est, 'LineWidth', lw_signal, 'DisplayName', 'Estimate');
    
    % Plot Bounds (Dashed)
    bound_label = sprintf('Uncertainty (%d\\sigma)', num_sigma);
    plot(time_kf, x_kf(:, i) + sigma_vec, '--', ...
        'Color', col_bound, 'LineWidth', lw_bound, 'DisplayName', bound_label);
    plot(time_kf, x_kf(:, i) - sigma_vec, '--', ...
        'Color', col_bound, 'LineWidth', lw_bound, 'HandleVisibility', 'off');
    
    ylabel(['State x_', num2str(i)], 'FontSize', font_size);
    set(gca, 'FontSize', font_size);
    
    % Legend only on the first subplot
    if i == 1
        legend('Location', 'best');
    end
    
    hold off;
end
xlabel('Time [s]', 'FontSize', font_size);
linkaxes(ax1, 'x');

%% 4. Figure 2: Estimation Error Summary (Only if Truth exists)
if has_truth
    fig2 = figure('Name', 'Estimation Errors');
    sgtitle('Estimation Error Summary by State', 'FontSize', font_size+2, 'FontWeight', 'bold');
    
    ax2 = gobjects(nx, 1);
    
    for i = 1:nx
        ax2(i) = subplot(nx, 1, i);
        hold on; grid on; box on;
        
        % Calculate sigma and error
        std_dev = sqrt(squeeze(P_kf(i, i, :)));
        sigma_vec = num_sigma * std_dev;
        
        % Ensure truth_kf matches x_kf size for subtraction
        % (Using simple sizing check to prevent crash if alignment failed)
        n_pts = min(size(x_kf, 1), size(truth_kf, 1));
        error_vec = x_kf(1:n_pts, i) - truth_kf(1:n_pts, i);
        
        % Plot Error
        plot(time_kf(1:n_pts), error_vec, ...
            'Color', col_err, 'LineWidth', lw_signal, 'DisplayName', 'Error');
            
        % Plot Bounds
        bound_label = sprintf('\\pm %d\\sigma Bound', num_sigma);
        plot(time_kf, sigma_vec, '--', ...
            'Color', col_bound, 'LineWidth', lw_bound, 'DisplayName', bound_label);
        plot(time_kf, -sigma_vec, '--', ...
            'Color', col_bound, 'LineWidth', lw_bound, 'HandleVisibility', 'off');
            
        ylabel(['Error x_', num2str(i)], 'FontSize', font_size);
        set(gca, 'FontSize', font_size);
        
        if i == 1
            legend('Location', 'best');
        end
        
        hold off;
    end
    xlabel('Time [s]', 'FontSize', font_size);
    linkaxes(ax2, 'x');
end

end