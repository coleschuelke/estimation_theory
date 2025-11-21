function [] = plot_kf(time_kf, x_kf, P_kf, time_truth, x_truth, num_sigma)
%PLOT_KF Plots KF results with shaded bounds and discrete alignment
%   
%   Usage:
%       plot_kf(t_kf, x_kf, P_kf)
%       plot_kf(t_kf, x_kf, P_kf, t_true, x_true)
%       plot_kf(t_kf, x_kf, P_kf, [], [], 3) 
%       plot_kf(t_kf, x_kf, P_kf, t_true, x_true, 3)

    %% 1. Input Parsing (Positional Implementation)
    p = inputParser;
    addRequired(p, 'time_kf');
    addRequired(p, 'x_kf');
    addRequired(p, 'P_kf');
    addOptional(p, 'time_truth', []);
    addOptional(p, 'x_truth', []);
    addOptional(p, 'num_sigma', 1, @isnumeric); % Changed from addParameter to addOptional
    
    parse(p, time_kf, x_kf, P_kf, time_truth, x_truth, num_sigma);
    
    % Extract results
    num_sigma_val = p.Results.num_sigma;
    has_truth = ~isempty(time_truth) && ~isempty(x_truth);
    nx = size(x_kf, 2);

    % Color Palette
    colors.truth = [0.4660 0.6740 0.1880]; % Green
    colors.est   = [0.0000 0.4470 0.7410]; % Blue
    colors.bound = [0.0000 0.4470 0.7410]; % Blue (used for fill)
    colors.err   = [0.4940 0.1840 0.5560]; % Purple
    alpha_val    = 0.15;                   % Transparency for bounds

    %% 2. Data Alignment (Discrete Expansion)
    x_truth_aligned = [];
    x_kf_plot = x_kf;
    t_plot = time_kf;
    P_plot = P_kf;

    if has_truth
        n_truth = size(x_truth, 1);
        
        % Duplicate and shift logic (Predict -> Update step alignment)
        idx = repelem(1:n_truth, 2);
        expanded_truth = x_truth(idx, :);
        
        if size(expanded_truth, 1) > 1
            expanded_truth = expanded_truth(2:end, :);
        end
        
        % Safe Truncation to match lengths exactly
        min_len = min(size(x_kf, 1), size(expanded_truth, 1));
        
        x_truth_aligned = expanded_truth(1:min_len, :);
        x_kf_plot       = x_kf(1:min_len, :);
        t_plot          = time_kf(1:min_len);
        P_plot          = P_kf(:, :, 1:min_len); 
    end

    %% 3. Figure 1: State Estimates
    fig1 = figure('Name', 'State Estimates');
    tlo1 = tiledlayout(nx, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    title(tlo1, 'State Estimate Summary', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(tlo1, 'Time [s]', 'FontSize', 12);

    ax_handles = gobjects(nx, 1);

    for i = 1:nx
        ax_handles(i) = nexttile;
        hold on; grid on; box on;
        
        % -- Prepare Data --
        sigma = sqrt(squeeze(P_plot(i,i,:))) * num_sigma_val;
        if size(sigma, 2) > size(sigma, 1), sigma = sigma'; end
        est_vec = x_kf_plot(:, i);
        
        % -- Plot Shaded Bounds --
        fill_t = [t_plot; flipud(t_plot)];
        fill_y = [est_vec + sigma; flipud(est_vec - sigma)];
        
        fill(fill_t, fill_y, colors.bound, ...
            'FaceAlpha', alpha_val, 'EdgeColor', 'none', ...
            'DisplayName', sprintf('%d\\sigma Bound', num_sigma_val));

        % -- Plot Truth --
        if has_truth
            plot(t_plot, x_truth_aligned(:, i), '-', ...
                'Color', colors.truth, 'LineWidth', 1.5, 'DisplayName', 'Truth');
        end

        % -- Plot Estimate --
        plot(t_plot, est_vec, ...
            'Color', colors.est, 'LineWidth', 1.5, 'DisplayName', 'Estimate');

        ylabel(['State x_', num2str(i)], 'FontWeight', 'bold');
        
        if i == 1
            legend('Location', 'best');
        end
    end
    linkaxes(ax_handles, 'x');

    %% 4. Figure 2: Estimation Errors (Only if Truth exists)
    if has_truth
        fig2 = figure('Name', 'Estimation Errors');
        tlo2 = tiledlayout(nx, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        title(tlo2, 'Estimation Error Summary', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel(tlo2, 'Time [s]', 'FontSize', 12);
        
        ax_err = gobjects(nx, 1);

        for i = 1:nx
            ax_err(i) = nexttile;
            hold on; grid on; box on;
            
            % -- Prepare Data --
            sigma = sqrt(squeeze(P_plot(i,i,:))) * num_sigma_val;
            if size(sigma, 2) > size(sigma, 1), sigma = sigma'; end
            error_vec = x_kf_plot(:, i) - x_truth_aligned(:, i);
            
            % -- Plot Shaded Bounds --
            fill_t = [t_plot; flipud(t_plot)];
            fill_y = [sigma; flipud(-sigma)];
            
            fill(fill_t, fill_y, colors.bound, ...
                'FaceAlpha', alpha_val, 'EdgeColor', 'none', ...
                'DisplayName', sprintf('\\pm%d\\sigma Bound', num_sigma_val));
            
            % -- Plot Error --
            plot(t_plot, error_vec, ...
                'Color', colors.err, 'LineWidth', 1.5, 'DisplayName', 'Error');
            
            ylabel(['Error x_', num2str(i)], 'FontWeight', 'bold');
            
            if i == 1
                legend('Location', 'best');
            end
        end
        linkaxes(ax_err, 'x');
    end
end