clear all;
close all;
clc;

% Root finding parameters
convergence_thresh = 1e-6;
maxiter = 100;
vector = -10:10;

% Initialize solution vectors
sol = zeros(21, 21, 2);
step_sols = zeros(3, 6, 2);

% Create x_g grid
[C, R] = meshgrid(vector, vector);
grid = cat(3, R, C);

for row=vector
    for col=vector
        x_g = squeeze(grid(row+11, col+11, :));
        saved_sol = 0;
        save = false;
        if isequal(x_g, [4; -4]) || isequal(x_g, [6; 0]) || isequal(x_g, [-5; 5])
            step_sol = zeros(2, 6);
            save = true;
            saved_sol = saved_sol + 1;
        end
        
        step = inf;
        iter = 1;
        while step > convergence_thresh && iter < maxiter
            
            if save == true && iter <= 6
                step_sol(:, iter) = x_g;
            end

            x = x_g;
            f = [x(1) + x(2) + x(1)*x(2) + 5;
            x(1)^2 + 2*x(2) - x(2)^2 - 2];
            dfdx = [1+x(2), 1 + x(1);
            2*x(1), 2-2*x(2)];
        
            x_g = x_g - dfdx\f;
            step = norm(x_g - x);

            iter = iter + 1;
        end
        if iter == maxiter
            x_g = [inf, inf];
        end

        % Save the solution
        sol(row+11, col+11, :) = x_g;
        if save == true
            step_sols(saved_sol, :, :) = step_sol.'; % TODO
        end

    end
end

% Reshape to columns
col_sol = reshape(permute(sol, [3 1 2]), size(sol, 3), []).';
col_guess = reshape(permute(grid, [3 1 2]), size(grid, 3), []).';
% Round for finding uniqueness
col_sol = round(col_sol, 4);
[unique_sols, ~, labels] = unique(col_sol, 'rows');

num_unique_sols = size(unique_sols, 1);

% --- 3. Plotting by Group ---

figure;
hold on;

% Define the symbols and colors we'll cycle through
symbols = ['o', 's', '^', 'v', 'P', '*', 'X', 'D'];
colors = lines(num_unique_sols); % 'lines' is a good colormap

plot_handles = zeros(num_unique_sols, 1); % To store handles for the legend
legend_entries = cell(num_unique_sols, 1);

for k = 1:num_unique_sols
    % Get the symbol and color for this group
    sym = symbols(mod(k-1, length(symbols)) + 1);
    col = colors(k, :);
    
    % --- A. Plot all initial guesses for this group ---
    
    % Find the indices (rows) of all guesses that belong to group 'k'
    indices_for_this_group = (labels == k);
    
    % Get the [x, y] coordinates of these guesses
    guesses_to_plot = col_guess(indices_for_this_group, :);
    
    % Plot the cloud of initial guesses
    plot(guesses_to_plot(:, 1), guesses_to_plot(:, 2), ...
        'Marker', sym, ...
        'Color', col, ...
        'LineStyle', 'none', ...
        'HandleVisibility', 'off'); % Hide from legend
        
    % --- B. Plot the unique solution itself ---
    
    % Get the [u, v] coordinates of this group's solution
    sol = unique_sols(k, :);
    
    % Plot the solution point to be highly visible
    % (Larger, with a black edge and solid face color)
    plot_handles(k) = plot(sol(1), sol(2), ...
        'Marker', sym, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12, ...
        'LineWidth', 2);
        
    % Store the text for the legend
    legend_entries{k} = sprintf('Solution %d: (%.2f, %.2f)', k, sol(1), sol(2));
end

% --- Finalize Plot ---
hold off;
axis equal; % Ensures x and y axes have the same scale
grid on;
xlabel('X-coordinate');
ylabel('Y-coordinate');
title('Initial Guesses Grouped by Convergent Solution');
legend(plot_handles, legend_entries, 'Location', 'best');