function [NEES, SNEES, RMSE] = filter_analysis(x_filter, P_filter, x_truth)
%FILTER_ANALYSIS Analyze a filter
%   Plot, calculate statistics, etc. Currently configured for a single run

nx = size(x_filter, 2);
num_meas = size(x_filter(1:2:end, :), 1);

x_est = x_filter(1:2:end, :); % Keep only the posteriors
P_est = P_filter(:, :, 1:2:end);

% Calculate NEES
x_tilde = x_est - x_truth;
x_tilde_temp = reshape(x_tilde.', nx, 1, num_meas);
x_tilde_mult = pagemldivide(P_est, x_tilde_temp);
NEES = sum(squeeze(pagemtimes(reshape(x_tilde.', 1, nx, num_meas), x_tilde_mult)))/num_meas;

% Calculate SNEES
SNEES = NEES/nx;

% Calculate RMSE
RMSE = sqrt(mean((x_est - x_truth).^2, 1));

% Calculate 99% bounds for true consistency
nees_ub = chi2inv(0.995, num_meas*nx)/num_meas;
nees_lb = chi2inv(0.005, num_meas*nx)/num_meas;

snees_ub = nees_ub / nx;
snees_lb = nees_lb / nx;

% Check against bounds
nees_cons = NEES < nees_ub && NEES > nees_lb;
snees_cons = SNEES < snees_ub && SNEES > snees_lb;

%% Visualization
figure('Name', 'Filter Performance Summary', 'Color', 'w', 'NumberTitle', 'off');

% --- Top Panel: RMSE per State ---
subplot(2, 2, [1 2]); 
b = bar(1:nx, RMSE, 'FaceColor', [0.2, 0.4, 0.6]);
grid on;
title('Root Mean Square Error (RMSE) per State');
xlabel('State Index');
ylabel('RMSE Magnitude');
xticks(1:nx);
% Add value labels on top of bars
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(round(b.YData, 3));
text(xtips, ytips, labels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% --- Bottom Left: NEES Consistency ---
subplot(2, 2, 3);
hold on;
% Plot Bounds (Capture handle of one for the legend)
h_nees_bound = yline(nees_ub, 'r--', 'LineWidth', 1.5);
yline(nees_lb, 'r--', 'LineWidth', 1.5);

% Plot Actual NEES
if nees_cons
    markerColor = 'g';
    statusText = 'PASSED';
else
    markerColor = 'r';
    statusText = 'FAILED';
end
h_nees_pt = plot(1, NEES, 'o', 'MarkerSize', 10, 'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'k');

% Styling
title(['NEES Consistency: ' statusText]);
ylabel('NEES Value');
% Legend: Only pass the handle for one bound line and the data point
legend([h_nees_bound, h_nees_pt], {'99% Bounds', 'Calculated NEES'}, 'Location', 'best');
set(gca, 'XTick', []); 
grid on;
hold off;

% --- Bottom Right: SNEES Consistency ---
subplot(2, 2, 4);
hold on;
% Plot Bounds
h_snees_bound = yline(snees_ub, 'b--', 'LineWidth', 1.5);
yline(snees_lb, 'b--', 'LineWidth', 1.5);

% Plot Actual SNEES
if snees_cons
    markerColor = 'g';
    statusText = 'PASSED';
else
    markerColor = 'r';
    statusText = 'FAILED';
end
h_snees_pt = plot(1, SNEES, 's', 'MarkerSize', 10, 'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'k');

% Styling
title(['SNEES Consistency: ' statusText]);
ylabel('SNEES Value');
% Legend: Only pass the handle for one bound line and the data point
legend([h_snees_bound, h_snees_pt], {'99% Bounds', 'Calculated SNEES'}, 'Location', 'best');
set(gca, 'XTick', []); 
grid on;
hold off;

end