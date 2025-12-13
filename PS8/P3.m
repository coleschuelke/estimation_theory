clear variables;
close all;
clc;

function [cartesian] = h_nl(r, theta)
    x = r*cos(theta);
    y = r*sin(theta);
    cartesian = [x; y];
end

function [H] = H(rbar, thetabar)
    H = [cos(thetabar), -rbar*sin(thetabar); sin(thetabar), rbar*cos(thetabar)];
end

function [cartesian] = H_lin(rbar, thetabar, r, theta)
    H(rbar, thetabar);
    cartesian = h_nl(rbar, thetabar) + H*[r-rbar; theta-thetabar];
end

% Case i
rbari = 76;
thetabari = -3*pi/180;
P_pi = diag([1^2, (pi/180)^2]);

cbari = h_nl(rbari, thetabari)
Pcci = H(rbari, thetabari)*P_pi*H(rbari, thetabari).'

% Case ii
rbarii = 76;
thetabarii = -3*pi/180;
P_pii = diag([1^2, (15*pi/180)^2]);

cbarii = h_nl(rbarii, thetabarii)
Pccii = H(rbarii, thetabarii)*P_pii*H(rbarii, thetabarii).'

%% Unscented Transform
function [xbar_t, Pxx_t] = ut(xbar, Pxx, nl_fun, alpha, beta, kappa, varargin)
    
    % varargin to take more points

    nx = length(xbar);
    npts = nx; % Doesn't include the mean point and includes symmetry
    Sx = chol(Pxx).';
    lambda = (alpha^2)*(nx+kappa) - nx;

    % Build points
    chi = zeros(nx, 2*npts+1);
    chi0 = xbar;
    chi(:, 1) = chi0;
    for i=1:npts
        chi(:, 2*i) = chi0 + sqrt(nx + lambda)*Sx(:,i);
        chi(:, 2*i+1) = chi0 - sqrt(nx + lambda)*Sx(:,i);
    end

    % Push points through function
    % chi_t = arrayfun(nl_fun, chi); % I know I found a way to do this the other day
    for m=1:2*npts+1
        chi_t(:, m) = nl_fun(chi(:, m));
    end

    w0 = lambda/(nx+lambda);
    wi = 1/(2*(nx+lambda));

    xbar_t = w0*chi_t(:, 1) + sum(wi*chi_t(:, 2:end), 2);
    Pxx_t = (lambda/(nx+lambda) + 1 - alpha^2 + beta)*(chi_t(:, 1) - xbar_t)*(chi_t(:, 1) - xbar_t).';
    for k=2:2*npts+1
        Pxx_t = Pxx_t + wi*(chi_t(:, k) - xbar_t)*(chi_t(:, k) - xbar_t).';
    end
end

nl_handle = @(x) h_nl(x(1), x(2));

[cbar_uti, Pcc_uti] = ut([rbari; thetabari], P_pi, nl_handle, 1e-3, 2, 0)
[cbar_utii, Pcc_utii] = ut([rbarii; thetabarii], P_pii, nl_handle, 1e-3, 2, 0)


% Validation with vectors
polar_points = zeros(2, 1e6);
for i=1:1e6
    polar_points(:, i) = draw_gaussian([76; -3*pi/180], diag([1^2, (15*pi/180)^2]));
end

cartesian_points = zeros(2, 1e6);
for i=1:1e6
    cartesian_points(:, i) = h_nl(polar_points(1, i), polar_points(2, i));
end

cartesian_mean = mean(cartesian_points, 2)
cartesian_cov = cov(cartesian_points.')

plot_cov_comparison([76; -3*pi/180], diag([1^2, (15*pi/180)^2]), cartesian_mean, cartesian_cov, cbarii, Pccii, cbar_utii, Pcc_utii);

function plot_cov_comparison(mu_prior, P_prior, mu_true, P_true, mu_lin, P_lin, mu_ut, P_ut)
% PLOT_COV_COMPARISON 
% Visualizes and compares 4 Gaussian distributions.

    % Create figure
    if isempty(get(groot,'CurrentFigure'))
        figure('Color', 'w', 'Name', 'Covariance Comparison');
    else
        set(gcf, 'Color', 'w');
    end
    
    hold on; grid on; axis equal;
    
    % --- 1. Prior (Grey) ---
    plot_gaussian(mu_prior, P_prior, [0.5 0.5 0.5], '--', 'o');
    
    % --- 2. True Posterior (Green) ---
    plot_gaussian(mu_true, P_true, [0 0.6 0], '-', 'p');
    
    % --- 3. Linearized/EKF (Red) ---
    plot_gaussian(mu_lin, P_lin, [0.8 0 0], '-.', 'x');
    
    % --- 4. Unscented/UKF (Blue) ---
    plot_gaussian(mu_ut, P_ut, [0 0 0.8], '-', '+');
    
    legend('Prior', 'True Posterior', 'Linearized (EKF)', 'Unscented (UKF)', 'Location', 'best');
    title('Comparison of Mean and Covariance Estimates');
    xlabel('State X'); ylabel('State Y');
    set(gca, 'FontSize', 12);
    hold off;
end

% ---------------------------------------------------------
% HELPER FUNCTION (ROBUST VERSION)
% ---------------------------------------------------------
function output = plot_gaussian(mu, P, col, style, marker)
    % --- 1. Handle Missing Arguments (Defaults) ---
    if nargin < 3, col = 'b'; end       % Default: Blue
    if nargin < 4, style = '-'; end     % Default: Solid Line
    if nargin < 5, marker = '.'; end    % Default: Dot

    % Settings
    n_sig = 2;    % 2-sigma
    lw    = 2;    % Line Width
    ms    = 10;   % Marker Size

    mu = mu(:);   % Ensure column vector
    
    % --- 2. Calculate Ellipse Geometry ---
    t = linspace(0, 2*pi, 100);
    unit_circle = [cos(t); sin(t)];
    
    [V, D] = eig(P);
    D = max(D, 0); % Safety for negative numerical noise
    scale_mat = sqrt(D) * n_sig; 
    
    ellipse_pts = V * scale_mat * unit_circle;
    final_pts = ellipse_pts + mu;
    
    % --- 3. Plotting ---
    % Only plot if we are holding the plot (standard usage) 
    % or if no output is requested (void call)
    if ~nargout || ishold
        plot(final_pts(1,:), final_pts(2,:), 'Color', col, 'LineStyle', style, 'LineWidth', lw);
        h_marker = plot(mu(1), mu(2), 'Marker', marker, 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', ms, 'LineWidth', 2);
    end

    % --- 4. Smart Output ---
    % If the user asks for a return value (e.g., pts = draw_gaussian...),
    % return the POINTS. If they want the handle, they usually don't assign it 
    % inside a math loop.
    if nargout > 0
        output = final_pts; % Returns 2x100 matrix of ellipse points
    else
        % If no output requested, just return the marker handle silently
        output = h_marker; 
    end
end