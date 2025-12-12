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

    xbar_t = (lambda/(nx+lambda))*chi_t(:, 1) + sum((1/(2*(nx+lambda)))*chi(:, 2:end), 2);
    Pxx_t = (lambda/(nx+lambda) + 1 - alpha^2 + beta)*(chi(:, 1) - xbar_t)*(chi(:, 1) - xbar_t).';
    for k=1:2*npts
        Pxx_t = Pxx_t + 2*nx*(chi(:, k) - xbar_t)*(chi(:, k) - xbar_t).';
    end
end

nl_handle = @(x) h_nl(x(1), x(2));

[cbar_uti, Pcc_uti] = ut([rbari; thetabari], P_pi, nl_handle, 1e-3, 2, 0)
[cbar_utii, Pcc_utii] = ut([rbarii; thetabarii], P_pii, nl_handle, 1e-3, 2, 0)