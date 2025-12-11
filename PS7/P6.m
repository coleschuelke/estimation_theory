clear variables;
close all;
clc;

% Set a seed 
rng(45);

% Building the functions
function [p] = h_nl(x1, x2)
    p = [sqrt(x1^2 + x2^2); atan2(x1, x2)];
end

function [x] = hinv_nl(rho, theta)
    x = [(rho^2 - ((rho^2) / (tan(theta)^2 + 1)))^0.5; ((rho^2) / (tan(theta)^2 + 1))^0.5];
end

h_handle = @(x1, x2) h_nl(x1, x2);
hinv_handle = @(rho, theta) hinv_nl(rho, theta);

% Building the linearization
function [dhdrho] = dhinv1_drho(rho, theta)
    dhdrho = (((rho^2 * tan(theta)^2) / (tan(theta)^2 + 1))^0.5) / rho;
end

function [dhdtheta] = dhinv1_dtheta(rho, theta)
    dhdtheta = (cot(theta) * csc(theta)^2 * ((rho^2) / (cot(theta)^2 + 1))^1.5) / rho^2;
end

function [dhdrho] = dhinv2_drho(rho, theta)
    dhdrho = (((rho^2) / (tan(theta)^2 + 1))^0.5) / rho;
end

function [dhdtheta] = dhinv2_dtheta(rho, theta)
    dhdtheta = - (tan(theta)*sec(theta)^2*((rho^2)/(tan(theta)^2 + 1))^1.5) / rho^2;
end

function [H] = H_lin(rho, theta)
    H = [dhinv1_drho(rho, theta), dhinv1_dtheta(rho, theta); dhinv2_drho(rho, theta), dhinv2_dtheta(rho, theta)];
end

% Setting up the problem
R = diag([100^2, (0.5*pi/180)^2]);
H = H_lin(10^5, 45*pi/180);

% Linearized cartesian covariance
Rc = H.'*R*H;

% Draw the realizations
polar_noise = mvnrnd([0;0], R, 1000).';
polar_meas = [10^5; 45*pi/180] + polar_noise;
cart_meas = zeros(size(polar_meas));

for i=1:length(cart_meas)
    cart_meas(:, i) = hinv_nl(polar_meas(1, i), polar_meas(2, i));
end

cart_avg = mean(cart_meas, 2)
cart_0 = hinv_nl(10^5, 45*pi/180)

% Calculate the cartesian noise
cart_noise = cart_meas - cart_0;

% Hypothesis test for the covariance matrix
W = 0;
for i=1:length(cart_meas)
    W = W + cart_0.'*Rc*cart_0;
end

ub = chi2inv(0.995, 2000)
lb = chi2inv(0.005, 2000)
W

residual_prob = chi2pdf(W, 2000)