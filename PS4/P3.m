clear variables;
close all;
clc;

% Load measurements
thist =[0; 0.1000; 0.2000; 0.3000; 0.4000; 0.5000;
0.6000; 0.7000; 0.8000; 0.9000; 1.0000];

zprime = [7.7969; 1.4177; -3.0970; -7.6810; -9.8749; -6.1828;
-0.8212; 4.5074; 8.2259; 9.5369; 6.2827];

% Set up the problem
off_diag = 0.5 * ones([1, length(zprime)-1]);
R = eye(size(zprime, 1)) + diag(off_diag, 1) + diag(off_diag, -1);
Ra = chol(R);
Rait = inv(Ra.');
function h = hprime(x, t)
    h = zeros(length(t), 1);
    for j=1:length(t)
        h(j) = x(1)*cos(x(2)*t(j) + x(3));
    end
end
function H = Hprime(x, t)    
    H = zeros(length(t), length(x));
    for j=1:length(t)
        H(j, :) = [cos(x(2)*t(j) + x(3)), -t(j)*x(1)*sin(x(2)*t(j) + x(3)), -x(1)*sin(x(2)*t(j))];
    end
end

% Gauss Newton Method
x_g = [1 6 0].';
convergence_thresh = 0.002;
dx = 1;
while dx > convergence_thresh
    % Calculate dx
    Hx = Rait*Hprime(x_g, thist);
    hx = Rait*hprime(x_g, thist);
    dx = (Hx.'*Hx)\Hx.'*(Rait*zprime-hx);
    % Update guess
    x_g = x_g + dx;
end

x_g

Hfinal = Rait*Hprime(x_g, thist);
Pxx = Hfinal.'*Hfinal