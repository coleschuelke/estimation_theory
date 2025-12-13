clear variables;
close all;
clc;

addpath('..\Functions\');
addpath('codeAndDataForStudents\');

% Load data
load('problem4data.mat');
num_meas = length(sonar);
nparts = 5000;

% Initialize output
x_out = zeros(num_meas, 3);
P_out = zeros(3, 3, num_meas);

% Set up problem matrices
Q = diag([0.1, 5*pi/180].^2);
R = eye(3);

% Draw initial particles;
x_values = unifrnd(minx, maxx, 1, nparts);
y_values = unifrnd(miny, maxy, 1, nparts);
heading = unifrnd(0, 2*pi, 1, nparts);
weights = 1/nparts*ones(1, nparts);
particles = [x_values; y_values; heading; weights]; % Each col is a particle

% Loop through the measurements
for k=1:num_meas
    % Extract the control input
    u = encoder(k).u; % 2x1
    process_noise = draw_gaussian(zeros(2, 1), Q, nparts);
    u_noisy = u + process_noise;

    % Propagate particles based on the noisy control input
    particles(1, :) = particles(1, :) + u_noisy(1, :) .* cos(particles(3, :) + u_noisy(2, :));
    particles(2, :) = particles(2, :) + u_noisy(1, :) .* sin(particles(3, :) + u_noisy(2, :));
    particles(3, :) = particles(3, :) + u_noisy(2, :);

    % Compute zbar for each particle
    % (x-y)^2 = x^2 + y^2 - 2xy
    p_sq = sum(particles(1:2, :).^2, 1).';
    b_sq = sum(beacons.^2, 2).'; 
    d_sq = max(p_sq + b_sq - 2 * (particles(1:2, :).' * beacons.'), 0);
    
    % Euclidean distance
    d = sqrt(d_sq);
    d = mink(d, 3, 2).'; % 3xnparts;

    % Update the weights and normalize
    z = sonar(k).z;
    resids = d - z;
    lw = 0.5*sum(resids .* (inv(R) * resids), 1) + particles(4, :);
    lw_max = max(lw);
    w = exp(lw - lw_max);
    w = w / sum(w);
    particles(4, :) = w;

    % Form estimate
    xhat = particles(1:3, :) * particles(4, :).';
    xhat(3, :) = mod(xhat(3, :), 2*pi); % Handle wraparound
    % Calculate covariance TODO

    % Save estimates
    x_out(k, :) = xhat;

    % Resample (if needed)
    N_eff = 1 / sum(particles(4, :).^2);

    if N_eff < nparts/2
        l=1;
        while l < nparts
            eta = unifrnd(0, 1);
            m = 1;
            w_sum = 0;
            while w_sum <= eta && m <= nparts
                w_sum = w_sum + particles(4, m);
                m = m+1;
            end
            xm = particles(1:3, m-1);
            particles(:, l) = [xm; 1/nparts];
            l = l+1;
        end
        
    end


end