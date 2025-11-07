function [x_out,P_out] = kalman_filter(F, Gamma, H, Q, R, xbar0, P0, z)
%KALMAN_FILTER The most basic KF ever to exist

    num_meas = size(z, 1);
    
    % Initialize output
    x_out = zeros(length(xbar0), num_meas*2 + 1);
    x_out(:, 1) = xhat0;
    P_out = repmat(zeros(size(P0)), 1, length(thist)*2 + 1);
    P_out(:, 1:2) = P0;
    
    % Inititialize priors
    x_post = xhat0;
    P_post = P0;
    
    % Recursive Estimation
    for k=1:num_meas
        
        % Prediction step
        x_prior = Fk*x_post;
        P_prior = Fk*P_post*Fk.' + Gammak*Qk*Gammak.';
        
        % Store the prediction 
        x_out(:, 2*k) = x_prior;
        P_out(:, 4*k-1:4*k) = P_prior;
    
        % Update step
        nu = z(k, :) - Hk*x_prior;
        S = Hk*P_prior*Hk.' + Rk;
        K = P_prior*Hk.'/inv(S);
        x_post = x_prior + K*nu;
        P_post = P_prior - K*S*K.';
    
        % Store the correction
        x_out(:, 2*k + 1) = x_post;
        P_out(:, 4*k+1:4*k+2) = P_post;
    end


end