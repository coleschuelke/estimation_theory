function [x_out, P_out, t_out, nis] = kalman_filter(F, Gamma, H, Q, R, xhat0, P0, z)
%KALMAN_FILTER The most basic KF ever to exist


    % Pass out a time vector
    t_out = [0 repelem(1:length(z), 2)].';

    num_meas = size(z, 1)/size(H, 1);
    
    % Initialize output
    x_out = zeros(num_meas*2 + 1, length(xhat0));
    x_out(1, :) = xhat0;
    P_out = zeros([size(P0), num_meas*2 + 1]);
    P_out(:, :, 1) = P0;
    nis = 0;
    
    % Inititialize priors
    x_post = xhat0;
    P_post = P0;
    
    % Recursive Estimation
    for k=1:num_meas

        % Prediction step
        x_prior = F*x_post;
        P_prior = F*P_post*F.' + Gamma*Q*Gamma.';
        
        % Store the prediction 
        x_out(2*k, :) = x_prior;
        P_out(:, :, 2*k) = P_prior;
    
        % Update step
        nu = z(k, :).' - H*x_prior;
        S = H*P_prior*H.' + R;
        K = P_prior*H.'/S;
        x_post = x_prior + K*nu;
        P_post = P_prior - K*S*K.';
    
        % Store the correction
        x_out(2*k + 1, :) = x_post;
        P_out(:, :, 2*k + 1) = P_post;

        if nargout > 3
            nis = nis + nu.'*inv(S)*nu;
        end
    end


end