function [xkp1, Fk, GammaK] = discrete_prop(xk, uk, vk, tk, dt, propfun)
%DISCRETE_PROP Propagate one discrete time step forward
    nx = size(xk, 1);
    nu = size(uk, 1);
    nv = size(vk, 1);

    ode_options = odeset('RelTol', 1e-12); % Not positve this will work

    % Apply the controls and process noise ZOH
    odefun = @(t, x) propfun(t, x, uk, vk);
   
    % Form initial conditions
    x0 = xk;
    F0 = reshape(eye(nx), [], 1);
    Gamma0 = zeros([nx*nv, 1]);
    xbig0 = [x0; F0; Gamma0];

    % Propagate with dynamics
    [t, xbig] = ode113(odefun, [tk, tk+dt], xbig0, ode_options);

    % Pass out the quantities
    xkp1 = xbig(end, 1:nx).';
    Fk = reshape(xbig(end, nx+1:nx+nx^2).', nx, nx);
    GammaK = reshape(xbig(end, end-nx*nv+1:end).', nx, nv);


end