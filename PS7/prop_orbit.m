function [xbigdot] = prop_orbit(t, x, u, v)
%PROP_ORBIT odefun to propagate x, f, and gamma
    % Break out useful quantities
    rs = x(1:3);
    vs = x(4:6);
    r1 = rs(1);
    r2 = rs(2);
    r3 = rs(3);

    % Constants and shorthand
    nx = 6;
    nu = 3;
    nv = 3;
    mu = 3.986005e14;
    rs3 = norm(rs)^3;
    rs5 = norm(rs)^5;

    % Dynamics function
    rsdot = vs;
    vsdot = (-mu/rs3) * rs + u + v;
    xdot = [rsdot; vsdot];

    % Partial wrt x
    ddx = zeros(nx);
    ddx(1:3, 4:6) = eye(3);
    ddx(4:6, 1:3) = [mu*(2*r1^2 - r2^2 - r3^2)/rs5, 3*mu*r1*r2/rs5, 3*mu*r1*r3/rs5;
        3*mu*r1*r2/rs5, mu*(2*r2^2 - r1^2 - r3^2)/rs5, 3*mu*r2*r3/rs5;
        3*mu*r1*r3/rs5, 3*mu*r3*r3/rs5, mu*(2*r3^2 - r1^2 - r2^2)/rs5];
    
    % Partial wrt v
    ddv = zeros(nx, nv);
    ddv(4:6, :) = eye(3);

    % Reshape and combine for propagation
    xbigdot = [xdot; reshape(ddx, [], 1); reshape(ddv, [], 1)];
end