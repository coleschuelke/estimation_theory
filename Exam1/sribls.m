function [xhat,Rotilde] = sribls(Hprime,zprime,R)
    % sribls : Square-root-information-based least squares routine. Given the
    % linear measurement model
    % zprime = Hprime*x + w, w ~ N(0,R)
    %
    % sribls returns xhat, the least squares (and maximum likelihood)
    % estimate of x, and Rotilde, the associated square root
    % information matrix.
    %
    % INPUTS
    % Hprime ----- nz-by-nx measurement sensitivity matrix.
    % zprime ----- nz-by-1 measurement vector.
    % R ---------- nz-by-nz measurement noise covariance matrix.
    %
    % OUTPUTS
    % xhat ------- nx-by-1 maximum likelihood estimate of x.
    % Rotilde ---- nx-by-nx square root information matrix, where P =
    % inv(Rotilde'*Rotilde) is the estimation error covariance matrix.
    %
    %+----------------------------------------------------------------------------
    Ra = chol(R);
    RaInv = inv(Ra);
    RaT = Ra';
    RaInvT = inv(Ra');

    z = inv(Ra')*zprime;
    H = RaInvT*Hprime;
    
    [Qtilde, Rtilde] = qr(H);

    ztilde = Qtilde'*z;

    dim = min(size(H));

    Rtilde0 = Rtilde(1:dim, 1:dim);
    ztilde0 = ztilde(1:dim);

    xhat = ztilde0\Rtilde0;
    Rotilde = Rtilde0;
end