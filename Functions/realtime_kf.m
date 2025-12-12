function [xbarkp1, Pbarkp1, xhatkp1, Pkp1, nu, S] = realtime_kf(Fk, Gk, Gammak, Hk, Qk, Rk, uk, zkp1, xhatk, Pk)
%REALTIME_KF Single KF step
%   Detailed explanation goes here

% Predict
xbarkp1 = Fk*xhatk + Gk*uk;
Pbarkp1 = Fk*Pk*Fk.' + Gammak*Qk*Gammak.';

% Update
nu = zkp1.' - Hk*xbarkp1;
S = Hk*Pbarkp1*Hk.' + Rk;
K = Pbarkp1*Hk.'/S; % Kalman gain
xhatkp1 = xbarkp1 + K*nu; % Updated state estimate
Pkp1 = (eye(nx) - K*Hk)*Pbarkp1*(eye(nx) - K*Hk).' + K*Rk*K.'; % Joseph form cov update
end