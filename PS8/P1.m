clear variables;
close all;
clc;

% Model Parameters
c1 = 0.1; c2 = 0.5; c3 = 1;
hz1 = 1; hz2 = 2; hz3 = 3;

a = c1/hz1; b = 1/hz1; % CHOOSE YOUR MODEL HERE

% System Matrices
A = zeros([0, 1; 0, -a]);
B = [0, b];
Gamma = eye(2);
H = [1, 0];
R = 0.1;
Q = 0.001*diag([0.1, 1]);
dt = 0.1;

Fk1 = expm(A1*dt);
Gk1 = ode45(@(t, y) expm(A1*(dt - t)), [0, dt], zeros()); % TODO
Fk2 = expm(A1*dt);
Gk2 = ode45(@(t, y) expm(A1*(dt - t)), [0, dt], zeros()); % TODO
Fk3 = expm(A1*dt);
Gk3 = ode45(@(t, y) expm(A1*(dt - t)), [0, dt], zeros()); % TODO

% Stack the multiple models
Fk = cat(3, Fk1, Fk2, Fk3);
Gk = cat(3, Gk1, Gk2, Gk3);
Gammak = cat(3, Gamma, Gamma, Gamma);
Hk = cat(3, H, H, H);
Rk = cat(3, R, R, R);
Qk = cat(3, Q, Q, Q);


% Simulate Truth

% Run the MM filter
[t_out, x_out, P_out] = static_mm_filter(Fk, Gammak)
