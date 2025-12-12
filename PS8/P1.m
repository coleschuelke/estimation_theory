clear variables;
close all;
clc;

% Sim Parameters
kmax = 5000;

% Model Parameters
c1 = 0.1; c2 = 0.5; c3 = 1;
hz1 = 1; hz2 = 2; hz3 = 3;

a1 = c1/hz1; b1 = 1/hz1; 
a2 = c2/hz2; b2 = 1/hz1;
a3 = c3/hz3; b3 = 1/hz3;

% System Matrices
A1 = [0, 1; 0, -a1];
A2 = [0, 1; 0, -a2]; 
A3 = [0, 1; 0, -a3];
B1 = [0; b1];
B2 = [0; b2];
B3 = [0; b3];

Gamma = eye(2);
H = [1, 0];
R = 0.1;
Q = 0.001*diag([0.1, 1]);
dt = 0.1;

Fk1 = expm(A1*dt);
[~, y] = ode45(@(t, y) reshape(expm(A1*(dt - t)), 4, 1), [0, dt], zeros(4, 1));
Gk1 = reshape(y(end, :), 2, 2)*B1;
Fk2 = expm(A2*dt);
[~, y] = ode45(@(t, y) reshape(expm(A2*(dt - t)), 4, 1), [0, dt], zeros(4, 1));
Gk2 = reshape(y(end, :), 2, 2)*B2;
Fk3 = expm(A3*dt);
[~, y] = ode45(@(t, y) reshape(expm(A3*(dt - t)), 4, 1), [0, dt], zeros(4, 1));
Gk3 = reshape(y(end, :), 2, 2)*B3;

% Stack the multiple models
Fk = cat(3, Fk1, Fk2, Fk3);
Gk = cat(3, Gk1, Gk2, Gk3);
Gammak = cat(3, Gamma, Gamma, Gamma);
Hk = cat(3, H, H, H);
Rk = cat(3, R, R, R);
Qk = cat(3, Q, Q, Q);

% Create input
u = zeros(kmax, 1);
u(501:600, :) = 7;

xbar0 = [0;0.1];
P0 = Q*100;


% Simulate Truth
[t_truth, x_truth, z] = mcltisim(Fk1, Gk1, Gamma, H, Q, R, u, xbar0, P0, kmax);

% Run the MM filter
[t_mm, x_mm, P_mm, mu_mm] = static_mm_filter(Fk, Gk, Gammak, Hk, Qk, Rk, u, z, xbar0, P0);


%% Plotting
figure;
hold on;
plot(t_truth, x_truth(:, 1), 'g');
x1_think = squeeze(x_mm(1, end, :));
plot(t_mm, x1_think, 'r')

figure;
plot(t_mm, mu_mm);