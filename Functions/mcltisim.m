function [t_out, x_out, z_out] = mcltisim(F,G,Gamma,H,Q,R,u,xbar0,P0,kmax)

% Initialize the output vectors
t_out  = 0:kmax;
x_out = zeros(kmax+1, length(xbar0));
z_out = zeros(kmax, size(H, 1));


% Draw initial condition and save
x0 = mvnrnd(xbar0, P0, 1);
x_out(1, :) = x0;

% Generate noise
vk_vec = mvnrnd(zeros(size(Q, 1), 1), Q, kmax).';
wk_vec = mvnrnd(zeros(size(R, 1), 1), R, kmax).';

% Run the sim
x_last = x0.';
for k=1:kmax
    % Propagate state
    vk = vk_vec(:, k);
    x = F*x_last + G*u(k, :).' + Gamma*vk;
    x_out(k+1, :) = x.';

    % Take measurement
    wk = wk_vec(:, k);
    zk = H*x + wk;
    z_out(k, :) = zk;

    % For next iteration
    x_last = x;
end
