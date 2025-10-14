clear all;
close all;
clc;

N_req = -1;
for N=1:10000
    sigma = sqrt(2*10^4/N);
    if normcdf(80, 100, sigma) < 0.05
        N_req = N;
        break
    end
end

N_req