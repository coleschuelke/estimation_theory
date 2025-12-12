clear variables;
close all;
clc;

addpath('..\Functions\');

xbar = [5; 1; 2];
Pxx = diag([7, 4, 9]);
samples = 1000000;

x = zeros(length(xbar), samples);

for i=1:samples
    x(:, i) = draw_gaussian(xbar, Pxx);
end

xbar_exp = mean(x, 2)

P = zeros(3, 3, samples);
for j=1:samples
    P(:, :, j) = (x(:, j)-xbar)*(x(:, j)-xbar).';
end

Pxx_exp = mean(P, 3)

% Can do a chi-squared test to check cov, but looks good from here