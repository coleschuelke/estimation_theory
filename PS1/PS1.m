close all;
clc;
clear;

vec_1 = rand(100000, 1);
vec_2 = (rand(100000, 1) + rand(100000, 1))/2;

vec_100 = zeros(100000, 1);

for i=1:100
    sample = rand(100000, 1);
    vec_100 = vec_100 + sample;
end
vec_100 = vec_100 / 100;

subplot(3, 1, 1)
histogram(vec_1, 'Normalization','probability')

subplot(3, 1, 2)
histogram(vec_2, 'Normalization','probability')

subplot(3, 1, 3)
histogram(vec_100, 'Normalization','probability')

disp('Done')