clear variables;
close all;
clc;

format long;

case_1 = scalar_ratio(10^4, 10^5)
case_2 = scalar_ratio(10^6, 10^5)

% Vector case
x11 = -10^4;
x12 = 10^5;
x21 = 10^4;
x22 = 10^5;

H = [-x11 / (sqrt(x11^2 + x12^2)), -x12 / (sqrt(x11^2 + x12^2)); -x21 / (sqrt(x21^2 + x22^2)), -x22 / (sqrt(x21^2 + x22^2))];

GDOP = trace(inv(H'*H)*(H.'*H)*inv(H.'*H).')

function [ratio] = scalar_ratio(xbar, y)
    ratio = sqrt((xbar^2 + y^2)/(xbar^2));
end