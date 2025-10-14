clear all;
close all;
clc;

% Pick an initial x_g
x_g = [-5; 5];

% Update x_g 
update = 1;
while update > 0.001
    x = x_g;
    f = [x(1) + x(2) + x(1)*x(2) + 5;
    x(1)^2 + 2*x(2) - x(2)^2 - 2];
    dfdx = [1+x(2), 2*x(1);
    1+x(1), 2-2*x(2)];

    x_g = x_g - dfdx\f;
    update = abs(x_g - x);
end

x_g