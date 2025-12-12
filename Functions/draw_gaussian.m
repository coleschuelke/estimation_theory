function [x] = draw_gaussian(xbar, Pxx)
    T = chol(Pxx);

    x = T * randn(length(xbar), 1) + xbar;
end