function [x] = draw_gaussian(xbar, Pxx, nsamples)
    T = chol(Pxx);
    x = T * randn(length(xbar), nsamples) + xbar;
end