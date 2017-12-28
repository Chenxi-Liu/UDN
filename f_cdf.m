function [Fx] = f_cdf(X, x)

for cntr = 1:length(x)
    Fx(cntr) = sum(X <= x(cntr))/length(X);
end