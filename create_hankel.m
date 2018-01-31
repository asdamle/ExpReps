function [H] = create_hankel(f)
%length of f must be odd
N = length(f);
H = hankel(f(1:(N+1)/2),f((N+1)/2:end));
