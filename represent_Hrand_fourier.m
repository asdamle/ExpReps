function [nodes, weights, xloc] = represent_Hrand_fourier(f,acc,l)
%Anil Damle
%develop hankel based representation for f through use of fft 
%
%   f: signal to be represented, vector length 2n+1 w/ n even
%   acc: accuracy of representation, scalar

%outputs 
%     dc: dc term of fft of f
%     nodes: nodes of representation, vector of necessary length
%     weights: weights of representation, vector of necessary length

N = length(f);
temp = fft(f);
[nodes, weights, xloc] = represent_Hrand(temp(1:(N+1)/2),acc,l);


