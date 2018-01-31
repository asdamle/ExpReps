function [est] = reconstruct_H_fourier(nodes,weights,xloc)
%Anil Damle
%reconstruct function from nodes and weights and dc term

%inputs 
%     dc: dc term of fft of f
%     nodes: nodes of representation, vector of necessary length
%     weights: weights of representation, vector of necessary length
%     xloc: locations for evaluation on the fourier side

%outputs
%     est: vector of length n constructed from nodes and weights

temp = reconstruct_H(nodes,weights,xloc);
est = ifft([real(temp(1)) temp(2:end) conj(flip(temp(2:end)))]);

