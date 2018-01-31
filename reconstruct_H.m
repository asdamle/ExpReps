function [est] = reconstruct_H(nodes,weights,xloc)
%Anil Damle
%reconstruct function from nodes and weights
%do not use if representation of fft exists
%then use reconstruct_H_fourier

%inputs 
%     nodes: nodes of representation, vector of necessary length
%     weights: weights of representation, vector of necessary length
%     xloc: locations for evaluation (normally linspace(0,1,N) or
%     linspace(0,N-1,N)

%outputs
%     est: vector of length n constructed from nodes and weights

V = exp((ones(length(xloc),1)*nodes.').*(xloc(:)*ones(1,length(nodes))));
est = V*weights(:);
est = est.';



% %%% just in case...
% N = length(xloc);
% 
% V = (ones(N,1)*nodes.').^(xloc(:)*ones(1,length(nodes)));
% 
% est = V*weights(:);
% est = est.';