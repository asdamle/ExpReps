function [nodes, weights, xloc] = represent_H(f,acc)
% Anil Damle
% develop hankel based representation for f 
% do not use directly on fft of signals, for that use represent_H_fourier
% inputs:
%   f: signal to be represented, vector length 2n-1
%   acc: accuracy of representation, scalar
%
% outputs: 
%     nodes: nodes of representation, vector of necessary length
%     weights: weights of representation, vector of necessary length


%compute coneigensystem
f = f(:).';
N = length(f);
xloc = linspace(0,1,N);

%svd method for finding coneigenvalue/vector
[x, ~] = svd_coneigen(create_hankel(f),acc);

%find roots of con-eignepolynomial inside unit circle
rx = roots(flip(x(:).'));
rx = rx(abs(rx)<1);

% add if only real poles are desired...
% rx = rx(abs(imag(rx))<1e-8);
% rx = real(rx);

% maybe clean up with Newtons method...
% rxN = newton_vector(flip(x(:).'),polyder(flip(x(:).')),rx,acc/10);
% rx = rxN(abs(rxN)<1);

nodes = (N-1)*log(rx(:));


%construct vandermonde matrix
V = exp((ones(length(xloc),1)*nodes.').*(xloc(:)*ones(1,length(nodes))));
weights = V\f(:);
if ~isempty(find(abs(weights)<acc,1))
	nodes = nodes(abs(weights)>=acc);
	nodes = nodes(:);
	%recompute without un-necessary nodes
	V = exp((ones(length(xloc),1)*nodes.').*(xloc(:)*ones(1,length(nodes))));
    weights = V\f(:);
end
    
