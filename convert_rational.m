function [zm, wm] = convert_rational(nodes,weights,N)
%converts a decaying exponential representation in Fourier as computed by
%representH_fourier into a rational approximation valid on [0 1];
%need to know N to deal with DFT stuff

% 2N+1 samples (with N even) should have been used on [0,1) to generate the
% initial signal

zm = -1*nodes/N;
wm = weights/(2*N+1);
zm(imag(zm)<0)=zm(imag(zm)<0)+2*pi*1i;

wm = -1*wm/(2*pi*1i);
zm = zm/(2*pi*1i);



