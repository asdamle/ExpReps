function [U,S,V] = rand_svd_hankel(h,l)

%computes a small svd "rank l" randomly using algorithm 5.1 from 
% Halko, Martinsson and Tropp 2011. Here it is written specifically for
%Hankel matrices

%could be faster by 

[Q] = rand_range_hankel(h,l);

N = (length(h)+1)/2;
h = conj(h(:).');
c = fft(flip([h(N+1:end), 0, h(1:N)],2));

B = fast_H_vec(c,Q)';

[~, S, V] = svd(B,'econ');
U = conj(V); %this is to deal with the complex case