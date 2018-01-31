function [Q] = rand_range_hankel(h,l)

%generate a n x l matrix that gives the range of the hankel matrix described by
%the vector h

%algorithm 4.1 from Halko, Martinsson and Tropp 2011 written specifically for
%Hankel matrices

N = (length(h)+1)/2;
Omega = randn(N,l);
h = h(:).';
c = fft(flip([h(N+1:end), 0, h(1:N)],2));

Y = zeros(N,l);
for k = 1:l
	Y(:,k) = fast_H_vec(c,Omega(:,k));
end

[Q , ~] = qr(Y,0);