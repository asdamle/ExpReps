function [ b ] = fast_H_vec(c,x)
%fast application of a hankel matrix to a set of vectors x
x0 = [flip(x);zeros(size(x))];

% to save time c must be defined in the calling function as
% c = fft(flipdim([H((size(H,2)+1)/2+1:end), 0, H(1:(size(H,2)+1)/2)],2));
c = c(:);
b = fft((c*ones(1,size(x,2))).*(ifft(x0,[],1)));
b = b(1:size(x,1),:);

