%fermi dirac test
close all

beta = 10000;
f = @(x) 1./(1+exp(beta*(x)));
acc = 1e-12;
% find where we have to sample the function
x_fine = linspace(-10,10,100000);
f_fine = f(x_fine);
x_min = find(abs(f_fine-1) > acc/10,1);
x_max = find(abs(f_fine) < acc/10,1);
xl = max(abs([x_fine(x_min) x_fine(x_max)]));

x  = linspace(-xl,xl,1000);
fd = [zeros(1,250) flip(f(x)) ones(1,251) f(x) zeros(1,250)];
fd = [flip(f(x)) 1 f(x)];
% fd = [1 f(x)];
% fd = [1 f(x)];

N = length(fd);
temp = fft(fd);
temp(1:(N+1)/2);
% semilogy(abs(temp))

[nodes, weights, xloc] = represent_Hrand_fourier(fd,acc,30);

% idx = imag(nodes)>0;
% nodes = nodes(idx);
% weights = weights(idx);

[est] = reconstruct_H_fourier(nodes,weights,xloc);


plot_poles(nodes,weights,xloc)
figure
semilogy(abs(est-fd))


[zm, wm] = convert_rational(nodes,weights,N);
fr = eval_rational(zm,wm,x);
