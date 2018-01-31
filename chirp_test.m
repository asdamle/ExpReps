clc
close all
N = 500;
x = linspace(0,1,2*N+2);
x = x(1:end-1);
acc = 1e-5;

f = sin(73*pi*x.^2) + sin(29*pi*x);

w = erfc(50*(x-0.8)).*erfc(-50*(x-0.2))/4;
f = f.*w;

[nodes, weights, xloc] = represent_H_fourier(f,acc);

[est] = reconstruct_H_fourier(nodes,weights,xloc);
norm(est-f,'inf')
figure();
plot_poles(nodes,weights,xloc)

figure; semilogy(abs(f-est));

% for some reason, 
% the conversion to a rational representation does odd things here...
% [zm, wm] = convert_rational(nodes,weights,N);
% fr = eval_rational(zm,wm,x);
% 
% figure; plot(fr); hold on; plot(f,'r');
% figure; semilogy(abs(f-fr));

