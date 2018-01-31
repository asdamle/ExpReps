N = 250;
x = linspace(0,1,2*N+2);
x = x(1:end-1);

% discontinuity if viewed as periodic makes things a bit odd
nx = 0.57+0.05*1i;
f = 1./(x-nx)+1./(x-conj(nx));
w = exp(-70*(x-.5).^2);
f = f.*w;

[nodes, weights, xloc] = represent_H_fourier(f,1e-12);
[est] = reconstruct_H_fourier(nodes,weights,xloc);
plot_poles(nodes,weights,xloc)

[zm, wm] = convert_rational(nodes,weights,N);
fr = eval_rational(zm,wm,x);

figure; plot(fr); hold on; plot(f,'r');
figure; semilogy(abs(f-fr));