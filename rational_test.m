N = 1000;
x = linspace(0,1,2*N+2);
x = x(1:end-1);

f = zeros(1,length(x));
f(x>=1/5) = 1;
f(x>=2/5) = (3-5*x(x>=2/5)).^2;
f(x>=3/5) = 40*((3-5*x(x>=3/5)).^2).*((4-5*x(x>=3/5)).^3);
f(x>=4/5) = 0;


[nodes, weights, xloc] = represent_H_fourier(f,1e-8);
[est] = reconstruct_H_fourier(nodes,weights,xloc);
semilogy(abs((est-f)))
plot_poles(nodes,weights,xloc)


[zm, wm] = convert_rational(nodes,weights,N);
fr = eval_rational(zm,wm,x);

figure; plot(fr); hold on; plot(f,'r');
figure; semilogy(abs(f-fr));