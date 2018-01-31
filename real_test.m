x = linspace(0,1,201);
node1 = exp(-5);
node2 = exp(-11);
test_exp = 3*node1.^x-2*node2.^x;
[nodes, weights] = represent_Hrand_real(test_exp,1e-10,2);
[est_exp] = reconstruct_H(nodes,weights,linspace(0,1,201));

N = 250;
x = linspace(0,1,2*N+2);
x = x(1:end-1);

% discontinuity if viewed as periodic makes things a bit odd
nx = .57;
f = 1./(x-nx) + 1./(x-(1-nx));


[nodes2, weights2, xloc] = represent_Hrand_fourier(f,1e-10,2);
[est] = reconstruct_H_fourier(nodes2,weights2,xloc);
plot_poles(nodes2,weights2,xloc)

[zm, wm] = convert_rational(nodes2,weights2,N);
fr = eval_rational(zm,wm,x);
% 
% figure; plot(fr); hold on; plot(f,'r');
% figure; semilogy(abs(f-fr));