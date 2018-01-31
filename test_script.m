N = 750;
J = besselj(0,100*pi*linspace(0,1,N+1));

[ns, ws, xs] = represent_Hrand(J,1e-12,30);

% this is actually a bit less accurate in testing than the randomized version: 
% [ns, ws, xs] = represent_H(J,1e-12); 

[est] = reconstruct_H(ns,ws,xs);
norm(est-J,inf)

J2 = besselj(0,linspace(0,100*pi,3*N));
[est2] = reconstruct_H(ns,ws,linspace(0,1,3*N));
norm(est2-J2,inf)

x = linspace(0,1,201);
node1 = exp((-.05+.15*1i));
node2 = exp((-.1+.15*1i));
test_exp = (1+1i)*node1.^x+(-1i)*node2.^x;
[nodes, weights] = represent_Hrand(test_exp,1e-12,2);
[est_exp] = reconstruct_H(nodes,weights,linspace(0,1,201));


