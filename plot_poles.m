function [] = plot_poles(nodes,weights,xloc)
%Anil Damle
%generates a plot of the reconstructed signal with poles on a plot below
%inputs
%   dc: dc term output from represent_H_fourier
%   nodes: nodes vector as output from represent_H_fourier
%   weights: weights vector as output from represent_H_fourier
%   xloc: sample locations as output from represent_H_fourier
w1_orig = abs(weights);
w1 = ceil(20*w1_orig/max(w1_orig))+3;
w1_color = 1-w1_orig/max(w1_orig);

nodes1 = -nodes/(length(xloc)-1);
nodes1(imag(nodes1)<0)=nodes1(imag(nodes1)<0)+2*pi*1i;

x1 = abs(real(nodes1)).^2;
y1 = (imag(nodes1));
y1 = y1/(2*pi);
est = reconstruct_H_fourier(nodes,weights,xloc);
figure
n = 2*length(xloc)-1;
subplot(2,1,1)
plot(linspace(0,1,n),est)
xlim([0 1])

subplot(2,1,2)
hold on
for t = 1:length(x1)
    if 1-w1_color(t) < 1e-3/max(w1_orig)
        semilogy(y1(t),x1(t),'s','MarkerSize',w1(t),'MarkerFaceColor','r')
    else
        semilogy(y1(t),x1(t),'s','MarkerSize',w1(t),'MarkerFaceColor',[w1_color(t) w1_color(t) w1_color(t)])
    end
end
set(gca,'YDir','reverse','Yscale','log')
xlim([0 1])