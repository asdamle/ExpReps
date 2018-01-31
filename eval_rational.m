function fr = eval_rational(zm,wm,x)

% evaluates a rational approximation at points in x, should be valid on [0,1)
% because these are constructed using the DFT, a periodic assumption is made


fr = zeros(1,length(x));
for k = 1:length(zm)
    fr = fr + wm(k)./(x-zm(k)) + conj(wm(k))./(x-conj(zm(k)));
end