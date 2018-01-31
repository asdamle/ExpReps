function [x, err, idx] = svd_coneigen(H,acc)
% svd method for finding con-eigen pairs
% Anil Damle
% 
% Specifically, this function solves Hx = s*conj(x) where s is picked to be 
% below the target relative accuracy
% inputs:
% H: a NxN Hankel matrix
% acc: the desired relative accuracy
% 
% outputs:
% x: the Nx1 coneigen vector
% err: The error in the solution to the coneigenvalue problem i.e. 
[~,s,v] = svd(H);
u = conj(v); %to deal with complex case?
s = diag(s);
s2 = abs(s)/max(abs(s));
idx = find(s2<=acc, 1 );

if isempty(idx)
	fprintf('No coneigen value is small enough, acc will be %e\n',s2(end))
    x = v(:,end)+conj(u(:,end));
    idx = length(s);
else
    x = v(:,idx)+conj(u(:,idx));
end

%remove excess computation, the error should be fine, set to null to avoid
%compatibility issues with existing code. 
c = fft(flip([H(end,2:end), 0, H(1,:)]));
err = norm(fast_H_vec(c,x)-s(idx)*conj(x));
