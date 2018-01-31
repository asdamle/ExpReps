function [x, err, idx] = svd_coneigen_rand(f,acc,l)
% svd method for finding con-eigen pairs
% Anil Damle
% 
% Specifically, this function solves Hx = s*conj(x) where s is picked to be 
% below the target relative accuracy
% inputs:
% f:   signal to be represented, vector length 2n-1
% acc: the desired relative accuracy, or l the expected rank...
% 
% outputs:
% x: the Nx1 coneigen vector
% err: The error in the solution to the coneigenvalue problem i.e. 
N = length(f);
p = 10;
N2 = (N+1)/2;
while l+p < (N+1)/2;
    [u,s,v] = rand_svd_hankel(f,min(l+p,N2));
    s = diag(s);
    s2 = abs(s)/max(abs(s));
    idx = find(s2<acc,1);
    if isempty(idx)
        l = l+10;
    else
        break;
    end
end
if isempty(idx)
	fprintf('No coneigen value is small enough, acc will be %e\n',s2(end))
    x = v(:,end)+conj(u(:,end));
    idx = length(s);
else
    x = v(:,idx)+conj(u(:,idx));
end

c = fft(flip([f((N+1)/2+1:end), 0, f(1:(N+1)/2)],2));
err = norm(fast_H_vec(c,x)-s(idx)*conj(x));