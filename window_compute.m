%compute windows for use with exponential approximation
function [window] = window_compute(n,order)
%input: n is desired window length
%input: order defines sharpness of window transition
%output: window is a vector of length n that defines window

%computes |m_0| from wavelet analysis and scales to use as window

%create window, sourced from wavelet class notes 
moments = order;
syms y z 
p = 0;
tol = 10e-9;
for k = 0:moments-1
    p = p+nchoosek(moments-1+k,k)*(y^k);
end

g = subs(p,y,(-(1-z)^2)/(4*z));

GG = simplify((z^(moments-1))*g);

GG_poly = sym2poly(GG);

GG_roots = roots(GG_poly);

used = ones(size(GG_roots));    %used is 1 if root has not been used
TT = 1;
while sum(used)>0
    a = find(used==1);
    temp = GG_roots(a(1));
    
    if imag(temp) == 0;
        used(a(1)) = 0;
        
        temp_compare = abs(GG_roots-1/temp);
        used(temp_compare<tol) = 0;
        
        TT = TT*(z-temp); 
    else
        used(a(1)) = 0;
        
        temp_compare = abs(GG_roots-conj(temp));
        used(temp_compare<10e-9) = 0;
        
        temp_compare = abs(GG_roots-1/temp);
        used(temp_compare<10e-9) = 0;
        
        temp_compare = abs(GG_roots-conj(1/temp));
        used(temp_compare<10e-9) = 0;
        
        TT = TT*(z^2-2*z*real(temp)+abs(temp)^2);
    end
    
end

TT1 = subs(TT,z,1);
T = ((((z+1)^moments)*sqrt(2)*TT)/(TT1*(2^moments)));
coeffs = sym2poly(T);

syms s;
Test = subs(T,z,exp(-2*pi*1i*s));
temp = linspace(-.5,.5,n);
window = (abs(double(subs(Test,s,temp)))./sqrt(2)).^2;


