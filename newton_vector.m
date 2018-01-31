function [x] = newton_vector(f,df,x,tol)

max_iter = 1000;
fx = polyval(f,x);
for k = 1:max_iter    
    dfx = polyval(df,x);
    x = x - fx./dfx;
    fx = polyval(f,x);
    if max(abs(fx)) < tol
        break;
    end
end
x = x(abs(fx) < tol);

%filter out duplicates, in a crude manner
for k = 1:length(x)
    x = [x(k) x(abs(x(k)-x)>tol)];
    if k == length(x)
        break;
    end
end