%power method for hankel matricies
function [x u E iter] = power_method_H(A, acc, Tol)
%compute the con-eigenvectors/con-eigenvalues of a hankel matrix defined by signal
%A, with Tol and compute as many as necessary such that u/umax<acc where
%umax is the largest eigenvalue and u is the current

%inputs: A is the signal that defines the hankel matrix dimension must be
%(2^n) - 1
%        acc defines the number of eigenvalues to compute, i.e compute
%        eigenvalue/vector pairs until u/umax<acc where umax is the largest
%           eigenvalue and u is the current 
%       Tol gives the error to use when computing eignevalue/vector pairs

%outputs: E is the error at the last step
%         iter is the number of iterations needed
%         u is a sorted list of con-eigenvalues from largest in magnitude to
%         smallest in magnitude
%         x are the con-eigenvectors, if A is real, then c is the
%         con-eigenvector. If A is complex then x = [conj(v); v] where v is
%         the con-eigenvector, computed by finding eigenvectors of [0 H;
%         conj(H) 0]



%needed constants and set max iterations
n = (size(A,2)+1)/2;
max_iter = 10^12;

%compute fft of A from circulant extened hankel matrix to use in matrix
%multiply, also compute the fft of conj(A) to use when multiplying by
%\conj(H)
c = fft(flipdim([A((size(A,2)+1)/2+1:end), 0, A(1:(size(A,2)+1)/2)],2));
A_conj = conj(A);
c_conj = fft(flipdim([A_conj((size(A_conj,2)+1)/2+1:end), 0, A_conj(1:(size(A_conj,2)+1)/2)],2));
if isreal(A) %for real H
    
    %set up initial values
    x0 = eye(n,1);
    E = 0;
    k = 1;
    x = x0/norm(x0,inf);
    u = NaN;
    
    for k = 2:max_iter
        
        %save current u and x for use
        old_u = u;
        p_previous = x;
        
        %A*x
        for temp = 1:size(x,2)
            y(:,temp) = real(fast_H_vec(c,x(:,temp)));
        end
        
        %qr on x
        [x not_needed] = qr(y,0);        
        
        %compute eigenvalues
        for temp = 1:size(x,2)
            u(temp) = x(:,temp).'*y(:,temp);
        end
        %compute error
%         E = max(abs(old_u-u));
         E = norm(real(fast_H_vec(c,x(:,end)))-u(end)*x(:,end));

%check break condition
        if E < Tol
            
            %are more vectors needed to meet acc constraint? if yes add a
            %new vector and continue, else break
            if abs(u(end)/u(1))<acc || length(u) == n
                E2 = 0;%max(abs(p_previous(:,end)-x(:,end)));
                if E2 <Tol
                     break;
                end
            else
                x = [x,randn(n,1)];
                u = [u NaN];
            end
        end
    end
    iter = k;
    
    z = x./p_previous;
    theta = angle(z(1,:));
    for temp = 1:size(x,2)
        x(:,temp) = exp(-i*theta(temp)/2)*x(:,temp);
    end
    
    
    
    
else %complex H
    %initial set up
    x0 = randn(n,1)+i*randn(n,1);
    x0 = [conj(x0);x0];
    E = 0;
    k = 1;
    x = x0/norm(x0,inf);
    u = NaN;

    %loop
    for k = 2:max_iter
        
        
        old_u = u;
        p_previous = x;
        x1 = x(1:n,:);
        x2 = x(1+n:2*n,:);
        
        %multiply by matrix
        for temp = 1:size(x,2)
            y(:,temp) = [fast_H_vec(c,x2(:,temp));fast_H_vec(c_conj,x1(:,temp))]; 
        end
        
        %qr on vectors
        tic
        [x not_needed] = qr(y,0);
        timeqr = toc;
        %maintain [conj(v); v] format
        x1 = x(1:n,:);
        x = [x1;conj(x1)];

        %compute eigenvalues
        for temp = 1:size(x,2)
            u(temp) = x(:,temp)'*y(:,temp)/(x(:,temp)'*x(:,temp)); %check .' vs '
        end
        
        %apply adjustment for coneigenvectors
        x1 = x(1:n,end);
        x2 = x(1+n:2*n,end);
        z = (x./p_previous);
        %angle rotation
        theta = angle(z(1,end));

        %use angle and adjust last con-eigenvectors (all that is needed for
        %error computation)
        %
            x(:,end) = [exp(-i*theta/2)*x1(:,end);exp(i*theta/2)*x2(:,end)];
        
        x1 = x(1:n,end);
        x2 = x(1+n:2*n,end);
        %check error
%         E2 = sum(abs(old_u-u));
        E = norm((fast_H_vec(c,x2(:,end)))-u(end)*x1(:,end),inf);
        p_end = x;
        
        %are more vectors needed to meet acc constraint? if yes add a
        %new vector and continue, else break
        if E < Tol
            if abs(u(end)/u(1))<acc || length(u) == 2*n-1
                break;
            else
                sprintf('%d singular values computed, acc %d. QR took %d and %d iterations needed',length(u),abs(u(end)/u(1)),timeqr,k)
                xnew = randn(n,1)+i*randn(n,1);
                xnew = [conj(xnew);xnew];
                num_new = min(20,2*n-length(u)-1);
                x = [x,xnew*ones(1,num_new)];
                u = [u NaN*ones(1,num_new)];
            end
        end
    end
    
    %now adjust to get actual con-eigenvectors
    %of H
    iter = k;
    x1 = x(1:n,:);
    x2 = x(1+n:2*n,:);
    z = (x./p_previous);
    %angle rotation
    theta = angle(z(1,:));
    
    %use angle and adjust con-eigenvectors 
    for temp = 1:size(x,2)
        x(:,temp) = [exp(-i*theta(temp)/2)*x1(:,temp);exp(i*theta(temp)/2)*x2(:,temp)];
    end
end


