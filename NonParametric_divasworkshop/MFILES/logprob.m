function logp=logprob(y0,traincv,x,basisX)

    errorflags=0;
    
    tcv=traincv(x);
    try
        L=chol((tcv),'lower');
        %	catch
        %			disp('fall back');
        %			L=cholcov(tcv,'lower');
        %			keyboard
        %		end
        alfa=L'\(L\y0);
        doChol=1;
    catch
        disp('Not positive definite!')
        errorflags=-1;
        doChol = 0;
        [m,n] = size(tcv);
        [U,S,V] = svd(tcv,0);
        s = diag(S);
        tol = max(m,n) * eps(max(s));
        r = sum(s > tol);
        invtraincv = V(:,1:r)*diag(s(1:r).^-1)*U(:,1:r)';      
        alfa = invtraincv * y0;
        %	            disp('SVD')
    end

    logpterms(1) = -.5*abs(y0'*alfa);
    logpterms(3) = -.5*length(y0)*log(2*pi);
    % logpterms(2) = -.5 * log(det(traincv));
    if doChol
        logpterms(2) = - sum(log(diag(L)));
    else
        logpterms(2) = -.5 * sum(log((s(1:r))));
    end
    logp=sum(logpterms);
    if errorflags == -1
        logp = -1e20;
    end
    
end
