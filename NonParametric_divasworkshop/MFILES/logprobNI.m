function logp=logprobNI(x0,y0,dx0,dy0,~,modelspec,thet,spacex)
dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
        mspec.cvfunc = @(x1,x2) modelspec.cvfunc(x1,x2,thet);
        mspec.dcvfunc = @(x1,x2) modelspec.dcvfunc(x1,x2,thet);
        mspec.ddcvfunc = @(x1,x2) modelspec.ddcvfunc(x1,x2,thet);
        cvfunc=modelspec.cvfunc;

        [dK,df,d2f,yoffset] = GPRdx(x0,y0,dx0,dy0,mspec,1);

        logp = logprob(y0-yoffset,@(theta) cvfunc(x0,x0,theta)+diag(dy0.^2)+diag(dK),thet,[]);
end
