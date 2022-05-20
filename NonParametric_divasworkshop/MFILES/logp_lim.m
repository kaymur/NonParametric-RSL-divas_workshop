function [logp] = logp_lim(y0,f0,stdev,obs,stY,limiting,min_prob,dt)
% finds the log-likelihood of y0, given observed height obs (Z), and posterior
% GP, conditional on all other y (which has mu=f0, and s=stdev)
%
% INPUT
%	y0			sampled point we are evaluating the likelihood of
%   f0          mean of the predicted sea level, given all other data
%	stdev		standard deviation of the predicted sea level with augmented Gaussian measurement error
%   obs         Z = the observed elevation of the data point
%   limiting    indicator of whether marine limiting or terrestrial limiting
% Last updated by Erica Ashe, Wed Jan 24 2017
%%%%%
%interval = [floor(min(y0)-3*max(dY)):.1:ceil(max(y0)+3*max(dY))+20];
eval=y0;
f_0=f0;
st=stdev+stY;
        f_x = @(y) lim_indic(y,obs,limiting,stY);
        log_prob_ln = max(f_x(eval),10e-12);
            f_y = @(y) normpdf(y,f_0,st);
%        log_prob_norm(kk,:)=1;
        log_prob_norm = max(f_y(eval),10e-12);               
        log_probs = log_prob_ln*log_prob_norm;
        %H = @(x)f_x(eval(kk)-x).*f_y(x);
    %         log_probs(kk,:)= integral(H, ...
    %                       min([exp(lnmu),f_0])-3*max([exp(lnstd),st])  ,...
    %                       max([exp(lnmu),f_0])+3*max([exp(lnstd),st]));
            %logp=log(log_probs);
    logp=log(log_probs);
        %logp=log(log_probs);
end