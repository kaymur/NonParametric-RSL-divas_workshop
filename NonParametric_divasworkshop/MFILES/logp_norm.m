function [logp] = logp_norm(y0,f0,stdev,obs,lnmu,lnstd)
% finds the log-likelihood of y0, the proposed sample, given observed height obs (Z), and posterior
% GP, conditional on all other y (which has mu=f0, and s=stdev)
%
% INPUT
%
%	y0			sampled point we are evaluating the likelihood of
%   f0          mean of the predicted sea level, given all other data
%	stdev		standard deviation of the predicted sea level with augmented Gaussian measurement error
%   obs         Z = the observed data point
%	lnmu		normal mean of depth distribution
%	lnstd		normal standard deviation of depth distribution
% Last updated by Erica Ashe, Wed Jan 19 2017
%%%%%
%interval = [floor(min(y0)-3*max(dY)):.1:ceil(max(y0)+3*max(dY))+20];
eval=y0/1000-obs/1000; 
f_0=f0/1000;
st=stdev/1000;
    for kk=1:length(eval)
         f_x = @(y) normpdf(y,lnmu,lnstd);
        log_prob(kk,:) = f_x(eval(kk));
%        log_prob_ln(kk,:) = 1;
        f_y = @(y) normpdf(y,f_0,st);
%        log_prob_norm(kk,:)=1;
        log_prob_norm(kk,:) = f_y(y0(kk)/1000);
        log_probs(kk,:) = log_prob(kk,:)*log_prob_norm(kk,:);
%         [log_prob_ln log_prob_norm]
        %H = @(x)f_x(eval(kk)-x).*f_y(x);
    %         log_probs(kk,:)= integral(H, ...
    %                       min([exp(lnmu),f_0])-3*max([exp(lnstd),st])  ,...
    %                       max([exp(lnmu),f_0])+3*max([exp(lnstd),st]));
            %logp=log(log_probs);
    end
    logp=log(log_probs);
        %logp=log(log_probs);
end