function [logp] = logp_lim(y0,f0,stdev,obs,stY,limiting,min_prob,tri_max)
% finds the log-likelihood of y0, given observed height obs (Z), and posterior
% GP, conditional on all other y (which has mu=f0, and s=stdev)
%
% INPUT
%	y0			sampled point we are evaluating the likelihood of
%   f0          mean of the predicted sea level, given all other data
%	stdev		standard deviation of the predicted sea level with augmented Gaussian measurement error
%   obs         Z = the observed elevation of the data point
%   limiting    indicator of whether marine limiting or terrestrial limiting
%   min_prob    the minimum probability for the normal or the limiting 
% Last updated by Erica Ashe, May 19 2022
%%%%%
defval('min_prob', 10e-12);
defval('tri_max', 50e3);
defval('stY', 0);
defval('stdev', 0);

% eval=y0;
% f_0=f0;
st=stdev+stY;
        f_x = @(y) lim_indic(y,obs,limiting,st,tri_max);
        %log_prob_ln = max(f_x(eval),min_prob);
        log_prob_ln = max(f_x(y0),min_prob);
        f_y = @(y) normpdf(y,f0,st);
        %f_y = @(y) normpdf(y,f_0,st);
        log_prob_norm = max(f_y(y0),min_prob);               
        log_probs = log_prob_ln*log_prob_norm;
    logp=log(log_probs);
end