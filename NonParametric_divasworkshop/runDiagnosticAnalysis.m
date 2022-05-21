% visually inspect autocorrelation plots (both for thetas and the data) to determine the correct thinning factor
% ("nthin" and "thetachange" in the runSetUp.m script)

% test for convergence by visually inspecting the trace plots to determine
% the burn-in factor ("nburn" and "stepchange" in the "runSetUp.m" script)

% visually inspect the matrix of thetas to determine if they look reasonable (e.g., pushing up against a
% lower or upper bound)


%% analyse the (already thinned by given nthin) posterior samples of RSLs (ys)
    % as a matrix to plot several at a time
for ii=1:5:10     %size(thinned_ys,1)-7     
    plotmatrix(thinned_ys(ii:ii+5,:)')
    %plotmatrix(thinned_ys(iitime:iitime+5,:)')
    display(ii+5)
    subtitle(['Dist Samps-' num2str(ii)]);
    pdfwrite(['Posterior_Dist_Samp' num2str(ii)]);
end
    % one at a time
for ii=1:5          %size(thinned_ys,1)  % uncomment to go through every y or change the ii numbers
    plotmatrix(thinned_ys(ii:ii+5,:)')
    %plotmatrix(thinned_ys(iitime:iitime+5,:)')
    display(ii+5)
end
%     subtitle(['Dist Samps-' num2str(ii)]);
%     pdfwrite(['Posterior_Dist_Samp' num2str(ii)]);


%%%%%%%%%%%%%%%%%%%%%%%%%

%% trace of data
    % several on the same plot
clf;
for ii=1:5          % change ii=6:10 to do the next set
    plot(thinned_ys(ii,:));
    hold on;
end
% title(['Trace of Data ' num2str(meantime(trainsubz(ii)))]);
% pdfwrite(['Trace_Data_Samp' num2str(meantime(trainsubz(ii)))]);


%%%%%%%%%%%%%%%%%%%%%%%%%

%% look at the step sizes
clf; clear title;
for ii=1:5
    plot(steps_keep(ii,:));
    hold on;
end
% title(['Trace Steps ' num2str(meantime(trainsubz(ii)))]);
% pdfwrite(['Trace_Steps' num2str(meantime(trainsubz(ii)))]);


%%%%%%%%%%%%%%%%%%%%%%%%%

%% look at the log-likelihoods
clf; clear title;
for ii=1:5
    plot(exp(logp_keep(ii,:)));
    hold on;
end
% title(['Trace Likelihoods ' num2str(meantime(trainsubz(ii)))]);
% pdfwrite(['Trace_Likes_' num2str(meantime(trainsubz(ii)))]);

% %% trace of data

%%%%%%%%%%%%%%%%%%%%%%%%%

Nsamps=Nsamples;
lags=100;
if exist('all_ys')>0
    for jj=1:size(trainsubz)
        clf; clear title;
        autocorr(all_ys(jj,5000:100:end),lags);
        title(['Autocorr Data' num2str(cspecies(jj))]);
        pdfwrite(['ACF_all_ys_burn_5000' num2str(jj) '_type_' num2str(cspecies(jj))]);
    end
end
