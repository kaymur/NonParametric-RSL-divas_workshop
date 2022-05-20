% look at autocorrelations for thetas and thinned_ys
% calculate other metrics to make sure the model converges 
% can uncomment different parts for different checks

% analyse the y_samples:
    for iitime=1:22:size(thinned_ys,1)
        plotmatrix(thinned_ys(iitime:iitime+5,:)')
        title(['Posterior Dist Samples ' num2str(iitime)]);
        pdfwrite(['Posterior_Dist_Samp' num2str(iitime)]);
    end
    
    %% trace of data
%% look at trace of the orbicella data
c9 = find(cspecies(trainsubz)==9);
for ii=1:c9
    clf; clear title;
%    plot(y_samples(iithet,:));
     plot(thinned_ys(ii,:));
    title(['Trace of Data ' num2str(meantime(trainsubz(ii)))]);
    pdfwrite(['Trace_Data_Samp' num2str(meantime(trainsubz(ii)))]);
end

for ii=1:c9
    clf; clear title;
%    plot(y_samples(iithet,:));
     plot(steps_keep(ii,:));
    title(['Trace Steps ' num2str(meantime(trainsubz(ii)))]);
    pdfwrite(['Trace_Steps' num2str(meantime(trainsubz(ii)))]);
end

% for ii=c9
%     clf; clear title;
% %    plot(y_samples(iithet,:));
%      plot(exp(logp_keep(ii,:)));
%     title(['Trace Probs ' num2str(meantime(ii))]);
%     pdfwrite(['Trace_Probs' num2str(meantime(ii))]);
% end


% %% trace of data
%% look at trace of the orbicella data
c1 = find(cspecies(trainsubz)==1);
for ii=1:c1
    clf; clear title;
     plot(thinned_ys(ii,:));
    title(['Trace of Data Acrop' num2str(meantime(trainsubz(ii)))]);
    pdfwrite(['Trace_Data_AP' num2str(meantime(trainsubz(ii)))]);
end

for ii=1:c1
    clf; clear title;
     plot(steps_keep(ii,:));
    title(['Trace Steps Acrop' num2str(meantime(trainsubz(ii)))]);
    pdfwrite(['Trace_Steps_AP' num2str(meantime(trainsubz(ii)))]);
end

for ii=1:c1
    clf; clear title;
     plot(logp_keep(ii,:));
    title(['Trace Probs Acrop' num2str(meantime(trainsubz(ii)))]);
    pdfwrite(['Trace_Probs_AP' num2str(meantime(trainsubz(ii)))]);
end

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
% 
        clf;

        %% draw histograms of thetas from thinned samples
        plotmatrix(y_accepted(2:7,:)');
        %plotmatrix(thinned_ys(2:7,:)');
         title('Y Posteriors 1');
         pdfwrite('Y_posteriors1');
        plotmatrix(y_accepted(8:13,:)');
         title('Y Posteriors 2');
         pdfwrite('Y_posteriors2');

    sub9=find(cspecies(train_all)==9);
    for iii=1:6:length(sub9)-6
        clf;
        clear title;
        plotmatrix(y_accepted(sub9(iii:iii+5),:)');
         title('Y Posteriors Orbicella 1');
         pdfwrite(['Y_posteriorsOrb1' num2str(iii)]);
    end
    clf;
    for jj=1:7:size(sub9,1)
        clf; 
        %plot(thet_samples(jj,[500:1500 2000:3000 3500:4500]));
        %plot(thet_samples(jj,:));
        if jj+6<size(sub9,1)
            for ii=jj:jj+6
                plot(y_accepted(sub9(ii),:));
                hold on;
            end
            title('Trace Orb Samples');     
            pdfwrite(['Orb_Trace' num2str(jj)]);
        end
    end

clf;
plotmatrix(thets(:,:)');
title('All Hyperparam Posteriors');
pdfwrite('All_Hyper_thinned');

% clf;
% for jj=1:32:1000
% plot(testt,fs_all(jj,:)); hold on;
% end

thets=thinned_thets(:,nburn:end);
for jj=1:length(theta0)
     clf;clear title;
        autocorr(thets(jj,:));%,lags);
        %autocorr(thinned_thets(jj,50:end));%,lags);
        %autocorr(thinned_thets2(jj,:));%,lags);
        %title(['autocorr' num2str(i)]);% ' - ' num2str(Y(trainsubz(i))+meanSL)]);
        %pdfwrite(['autocorr_thinned_more_thets' num2str(jj)]);
       pdfwrite(['ACF_thets' num2str(jj)]);
end

figure;
for jj=1:length(theta0)
    %plot(thet_samples(jj,[500:1500 2000:3000 3500:4500]));
    %plot(thet_samples(jj,:));
    plot(thets(jj,:));
    hold on;
end
title('Trace Thetas Burned');
%title('Trace Thetas Sampled');
legend({'1','2','3','4','5','6','7','8','9'});%,'6'});
pdfwrite(['trace_thets']);

%% draw histograms of thetas
for iithet=1:size(thet_samples,1)
    clf; clear title;
%    histogram(thet_samples(iithet,end-2000:end));
     histogram(thets(iithet,:),20);
%     histogram(thet_samples(iithet,2000:end),20);
%    histogram(thet_samples(iithet,1:226),20);
    title(['Posterior Dist on Theta Samples ' num2str(iithet)]);
    pdfwrite(['Posterior_Dist_Theta_Samp' num2str(iithet)]);
end
%pdfwrite(['trace_thets_samped_thinned']);

st=size(thets,1);
plotmatrix(thets(1:st,:)');
title('Theta Posteriors');
pdfwrite('Thet_posteriors');


% figure;
% for jj=1:length(theta0)
%     plot(steps(:,jj+length(trainsubz)));
%     hold on;
% end
% title('Step Sizes end');
% legend({'Glob Amp','Temp Scale','Regional Amp','Temp Scale','Length Scale','White Noise'});
% pdfwrite(['Step_sizes_end']);
% 

%y_samps=y_samples;

% analyse the ages:
%     for iitime=1:51:500
%         plotmatrix(t_samples(iitime:iitime+5,:)')
%         title(['Posterior Dist on Age Samples ' num2str(iitime)]);
%         pdfwrite(['Posterior_Dist_Age_Samp' num2str(iitime)]);
%     end
%     for iitime=1:51:500
%          clf; clear title;
%          plot(t_samples(iitime,:));
%          title(['Trace of Ages ' num2str(iitime)]);
%          pdfwrite(['Trace_Ages_Samp' num2str(iitime)]);
%     end
% 
% % analyse the y_samples:
%     for iitime=1:51:size(thinned_ys,1)
%         plotmatrix(thinned_ys(iitime:iitime+5,:)')
%         title(['Posterior Dist Samples ' num2str(iitime)]);
%         pdfwrite(['Posterior_Dist_Samp' num2str(iitime)]);
%     end
% 
% %% trace of data
% for iithet=20:55:420
%     clf; clear title;
% %    plot(y_samples(iithet,:));
%      plot(thinned_ys(iithet,:));
%     title(['Trace of Data ' num2str(iithet)]);
%     pdfwrite(['Trace_Data_Samp' num2str(iithet)]);
% end
% 
% %% trace of thetas
% %thet_sams
% %for iithet=1:size(thet_samples,1)
% % for iithet=1:size(thet_samples,1)
% %     clf; clear title;
% %     plot(thet_samples(iithet,:));
% %     %plot(thet_samples(iithet,1:1000));
% %     title(['Trace of Theta ' num2str(iithet)]);
% %     pdfwrite(['Trace_Theta_Samp' num2str(iithet)]);
% % end
% %% trace of thetas after burning 5000
% % for iithet=1:size(thet_samples,1)
% %     clf; clear title;
% %     %plot(thet_samples(iithet,5000:end));
% %     plot(thet_samples(iithet,:));
% %     title(['Trace of Burned Thetas ' num2str(iithet)]);
% %     pdfwrite(['Trace_Burned_Theta_' num2str(iithet)]);
% % end
% 
% % 
% % if size(thinned_thets,2) >80
% %     for iithet=1:size(thinned_thets,1)
% %     clf; clear title;
% % %    hist(thet_samples(iithet,end-2000:end));
% % %     hist(thet_samples(iithet,:));
% % %    hist(thinned_thets(iithet,:));
% %         hist(thinned_thets(iithet,81:end));
% %         title(['Last 100 Post Dist on Theta ' num2str(iithet)]);
% %         pdfwrite(['Last_100_Posterior_Dist_Theta_' num2str(iithet)]);
% %     end
% % elseif size(thinned_thets,2) >5
% % %    for iithet=1:size(thinned_thets,1)
% % %    clf; clear title;
% %         for iithet=1:size(thinned_thets,1)
% %             clf; clear title;
% % %            figure;
% %         %    histogram(thet_samples(iithet,end-2000:end));
% %         %     histogram(thet_samples(iithet,:),20);
% %             if iithet==1
% %                 edges=[0:5e3:100e3];
% %             elseif iithet==2
% %                 edges=[0:1e3:13e3];
% %             elseif iithet==3
% %                 edges=[500:250:4500];
% %             elseif iithet==4
% %                 edges=[1000:500:11000];
% %             elseif iithet==5
% %                 edges=[1:0.25:7];
% %             elseif iithet==6
% %                 edges=[0:20:400];
% %             end
% %             histogram(thinned_thets(iithet,:),edges);
% %             title(['Posterior Distribution on Theta ' num2str(iithet)]);
% %             pdfwrite(['Post_Dist_Theta_Thinned' num2str(iithet)]);
% %         end
% % %    end
% % end
% % %         histogram(thinned_thets(iithet,5:end),edges);
% % %         histogram(thinned_thets(iithet,5:end),edges);
% % %         title(['Last X Post Dist on Theta ' num2str(iithet)]);
% % %         pdfwrite(['Last_Posterior_Dist_Theta_' num2str(iithet)]);
% % %     end
% % % end
% % %% draw histograms of thetas from total samples
% % % for iithet=1:size(thinned_thets,1)
% % %     clf; clear title;
% % % %    hist(thet_samples(iithet,end-2000:end));
% % %      hist(thet_samples(iithet,:));
% % %     %hist(thinned_thets(iithet,:));
% % %     title(['Posterior Dist on All Theta ' num2str(iithet)]);
% % %     pdfwrite(['Post_Dist_All_Theta_' num2str(iithet)]);
% % % end
% % 
% % % % plot last 5000
% % % for jj=1:size(sels,1)
% % %     plot(sels(jj,5000:10000));
% % %     pdfwrite(['sels_' num2str(jj)]);
% % % end
% % % for jj=1:size(props,2)
% % %     plot(props(5000:10000,jj));
% % %     pdfwrite(['props_' num2str(jj)]);
% % % end
% % % for jj=1:size(logps,1)
% % %     plot(logps(jj,5000:10000));
% % %     pdfwrite(['logps_' num2str(jj)]);
% % % end
% % 
% % % plot all 10k
% % % sels are the selectors: random numbers to compare to acceptance ratio
% % % for jj=1:size(sels,1)
% % %     plot(sels(jj,1:Nsamps));
% % %     pdfwrite(['sels_all' num2str(jj)]);
% % % end
% % % props are the proposed new value
% % % for jj=1:size(props,2)
% % %     plot(props(1:Nsamps,jj));
% % %     title(['props ' num2str(jj)]);
% % %     pdfwrite(['props_all' num2str(jj)]);
% % % end
% % % for jj=1:size(props,1)
% % %     plot(props(jj,:));
% % %     title(['props ' num2str(jj)]);
% % %     pdfwrite(['props_all' num2str(jj)]);
% % % end
% % % % logps are the likelihood value of the new proposed theta
% % % for jj=1:size(logps,1)
% % %     %plot(logps(jj,1:500));
% % %     plot(logps(jj,:));
% % %     title(['logps ' num2str(jj)]);
% % %     pdfwrite(['logps_all' num2str(jj)]);
% % % end
% % %% autocorrelation for the data:
% % % lags=30;
% % % for jj=1:length(theta0)
% % %      clf;
% % %         autocorr(t_samples(jj,:),lags);
% % %         %autocorr(thet_samples(jj,end-2000:end),lags);
% % %         title(['autocorr ages' num2str(jj)]);% ' - ' num2str(Y(trainsubz(i))+meanSL)]);
% % %         pdfwrite(['autocorr_ages_end' num2str(jj)]);
% % % end
% % 
% % % lags=150;
% % % for jj=1:length(theta0)
% % % %     for ii=1:length(thet_samples(1,1,:))
% % %          clf;
% % %         %autocorr(thinned_thets0(jj,800:1000),lags);
% % % %        autocorr(thet_samples(jj,3000:end),lags);
% % % %        autocorr(thet_samples(jj,end-500:end,ii),lags);
% % %         autocorr(thet_samples(jj,:),lags);
% % %         %autocorr(thet_samples(jj,end-2000:end),lags);
% % %         %title(['autocorr' num2str(i)]);% ' - ' num2str(Y(trainsubz(i))+meanSL)]);
% % % %        pdfwrite(['autocorr_thets_end_' num2str(jj)]);
% % %         pdfwrite(['autocorr_thets_all_' num2str(jj) '_' num2str(ii)]);
% % % %     end
% % % end
% % % lags=30;
% % % for jj=1:length(theta0)
% % %      clf;
% % %         autocorr(thinned_thets(jj,:));%,lags);
% % %         %title(['autocorr' num2str(i)]);% ' - ' num2str(Y(trainsubz(i))+meanSL)]);
% % %         pdfwrite(['autocorr_thinned_thets' num2str(jj)]);
% % % end
% % 
