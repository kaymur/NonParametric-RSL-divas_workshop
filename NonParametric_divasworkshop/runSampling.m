runSetUpSampling;
tic;

for nn=xx:Nsamples
    mm=mod(nn,Nthin);
    if mm~=0
        y_samples(:,mm+1,:) = y_samples(:,mm,:);
        thet_samples(:,mm+1,:) = thet_samples(:,mm,:);
    end    
    for pp=1:Ntemps
        wYs = y_samples(:,mm+1,pp);
        wY(trainsubz)=wYs;
        proposal = wYs;     % jump from last Y to proposed next Y
        for ii = randperm(length(wYs))
            proposal(ii)= wYs(ii) + randn*steps(ii);
            trainsub_null=find(limiting<10);       % all whether sampled or not    
            trainsub_temp = trainsub_null(trainsub_null~=trainsubz(ii));
            sub_temp = find(trainsubz~=trainsubz(ii));
            try
                [wfs(ii),wVs(ii)] = GaussianProcessRegression(meantime(trainsub_temp),wY(trainsub_temp),meantime(trainsubz(ii)),...
                        modelspec.traincv(t1(trainsub_temp),t1(trainsub_temp),dt1t1(trainsub_temp,trainsub_temp),wthet,Ycv(trainsub_temp,trainsub_temp),dy1y1(trainsub_temp,trainsub_temp),fp1fp1(trainsub_temp,trainsub_temp)),...
                        modelspec.cvfunc(t1(trainsub_temp),t1(trainsubz(ii)),dt1t1(trainsub_temp,trainsubz(ii)),wthet,dy1y1(trainsub_temp,trainsubz(ii)),fp1fp1(trainsub_temp,trainsubz(ii))),...
                        modelspec.cvfunc(t1(trainsubz(ii)),t1(trainsubz(ii)),dt1t1(trainsubz(ii),trainsubz(ii)),wthet,dy1y1(trainsubz(ii),trainsubz(ii)),fp1fp1(trainsubz(ii),trainsubz(ii))),[]);
                stds(ii)=sqrt(wVs(ii) + dYs(trainsubz(ii))^2);
            end
            if cspecies(trainsubz(ii))==0
                if limiting(trainsubz(ii))==-1||limiting(trainsubz(ii))==1
                    %logpproposed(ii)=logp_limiting(proposal(ii),wfs(ii),stds(ii),Z0(ii),dYs(trainsubz(ii)),limiting(trainsubz(ii)),1,dt(trainsubz(ii)));
                    %logpproposed(ii)=logp_lim(proposal(ii),wfs(ii),stds(ii),Z0(ii),dYs(trainsubz(ii)),limiting(trainsubz(ii)),1,dt(trainsubz(ii)));
                    %logpproposed(ii)=logp_lim(proposal(ii),wfs(ii),stds(ii),Z0(ii),dYs(trainsubz(ii)),limiting(trainsubz(ii)),,dt(trainsubz(ii)));
                    logpproposed(ii)=logp_lim(proposal(ii),wfs(ii),[],Z0(ii),[],limiting(trainsubz(ii)),10e-12,tri_max);
                else
                    logpproposed(ii)=logp_norm(proposal(ii),Z(ii),dYs(ii),Z0(ii),0,dYs(ii));
                end
%             elseif cspecies(trainsubz(ii))>11
%                     logpproposed(ii)=logp_norm(proposal(ii),Z(ii),dYs(ii),Z0(ii),0,dYs(ii));
%             else   
%                 if cspecies(trainsubz(ii))==9
%                      logpproposed(ii)= log(normpdf(proposal(ii),wfs(ii),stds(ii)))+...
%                         log(pdf(distKern{cspecies(trainsubz(ii))},(proposal(ii)-Z0(ii))/1000));       
% %                    logpproposed(ii)= log(pc*normpdf(proposal(ii),wfs(ii),stds(ii)))+...
% %                       log((1-pc)*pdf(distKern{cspecies(trainsubz(ii))},(proposal(ii)-Z0(ii))/1000));       
%                     %disp(fprintf('Num: %0.0f , Prop: %0.0f , Y:%0.0f , Diff: %0.0f , Prob: %0.0f , probs: %0.0f, f: %0.0f, st: %0.0f',[ii proposal(ii) Z0(ii) proposal(ii)-Z0(ii) log(pdf(distKern{cspecies(trainsubz(ii))},(proposal(ii)-Z0(ii))/1000)) log(normpdf(proposal(ii),wfs(ii),stds(ii))) wfs(ii) stds(ii)]));
            else
                logpproposed(ii)= log(normpdf(proposal(ii),wfs(ii),stds(ii))) + ...
                    log(pdf(distKern{cspecies(trainsubz(ii))},(proposal(ii)-Z0(ii))/1000));                                     
%                    logpproposed(ii)= log(pc*normpdf(proposal(ii),wfs(ii),stds(ii)))+...
%                       log((1-pc)*pdf(distKern{cspecies(trainsubz(ii))},(proposal(ii)-Z0(ii))/1000));                                           
            end
        end
    end
    selector=log(rand);
    cutoff = (logpproposed(ii)-wlogpold(ii))*temps(pp);
    boundsok = (sum(proposal(ii)<lb)==0).*(sum(proposal(ii)>ub)==0);
    if ~boundsok
        %disp(sprintf('%0.0f - %0.0f -- Round %0.0f - %0.0f -- out of bounds',[jj pp nn ii]));
    elseif (selector < cutoff)
        %disp(fprintf('%0.0f - %0.0f -- Round %0.0f - %0.0f -- accepted %0.2f -- %0.3f (%0.3f) -- %0.3f < %0.3f',[jj pp nn ii proposal(ii) logpproposed(ii) wlogpold(ii) selector cutoff]));
        wlogpold(ii) = logpproposed(ii);
        wYs(ii) = proposal(ii);
        wY(trainsubz)=wYs;
        acc_count(ii,pp) = acc_count(ii,pp)+1;
        wacc_count(ii,pp)=wacc_count(ii,pp)+1;
    else
%             disp(fprintf('%0.0f -- Round %0.0f, Prop: %0.2f, Obs: %0.2f -- rejected %0.2f-- %0.3f > %0.3f',[ii nn proposal(ii) Z0(ii) selector cutoff]));
    end                % change step size to get ideal 
            %   if none accepted, increase step-size 
            %   if accept_count too low, decrease step size
            %   if accept_count is too high, increase step-size
    y_samples(:,mm+1,pp) = wYs;
    for ii = randperm(length(wthet))
        proposalT=wthet;
        proposalT(ii) = proposalT(ii) + randn*steps(ii+length(wYs));  
        boundsok = (sum(proposalT<lb_thet')==0).*(sum(proposalT>ub_thet')==0);
        if ~boundsok
    %                %disp(sprintf('%0.0f -- Round %0.0f - %0.0f -- out of bounds',[pp nn ii]));
        else
            % basing the hyperparameters on all of the data, not just that being sampled
            if sampnorm==1
                 wlogpoldT=logprob(wYs,@(theta)modelspec.traincv(t1(train_all),t1(train_all),dt1t1(train_all,train_all),theta,Ycv(train_all,train_all),dy1y1(train_all,train_all),fp1fp1(train_all,train_all)),wthet);
                 logpproposedT=logprob(wYs,@(theta)modelspec.traincv(t1(train_all),t1(train_all),dt1t1(train_all,train_all),theta,Ycv(train_all,train_all),dy1y1(train_all,train_all),fp1fp1(train_all,train_all)),proposalT);
            else
                 wlogpoldT=logprob(wY(train_all),@(theta)modelspec.traincv(t1(train_all),t1(train_all),dt1t1(train_all,train_all),theta,Ycv(train_all,train_all),dy1y1(train_all,train_all),fp1fp1(train_all,train_all)),wthet);
                 logpproposedT=logprob(wY(train_all),@(theta)modelspec.traincv(t1(train_all),t1(train_all),dt1t1(train_all,train_all),theta,Ycv(train_all,train_all),dy1y1(train_all,train_all),fp1fp1(train_all,train_all)),proposalT);
            end
            %                 [~,~,logpproposed(ii)] = GaussianProcessRegression_ea(meantime(trainsub),wYs,testt,...
            %                         modelspec.traincv(t1,t1,dt1t1,proposal,Ycv(trainsub,trainsub),dy1y1,fp1fp1),...
            %                         modelspec.cvfunc(t1,t2,dt1t2,proposal,dy1y2,fp1fp2)',...
            %                         modelspec.cvfunc(t2,t2,dt2t2,proposal,dy2y2,fp2fp2),[]);
            selector=log(rand);
            cutoff = (logpproposedT-wlogpoldT);%*temps(pp);
            if (selector < cutoff)
                %disp(sprintf('%0.0f -- Round %0.0f - %0.0f -- accepted %0.2f -- %0.3f (%0.3f) -- %0.3f < %0.3f',[pp nn ii proposalT(ii) logpproposedT wlogpoldT exp(selector) exp(cutoff)]));
                wlogpoldT = logpproposedT;
                wthet(ii) = proposalT(ii);
                acc_count(ii+length(wYs),pp) = acc_count(ii+length(wYs),pp)+1;
                wacc_count(ii+length(wYs),pp)=wacc_count(ii+length(wYs),pp)+1;
            else
                %disp(sprintf('%0.0f -- Round %0.0f - %0.0f -- rejected %0.2f-- %0.3f > %0.3f',[pp nn ii proposalT(ii) exp(selector) exp(cutoff)]));
            end
        end  
    end
        
    acc_keep(:,nn)=wacc_count;
    prop_keep(:,nn)=proposal;
    logp_keep(:,nn)=logpproposed';
    steps_keep(:,nn)=steps;
    thet_samples(:,mm+1,pp)=wthet;
     
    if mod(nn,step_change)==0
        %   disp([sprintf('%0.0f -- ',pp) sprintf('%0.2f ',sum(acc_count))]);
        %   disp(sprintf('%0.3f -- Avg Accept, nn = %0.0f',[mean(acc_count)/nn, nn]));
        sub_big=intersect(find(wacc_count/step_change>.45),find(steps<maxstep));
        sub_small=intersect(find(wacc_count/step_change <.25),find(wacc_count>0));
        sub_none=intersect(find(wacc_count<1),find(steps<maxstep));
        if size(sub_none)>0
            steps(sub_none)= steps(sub_none)*1.8;                
        end
        if size(sub_big)>0
            steps(sub_big)= steps(sub_big)*1.15;
        end
        if size(sub_small)>0
            steps(sub_small)= steps(sub_small)/1.42;
        end
        sub_none=[]; sub_big=[]; sub_small=[];
        wacc_count=0*acc_count;
    end
    
    if mod(nn,Nthin)==0
        thinned_ys(:,nn/Nthin) = y_samples(:,end,pp);
        y_samples(:,1,pp)=y_samples(:,end,pp);
        thinned_thets(:,nn/Nthin) = thet_samples(:,end,pp);
        thet_samples(:,1,pp)=thet_samples(:,end,pp);
        timepassed=toc;
%         if length(savefile)>0
%             if exist('Y_fin')
%                 save(savefile,'wthet','thinned_ys','thinned_thets','y_samples','acc_count','logp','logpnorm','steps','thet_samples','jj','nn','dataset','modelspec','theta0',...
%                     'Nsamples','Nthin','steps','Nburnin','step_change','thet_change','Seed','Y_fin','dY','limiting','cspecies');
%             else
%                 save(savefile,'wthet','thinned_ys','thinned_thets','y_samples','acc_count','logp','logpnorm','steps','thet_samples','jj','nn','dataset','modelspec','theta0',...
%                     'Nsamples','Nthin','steps','Nburnin','step_change','thet_change','Seed','dY','limiting','cspecies');
%             end
%         end
    end
    if mod(nn,200)==0
        timepassed=toc;
        disp(sprintf('%0.0f seconds passed -- %0.0f estimated seconds remaining',[timepassed timepassed/(nn-1)*(Nsamples-nn)]));        
    end
end

step_size(trainsubz)=steps(1:length(trainsubz));
%save;