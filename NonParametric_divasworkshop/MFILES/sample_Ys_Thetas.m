
tic;

maxstep=max(8*dY(trainsubz),1000);
minstep=100;
for iii=1:length(modelspec.ub)
    maxstep=[maxstep; modelspec.ub(iii)-modelspec.lb(iii)];
end
if ~exist('thinned_ys')
    new=1;
    xx=1;
    wY=Y;
elseif ~exist('wthet')
    new=0;
    xx=size(thinned_ys,2)-1;
    theta0=wthet';
else
    new=0; xx=nn-1;
end
if exist('datasets')
    dataset=datasets{1};
end
if ~exist('wY')
    wY=Y;
end
    modelspec=modelspec(1);
    Nsamples=nsamps;
    Nburnin=nburn;
    lb=-150000; ub=20000;
    Ntemps=1;
    defval('deltaT',.01)
    defval('step_change',60);
    defval('thet_change', 100);
    defval('jj',1);
    temps = 1./[1+deltaT*(0:(Ntemps-1))];
if exist('Seed','var') 
    savefile=['intermed_samps' num2str(Seed)];
else
    savefile=['intermed_samps'];
end
%
% INPUT
%
%	dataset			the whole dataset
%	modelspec		defines the covariance function
%	theta0          fixed theta (hyperparameter) value at which to find likelihoods
%	Nsamples        number of samples to return
%   Nthin           number to thin samples by
%   coral_params    matrix of parameter values for log-normal distributions
%   step_size		size by which to jump from one y-value to another
%   thet_change     number of samples before re-optimizing theta
%   step_change     number of samples before changing the step size to get ideal acceptance
    
% OUTPUT
%
%   thinned_y       y-values of non-normally distributed data to be used in
%                       GP regression, afger being thinned by Nthin
%   acc_count       the number of samples that were accepted
%   step_size       the step_size for each sample at the end
%   samp_norm       indicates whether to sample the normally distributed data

% Last updated by Erica Ashe, Fri Nov 8 2019

%%%%%
    %disp(['Temperatures: ' sprintf('%0.2f ',temps)]);

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));
dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);
cspecies=dataset.cspecies;
datid = dataset.datid;
time1=dataset.time1;
time2=dataset.time2;
dt=dataset.dt;
limiting=dataset.limiting;
%indic=dataset.indic;
lat=dataset.lat;
long=dataset.long;
Y=dataset.Y;
Ycv=dataset.Ycv;
dYs=dataset.dY;
siteid=dataset.siteid;
sitenames=dataset.sitenames;
meantime=dataset.meantime;
sitecoords=dataset.sitecoords;
sitelen=dataset.sitelen;
if isfield(dataset,'obsGISfp')
    obsGISfp=dataset.obsGISfp;
else
    obsGISfp=ones(size(lat));
end

    defval('Nsamples',1000);
    defval('Nthin',10);
    defval('Nburnin',10);
    defval('savefile','');
    defval('plotfile','Sample_file');

    sub_all = intersect(find(limiting~=1),find(limiting~=-1));
    
    if exist('dat_type','var')
        if dat_type==3 
            trainsubz=find(limiting~=20);
        else
            trainsubz = find(limiting~=0&limiting~=20);
        end
    else
        trainsubz = find(limiting~=0&limiting~=20);
    end
    train_all = find(limiting~=20);
    Z = wY(train_all);
    if sampnorm
        trainsubz = find(limiting~=20);
        Z0 = Y(train_all);
    else
        Z0=wY(trainsubz);
    end
    %% Z is all of the data that is being used in the model, not just sampling
    ts0=intersect(find(limiting(train_all)==0),find(cspecies(train_all)~=0));
    trainsub0=find(limiting==0);
    sub02=find(limiting==0|limiting==2);
    testX = [lat(train_all) long(train_all) meantime(train_all)];
    lb_thet=modelspec(1).lb;
    ub_thet=modelspec(1).ub;
if new ==1
    stepT=(ub_thet-lb_thet)/120;
    stepT(4)=20;
    steps=[step_size(trainsubz); stepT'];
end
    mspec.cvfunc = @(x1,x2,thet) modelspec.cvfunc(x1,x2,dYears(x1,x2),thet,dy1y1',fp1fp1');
    mspec.dcvfunc = @(x1,x2,thet) modelspec.dcvfunc(x1,x2,dYears(x1,x2),thet,dy1y1',fp1fp1');
    mspec.ddcvfunc = @(x1,x2,thet) modelspec.ddcvfunc(x1,x2,dYears(x1,x2),thet,dy1y1',fp1fp1');
    Ys0=wY(trainsubz);

X1 = [lat long meantime];
dt1t1 = dYears(X1(:,3),X1(:,3));
dt1t2 = dYears(X1(:,3),testX(:,3));
dt2t2 = dYears(testX(:,3),testX(:,3));

dy1y1 = dDist(X1(:,1:2),X1(:,1:2));
dy1y2 = dDist(X1(:,1:2),testX(:,1:2));
dy2y2 = dDist(testX(:,1:2),testX(:,1:2));

testfp = obsGISfp(trainsubz);
fp1fp1=bsxfun(@times,obsGISfp-1,obsGISfp'-1)';
fp1fp2=bsxfun(@times,obsGISfp-1,testfp'-1)';
fp2fp2=bsxfun(@times,testfp-1,testfp'-1);

t1=X1(:,3);
t2=testX(:,3);
                    
%need logp values for each y being sampled

y_samples = repmat(Ys0,1,Nthin,Ntemps);
thet_samples= repmat(theta0',1,Nthin,Ntemps);
logp=[];
%fprintf(
logp0=-50*ones(length(trainsubz),1);
logp(1,:) = logp0*800;

if ~exist('acc_count')
    acc_count=zeros(size(trainsubz,1)+length(theta0),Ntemps);
    wacc_count=acc_count;
else
    wacc_count=0.*acc_count;
end

wlogpold = logp;
logpproposed = logp;
wlogpoldT = -4000;
logpproposedT = -4000;
wlogpoldT = logp(1);
logpproposedT = logp(1);
logpnorm=[];
wthet = theta0';
thet_samples= repmat(wthet,1,Nthin,Ntemps);

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
                    logpproposed(ii)=logp_limx(proposal(ii),wfs(ii),stds(ii),Z0(ii),dYs(trainsubz(ii)),limiting(trainsubz(ii)),1,dt(trainsubz(ii)));
                else
                    logpproposed(ii)=logp_norm(proposal(ii),Z(ii),dYs(ii),Z0(ii),0,dYs(ii));
                end
            elseif cspecies(trainsubz(ii))>11
                    logpproposed(ii)=logp_norm(proposal(ii),Z(ii),dYs(ii),Z0(ii),0,dYs(ii));
            else   
                if cspecies(trainsubz(ii))==9
                     logpproposed(ii)= log(normpdf(proposal(ii),wfs(ii),stds(ii)))+...
                        log(pdf(distKern{cspecies(trainsubz(ii))},(proposal(ii)-Z0(ii))/1000));       
%                    logpproposed(ii)= log(pc*normpdf(proposal(ii),wfs(ii),stds(ii)))+...
%                       log((1-pc)*pdf(distKern{cspecies(trainsubz(ii))},(proposal(ii)-Z0(ii))/1000));       
                    %disp(fprintf('Num: %0.0f , Prop: %0.0f , Y:%0.0f , Diff: %0.0f , Prob: %0.0f , probs: %0.0f, f: %0.0f, st: %0.0f',[ii proposal(ii) Z0(ii) proposal(ii)-Z0(ii) log(pdf(distKern{cspecies(trainsubz(ii))},(proposal(ii)-Z0(ii))/1000)) log(normpdf(proposal(ii),wfs(ii),stds(ii))) wfs(ii) stds(ii)]));
                else
                     logpproposed(ii)= log(normpdf(proposal(ii),wfs(ii),stds(ii)))+...
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
        end
        
    end
    
    y_samples(:,mm+1,pp) = wYs;
    %all_ys(:,nn) = wYs;

            % change step size to get ideal 
            %   if none accepted, increase step-size 
            %   if accept_count too low, decrease step size
            %   if accept_count is too high, increase step-size
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