% function [thetas,logps,acceptedcount,acceptmatrix,fs,sds,V2s]=SLNIGPSamp_1207(dataset,traincvaug,modelspec,thet0,lb,ub,spacex,prop_cov,Nsamples,Nthin,Nburnin,savefile,Sampsteps,subnotfixed,testt,meanSL)
if ~exist('nnn', 'var')
    nnn=1;
    dt0=datasets{1}.dt;
    dt=dt0;
    tic
    Ntemps=1;
    deltaT=.03;
    Nthin=40;
    Nburnin=100;
    %Nsamples=10000; Nsamps=Nsamples;
    Nsamples=nsamps;
    temps = 1./[1+deltaT*(0:(Ntemps-1))];
    lb_thet=modelspec(1).lb;
    ub_thet=modelspec(1).ub;
    %step_thet=steps(1:6);
    if modno==6&&spatial==0 || synthflag ==1
        steps=(ub_thet-lb_thet)/800;
        step_thet=steps;
    elseif modno~=7
        step_thet=(ub_thet-lb_thet)/800;
        %step_thet=[1700 500 400 450 .08 400 450 .08  30];
        steps=step_thet;
    else
        step_thet=[1700 500 400];
        steps=step_thet;
    end
    %step_thet=(ub_thet-lb_thet)/800;
    lbY=Y(trainsub)-100e3;
    ubY=Y(trainsub)+100e3;
    % ub=[ub_thet ubY'];
    % lb=[lb_thet lbY'];
    ub=[ub_thet];
    lb=[lb_thet];
    step_change=100;
    savefile=['intermed_samps' num2str(Seed) '.mat'];
    plotfile='plots';
    counts=[];

        obsGISfp=ones(size(lat));
    %         testreg=testloc.reg;
    %         testsites=testloc.sites;
        clear testts;
        
        %testloc.X(:,1:2);
        %testloc.X(:,3);

    dYs=dY;
    X1 = [lat long meantime];
    %trainsub=find(limiting>=-20);
    trainsub=find(limiting==0);
    trainsub0=trainsub;
    % trainsub=find(limiting==0);
    % trainsub0=trainsub;

    dt1t1 = dYears(X1(trainsub,3),X1(trainsub,3));
    %dt1t2 = dYears(X1(trainsub,3),testts(:,3));
    %dt2t2 = dYears(testts(:,3),testts(:,3));

    dy1y1 = dDist(X1(trainsub,1:2),X1(trainsub,1:2));
    %dy1y2 = dDist(X1(trainsub,1:2),testts(:,1:2));
    %dy2y2 = dDist(testts(:,1:2),testts(:,1:2));

%     testfp = ones(size(testts(:,3)));
     fp1fp1=bsxfun(@times,obsGISfp(trainsub)-1,obsGISfp(trainsub)'-1)';
%     fp1fp2=bsxfun(@times,obsGISfp(trainsub)-1,testfp'-1)';
%     fp2fp2=bsxfun(@times,testfp-1,testfp'-1);

    t1=X1(trainsub,3);
    wt1=t1;
    limY=limiting(trainsub);
    %t2=testts(:,3);
    %theta0=(lb_thet+ub_thet)./2;
    %theta0=[13600 4750 2730 3224 6.7 338 ];
    %[Ys0,~,logp0] = GaussianProcessRegression(meantime(trainsub0),Y(trainsub0),testt,modelspec.traincv(t1,t1,dt1t1,theta0,0*Ycv(trainsub0,trainsub0),dy1y1,fp1fp1),...
    %                       modelspec.cvfunc(t1,t2,dt1t2,theta0,dy1y2,fp1fp2)',...
    %                       modelspec.cvfunc(t2,t2,dt2t2,theta0,dy2y2,fp2fp2),[]);                     
    logp0=-4.3e4;
    logp = repmat(logp0,[1 1 Ntemps]);
    wlogpold = logp;
    logpproposed = logp;
    wthet = theta0;
    proposal=wthet;
    logps=repmat(logp,1,nsamps);

    % acc_count_t=zeros(size(theta0,2),Ntemps);
    % wacc_count_t=acc_count_t;
    % wlogpold = logp;
    % logpproposed = logp;
    % wthet = theta0;
    % proposal=wthet;
    % logps=repmat(logp,1,Nsamps);
    %logps=repmat(logp,size(theta0,2),Nsamps);
    trainsub_null=find(limY==0);           
    trainsubz = find(limiting==0);
    Z = Y(trainsubz);
    %    testX = [lat(trainsubz) long(trainsubz) meantime(trainsubz)];
X1 = [lat long meantime];
testX = [mean(lat)*ones(length(testt),1) mean(long)*ones(length(testt),1) testt'];
testXs = [X1; testX];
testfp = [ones(length(testt),1); obsGISfp(:)];
                mspec.cvfunc = @(x1,x2,thet) modelspec.cvfunc(x1,x2,dYears(x1,x2),thet,dy1y1',fp1fp1');
                mspec.dcvfunc = @(x1,x2,thet) modelspec.dcvfunc(x1,x2,dYears(x1,x2),thet,dy1y1',fp1fp1');
                mspec.ddcvfunc = @(x1,x2,thet) modelspec.ddcvfunc(x1,x2,dYears(x1,x2),thet,dy1y1',fp1fp1');

    %% using the step size from last run
    %steps=[step_thet (ubY'-lbY')/500];
    %logpY = repmat(logp0*10,[1 Ntemps]);
    logpY = repmat(logp0*10,[length(Z) 1 Ntemps]);
    y_samples = repmat(Z,1,1,Ntemps);
    t_samples = repmat(meantime(trainsub),1,1,Ntemps);
    thet_samples=repmat(theta0',1,1,Ntemps);
    wlogpY = logpY;
    logpropY = logpY;

    % acc_countY=zeros(size(Z,1),Ntemps);
    % wacc_countY=acc_countY;
    acc_count_Y=zeros(length(Z),Ntemps);
    wacc_count_Y=acc_count_Y;

    % acc_countThet=zeros(length(theta0),Ntemps);
    % wacc_countThet=acc_countThet;
    acc_count_thet=zeros(length(theta0),Ntemps);
    wacc_count_thet=acc_count_thet;

    % wlogpoldY = logpY;
    % logpropY = logpY;
    %propY=Y(trainsubz);
    propY=Z;%zeros(length(Z),1);
    wYs = propY;
                %propsY=zeros(length(Z),Nsamps);
                %selsY=zeros(length(Z),Nsamps);

    %acc_count=[acc_count_t; acc_countY];

    testcv = @(thetas) modelspec.cvfunc(t1,t1,dt1t1,thetas,dy1y1,fp1fp1);

            step_sizes=zeros(round(nsamps/step_change), length(theta0));%+length(Z));
                props=zeros(nsamps,length(theta0));
                sels=zeros(length(theta0),nsamps);
                thet_samples=repmat(theta0',1,1,Ntemps);
                counts=zeros(length(theta0)+length(Z),Ntemps,round(nsamps/step_change));
                thinned_thets=zeros(length(theta0),round(nsamps/Nthin));

    %             step_sizes=zeros(round(Nsamps/step_change), length(theta0));
    %             props=zeros(Nsamps,length(theta0));
    %             sels=zeros(length(theta0),Nsamps);
    %             thet_samples=repmat(theta0',1,Nsamps);
    %             counts=zeros(length(theta0),round(Nsamps/step_change));
    %             thinned_thets=zeros(length(theta0),round(Nsamps/Nthin));
    %   
    xx=2;
else
    xx=nnn;
    Ntemps=1;
    deltaT=.03;
    Nthin=40;
    Nburnin=100;
    logp0=-4.3e4;
    logp = repmat(logp0,[1 5000 Ntemps]);
    acc_count_thet=zeros(length(theta0),Ntemps);
    wacc_count_thet=acc_count_thet;
    lb_thet=modelspec(1).lb;
    ub_thet=modelspec(1).ub;
    thet_samples=repmat(theta0',1,1,Ntemps);
%     wlogpold = logp;
%     logpproposed = logp;
%     wthet = theta0;
%     proposal=wthet;
%     logps=repmat(logp,1,Nsamps);
end

if xx~=nsamps
    for nnn=xx:nsamps
    %         y_samples(:,nn,:)=y_samples(:,nn-1,:);
            thet_samples(:,nnn,:)=thet_samples(:,nnn-1,:);
            logp(:,nnn,:)=logp(:,nnn-1,:);
    %         logpY(:,nn,:)=logpY(:,nn-1,:);
    %    parfor pp=1:Ntemps
    pp=1; %when not doing paralellization
    acc_thet=acc_count_thet(:,pp);
        wacc_thet=wacc_count_thet(:,pp);

         wthet = thet_samples(:,nnn,pp)';
         wlogpold=logp(:,nnn,pp);
%         if nn>1
%             try
%             wlogpold=logprobNI(meantime(trainsub),wYs,dt(trainsub),dY(trainsub),modelspec.traincv(t1,t1,dt1t1,wthet,Ycv(trainsub,trainsub),dy1y1,fp1fp1),mspec,wthet);
%             catch
%                 wlogpold=logp(:,nn,pp);
%             end
%         end
        
        for ii = randperm(length(theta0))
            proposal=wthet;
            proposal(ii) = proposal(ii) + randn*steps(ii);  
            boundsok = (sum(proposal<lb_thet)==0).*(sum(proposal>ub_thet)==0);
            if ~boundsok
%                disp(sprintf('%0.0f -- Round %0.0f - %0.0f -- out of bounds',[pp nn ii]));
            else
                logpproposed=logprobNI(meantime(trainsub),wYs,dt(trainsub),dY(trainsub),modelspec.traincv(t1,t1,dt1t1,proposal,Ycv(trainsub,trainsub),dy1y1,fp1fp1),mspec,proposal);
                %logpproposed=logprobNI(meantime(trainsub),wYs,dt(trainsub),dY(trainsub),modelspec.traincv(t1,t1,dt1t1,proposal,Ycv(trainsub,trainsub),dy1y1,fp1fp1),[],proposal);
%                 [~,~,logpproposed(ii)] = GaussianProcessRegression_ea(meantime(trainsub),wYs,testt,...
%                         modelspec.traincv(t1,t1,dt1t1,proposal,Ycv(trainsub,trainsub),dy1y1,fp1fp1),...
%                         modelspec.cvfunc(t1,t2,dt1t2,proposal,dy1y2,fp1fp2)',...
%                         modelspec.cvfunc(t2,t2,dt2t2,proposal,dy2y2,fp2fp2),[]);
                selector=log(rand);
                cutoff = (logpproposed-wlogpold);%*temps(pp);
%                cutoff = (logpproposed(ii)-wlogpold(ii));%*temps(pp);
                if (selector < cutoff)
%                    disp(sprintf('%0.0f -- Round %0.0f - %0.0f -- accepted %0.2f -- %0.3f (%0.3f) -- %0.3f < %0.3f',[pp nn ii proposal(ii) logpproposed wlogpold exp(selector) exp(cutoff)]));
                    wlogpold = logpproposed;
                    wthet(ii) = proposal(ii);
                    acc_thet(ii) = acc_thet(ii)+1;
                    wacc_thet(ii)= wacc_thet(ii)+1;                
                else
%                    disp(sprintf('%0.0f -- Round %0.0f - %0.0f -- rejected %0.2f-- %0.3f > %0.3f',[pp nn ii proposal(ii) exp(selector) exp(cutoff)]));
                end
            end
        end

%     disp(wlogpold)
     logp(:,nnn,pp)=wlogpold;
     thet_samples(:,nnn,pp)=wthet;
     acc_count_thet(:,pp)=acc_thet;
     wacc_count_thet(:,pp)=wacc_thet;

        if mod(nnn,step_change)==0
            %(nn>50) && (mod(nn,step_change)==0)
            wacc_count=[wacc_count_thet];%; wacc_count_Y];
            disp([sprintf('%0.0f -- ',nnn) sprintf('%0.2f ',sum(wacc_count))]);
%            fprintf('%0.3f -- Avg Accept, nn = %0.0f\n',[mean(acc_count)/nn, nn]);
            step_sizes(round(nnn/step_change),:)=steps;
            sub_big=intersect(find(wacc_count/step_change>.65),find(steps./(ub-lb)<.3));
            sub_small=intersect(find(wacc_count/step_change <.35),find(steps./(ub-lb)>.000001));
            sub_none=intersect(find(wacc_count<1),find(steps./(ub-lb)<.02));
            if size(sub_none)>0
                steps(sub_none)= steps(sub_none)*1.25;                
            end
            if size(sub_big)>0
                steps(sub_big)= steps(sub_big)*1.22;
            end
            if size(sub_small)>0
                steps(sub_small)= steps(sub_small)/1.18;
            end
            sub_none=[]; sub_big=[]; sub_small=[];
%            counts(:,:,round(nn/step_change))=wacc_count;
            wacc_count_thet=0*acc_count_thet;
            %countsY(:,round(nn/step_change))=wacc_countY;
%            wacc_count_Y=0*acc_count_Y;
        end
%             timepassed=toc;
%             fprintf('%0.0f seconds passed -- %0.0f estimated seconds remaining\n',[timepassed timepassed/(nn-1)*(Nsamples-nn)]);        
%             timepassed=toc;
%             fprintf('%0.0f seconds passed -- %0.0f estimated seconds remaining\n',[timepassed timepassed/(nn-1)*(Nsamples-nn)]);        
        if mod(nnn,Nthin)==0
            thinned_thets(:,nnn/Nthin) = thet_samples(:,end,pp);
%            y_samps(:,nn/Nthin,pp)=y_samples(:,end,pp);
            acc_count=[acc_count_thet];%; acc_count_Y];
            if ~isempty(savefile)
                save(savefile,'thinned_thets','acc_count','logps','steps','thet_samples','ii','nnn','modelspec','theta0',...
                    'Nsamples','Nthin','step_sizes','Nburnin','step_change','y_samples','sels','props');
            timepassed=toc;
            fprintf('%0.0f seconds passed -- %0.0f estimated seconds remaining\n',[timepassed timepassed/(nnn-1)*(Nsamples-nnn)]);        
            end
        end
%   if mod(nn,step_change)==0
%         mixers=randperm(Ntemps);
%         for ppp=1
%             m1=mixers(ppp); m2=mixers(ppp+1);
%             logpmix = temps(m1)*logp(:,nn,m2) + temps(m2)*logp(:,nn,m1) - temps(m1)*logp(:,nn,m1) - temps(m2)*logp(:,nn,m2);
%             if log(rand) < logpmix
%                 mixed(nn,[0 1]+ppp) = [m1 m2];
%                 % no longer sampling thetas, but ys
%                 thm1 = thet_samples(:,nn,m1);
%                 thet_samples(:,nn,m1)=thet_samples(:,nn,m2);
%                 thet_samples(:,nn,m2) = thm1;
%                 logpm1 = logp(1,nn,m1);
%                 logp(1,nn,m1)=logp(1,nn,m2);
%                 logp(1,nn,m2)=logpm1;
%                % disp(sprintf('Mixed chains %0.0f and %0.0f',[m1 m2]));
%             else
%                 mixedThet(nn,[0 1]+ppp)=[0 0];
%             end
%         end
%   end
    end
end
