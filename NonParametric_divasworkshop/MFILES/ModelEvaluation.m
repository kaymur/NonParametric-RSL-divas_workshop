%% the probability of true sea level, given the posterior model prediction (for Gaussian and non-parametric) and given the data distribution (for uncorrelated)
%% define each point to find a posterior prediction on sea level
% trueYs are the true sea level at each true time age and corresponding lat/long

mod_eval=1;
if ~exist('noiseMasks')
    if exist('thetPXholo')
        noiseMasks = ones(1,length(thetPXholo));
    else
        noiseMasks = ones(1,length(theta0));
    end
end
noiseMasks(end)=0;
if exist('GIA0')>0 % then these tests are synthetic, and we need to read the "true RSLs" and compare
    RSLs.Y0=GIA0;
    RSLs.Lat=GIALat;
    RSLs.Lon=GIALon;
    RSLs.Time=GIAtime;    
end
if run_type == 2 % validation tests
    testt = 1950-[0:100:12000];
end

CalculatePosteriors;

if ~exist('wdataset')
    wdataset=datasets{1};
end
clear fp sdp testtsp;
clear testsitedefp;
if run_type == 2
    u=unique(wdataset.datid);
    subps=[];
    for pp=1:length(u)
        subp=find(wdataset.datid==u(pp));
        subq=find(wdataset.siteid==u(pp));
        subq=subq(1);
        if length(subp)>0
            testtsp{pp}=age(subp);
            testsitedefp.sites(pp,:)=[wdataset.siteid(subq) ...
                                wdataset.sitecoords(subq,:)];
            testsitedefp.names(pp)=wdataset.sitenames(subq);
            testsitedefp.names2(pp)=wdataset.sitenames(subq);
            testsitedefp.firstage(pp)=min(wdataset.meantime(subp))-1950;
        end
        subps=[subps ; subp];
    end
    [fp(subps),sdp(subps),~,testlocp]=RegressHoloceneDataSets(wdataset,testsitedefp,modelspec(1),thets(:,1)',trainsubx,noiseMasks(1,:),testtsp,refyear,[]);

    if size(thets,2)<2
        fs=zeros(length(testlocp.reg),size(ys,2));
        sds=zeros(length(testlocp.reg),size(ys,2));
        fsp_all0=zeros(NpostSamps*size(thets,2),length(testlocp.reg));
        for iii=1:size(ys,2)
            if ~isempty(thinned_ts)
                [fsp(:,iii),sdsp(:,iii),V2sp,testlocp,logp_test]=RegressHoloceneDataSets_sampled(wdataset,testsitedefp,modelspec,thets(:,iii)',trainsubx,noiseMasks(:,:),testtsp,refyear,collinear,[],[],ys(:,iii),thinned_ts(:,iii),trainsubx);
            else
                [fsp(:,iii),sdsp(:,iii),V2sp,testlocp,logp_test]=RegressHoloceneDataSets_sampled(wdataset,testsitedefp,modelspec,thets(:,iii)',trainsubx,noiseMasks(:,:),testtsp,refyear,collinear,[],[],ys(:,iii),[],trainsubx);
            end
            try
                L=chol(V2sp,'lower');
            catch
                logp_test=-1e12;
            end
            if logp_test < -1e10
            else
                rng(Seed*iii);
                fsp_post0=mvnrnd(fsp(:,iii),V2sp,NpostSamps);
                hold on;
                fsp_all0((iii-1)*NpostSamps+1:iii*NpostSamps,:)=fsp_post0;
            end
        end    
    else
        [fp,sp,~,testlocp,logp_test]=RegressHoloceneDataSets_sampled(wdataset,testsitedefp,modelspec,thets(:,1)',trainsubx,noiseMasks(1,:),testtsp,refyear,[],[],[],[],[],trainsubx);
        fsp=zeros(length(testlocp.reg),size(thets,2));%*size(ys,2));
        sdsp=zeros(length(testlocp.reg),size(thets,2));%*size(ys,2));
        fsp_all0=zeros(NpostSamps*size(thets,2),length(testlocp.reg));
        for iii=1:size(thets,2)
            if ~isempty(thinned_ts)
                [fsp(:,iii),sdsp(:,iii),V2sp,testlocp,logp_test]=RegressHoloceneDataSets_sampled(wdataset,testsitedefp,modelspec,thets(:,iii)',trainsubx,noiseMasks(1,:),testtsp,refyear,collinear,[],[],ys(:,iii),thinned_ts(:,iii),trainsubx);
            elseif ~isempty(ys)
                [fsp(:,iii),sdsp(:,iii),V2sp,testlocp,logp_test]=RegressHoloceneDataSets_sampled(wdataset,testsitedefp,modelspec,thets(:,iii)',trainsubx,noiseMasks(1,:),testtsp,refyear,collinear,[],[],ys(:,iii),[],trainsubx);
            else
                [fsp(:,iii),sdsp(:,iii),V2sp,testlocp,logp_test]=RegressHoloceneDataSets_sampled(wdataset,testsitedefp,modelspec,thets(:,iii)',trainsub,noiseMasks(1,:),testtsp,refyear,collinear,[],[],[],[],trainsub);
            end
            if logp_test < -1e11
            else
                try
                    L=chol(V2sp,'lower');
                catch
                    logp_test=-1e10;
                end
                if logp_test < -1e9
                    rng(Seed*iii);
                    fsp_post0=mvnrnd(fsp(:,iii),sdsp(:,iii)',NpostSamps);
                    hold on;
                    fsp_all0((iii-1)*NpostSamps+1:iii*NpostSamps,:)=fsp_post0;
                else
                    rng(Seed*iii);
                    fsp_post0=mvnrnd(fsp(:,iii),V2sp,NpostSamps);
                    hold on;
                    fsp_all0((iii-1)*NpostSamps+1:iii*NpostSamps,:)=fsp_post0;
                end
            end
        end
    end
    toc
    sub_bad=find(fsp_all0(:,1)==0);
    sub_good=find(fsp_all0(:,1)~=0);
    
    %% this is the distribution on which the true RSL will be evaluated
    fsp_all=fsp_all0(sub_good,:);
    maxf1=(quantile(fsp_all,.975)+meanSL);
    minf1=(quantile(fsp_all,.025)+meanSL);
    mxf1=max(fsp_all);
    mnf1=min(fsp_all);
    medf1=(median(fsp_all)+meanSL);
    meanf1=(mean(fsp_all)+meanSL);
    f95u1=(quantile(fsp_all,.975)+meanSL);
    f95l1=(quantile(fsp_all,.025)+meanSL);
    f67u1=(quantile(fsp_all,.833)+meanSL);
    f67l1=(quantile(fsp_all,.167)+meanSL);
end
if run_type == 2 % validation tests
    truth=MakePlotsSamplingVal(wdataset,medf1,maxf1,minf1,testlocp,['ModEval_modno_' num2str(modno) '_'],1000,[0 11000],meanSL,testt',coral_params,RSLs,Seed,modno);
    for jj = 1:length(testlocp.reg)
        % find true sea level at points where there are data
        scoord(1) = testlocp.X(jj,1);                
        scoord(2) = mod(testlocp.X(jj,2),360);
        timeGIA = testlocp.X(jj,3);
        trueSLs(jj)=interp3(GIALat,GIALon,GIAtime,GIA0,scoord(1),scoord(2),timeGIA);
    end
    allY=find(Y); 
    likes=zeros(size(Y(allY),1),1);
    for jj = 1:length(likes)
        likes(jj)=ksdensity(fsp_all(:,allY(jj))',trueSLs(allY(jj)));
    end
    badll=find(likes<1e-300);
    subll=find(likes>1e-300);
    like=likes;
    like(badll)=1e-300;
    like(subll)=likes(subll);
    fin_ll=sum(log(like));
    clf;

    resids = median(fsp_all(:,allY))-trueSLs(allY);
    m_resid=mean(resids);
    med_resid=median(resids);
    abs_mean_resid=mean(abs(resids));
    srmse = sqrt(sum(resids.^2)/length(resids));
    plot(testlocp.X(allY,3),resids,'.');
    pdfwrite('residuals_v_time');
    plot(trueSLs(allY),resids,'.');
    pdfwrite('residuals_v_RSL');
    
    medf1=(median(fsp_all)+meanSL);
    meanf1=(mean(fsp_all)+meanSL);
    f95u1=(quantile(fsp_all,.975)'+meanSL);
    f95l1=(quantile(fsp_all,.025)'+meanSL);
    f67u1=(quantile(fsp_all,.833)'+meanSL);
    f67l1=(quantile(fsp_all,.167)'+meanSL);
    if ~exist('Y_true')
        Y_true = trueYs;
    end
    resids = median(fsp_all(:,allY))'-Y_true(allY);
    m_resid=mean(resids);
    med_resid=median(resids);
    abs_mean_resid=mean(abs(resids));
    srmse = sqrt(sum(resids.^2)/length(resids));
else
    clf;
    
    resids = medf-truths';
    m_resid=mean(resids);
    med_resid=median(resids);
    abs_mean_resid=mean(abs(resids));
    srmse = sqrt(sum(resids.^2)/length(resids));
    
    plotdat.x=1950-testt';
    plotdat.y=medf;
    plotdat.dy=[f95u-medf medf-f95l] ;
    clf; 
    PlotWithShadedErrors(plotdat,[0 0 0],'none');
    hold on;
    if truth_flag==1
        ylim([-50e3 25e3]);
    else
        ylim([-180e3 15e3]);
    end
    labl = ['Post RSL ' num2str(jjj) 'Seed ' num2str(Seed) ' RSL Curve ' num2str(truth_flag)];
    title=labl;

    MakePlotsSampling(PXholo,fs(:,1),sds(:,1),V2s,[],labl,[],1000,[],[],meanSL,1950-testt',coral_params);

    if truth_flag==1
        ylim([-65e3 25e3]);
    else
        ylim([-180e3 15e3]);
    end
    testX = [mean(lat)*ones(length(testt),1) mean(long)*ones(length(testt),1) 1950-testt'];
    
    plot(testX(:,3),truths,'r','LineWidth',1);
    plot(testX(:,3),f95u,'--','Color',[.5, .5, .5]);
    plot(testX(:,3),f95l,'--','Color',[.5, .5, .5]);
    if ~exist('jj', 'var')
        jj=jjj;
    end
    pdfwrite(['Synth_post_RSL_' df '_' num2str(jj) '_' num2str(Seed) '_' num2str(truth_flag)]);        
   
    PredictRates;
end

% RunStepBurnThinAnalysis;
