%times over which you are making predictions of SL
testt0 = [0:100:12000];
testt = 1950-testt0;
datt = [meantime; testt'];

X1 = [lat long meantime];
testX = [mean(lat)*ones(length(testt),1) mean(long)*ones(length(testt),1) testt'];
testXs = [X1; testX];
testfp = [ones(length(testt),1); obsGISfp(:)];

dt1t1 = dYears(X1(:,3),X1(:,3));
dt1t2 = dYears(X1(:,3),testXs(:,3));
dt2t2 = dYears(testXs(:,3),testXs(:,3));

dy1y1 = dDist(X1(:,1:2),X1(:,1:2));
dy1y2 = dDist(X1(:,1:2),testXs(:,1:2));
dy2y2 = dDist(testXs(:,1:2),testXs(:,1:2));

fp1fp1=bsxfun(@times,obsGISfp(:)-1,obsGISfp(:)'-1)';
fp1fp2=bsxfun(@times,obsGISfp(:)-1,testfp'-1)';
fp2fp2=bsxfun(@times,testfp-1,testfp'-1);


%% posteriors only for test points testt
% dt1dt2 = dYears(X1(:,3),testX(:,3));
% dt2dt2 = dYears(testX(:,3),testX(:,3));
% dy1dy2 = dDist(X1(:,1:2),testX(:,1:2));
% fp2p2 = dDist(testX(:,1:2),testX(:,1:2));
testfp2 = [ones(length(testt),1)];
% fp1p2=bsxfun(@times,obsGISfp(:)-1,testfp2'-1)';
% fp2p2=bsxfun(@times,testfp2-1,testfp2'-1);
% t22=testX(:,3);
% dy2dy2=dDist(testX(:,1:2),testX(:,1:2));

%% posteriors only for test points testt
dt1dt2 = dYears(X1(:,3),testX(:,3));
dt2dt2 = dYears(testX(:,3),testX(:,3));
dy1dy2 = dDist(X1(:,1:2),testX(:,1:2));
fp2p2 = dDist(testX(:,1:2),testX(:,1:2));
testfp2 = [ones(length(testt),1)];
fp1p2=bsxfun(@times,obsGISfp(:)-1,testfp2'-1)';
fp2p2=bsxfun(@times,testfp2-1,testfp2'-1);
t22=testX(:,3);
dy2dy2=dDist(testX(:,1:2),testX(:,1:2));



logps=logp;
theta_samp=thet_samples;
t0 = meantime(train_all);
dt0t0=dYears(t0,t0);
dt0t2 = dYears(meantime(train_all),testX(:,3)); 
dy0dy0 = dDist(X1(train_all,1:2),X1(train_all,1:2));
dy0dy2 = dDist(X1(train_all,1:2),testX(:,1:2));
fp0p0 = bsxfun(@times,obsGISfp(train_all)-1,obsGISfp(train_all)'-1)';
fp0p2 = bsxfun(@times,obsGISfp(train_all)-1,testfp2'-1)';
alfa_flag=0;
testcv0 = @(thet) modelspec.cvfunc(meantime(train_all),testX(:,3),dYears(meantime(train_all),testX(:,3)),thet,dy0dy2,fp0p2);
testcv02 = @(thet) modelspec.cvfunc(testX(:,3),testX(:,3),dYears(testX(:,3),testX(:,3)),thet,dy2dy2,fp2p2);

%% when normal is not sampled, add those in to be used in the regression
y_accepted=[];
% if jjj>48 || (max(cspecies) == 0 && isempty(trainsublim) && sampnorm == 0)
%     ys0=Y(trainsub);
%     ys=[];
%     thinned_ts=[];
%     thets=thinned_thets(:,nburn:end);
% else
    %% includes Ys for normally distributed and sampled values for sampled thinned_ys;
    y_int=repmat(Y,1,size(thinned_ys,2));
    y_int(trainsubz,:)=thinned_ys;
    y_accepted=y_int(train_all,nburn:end);
    trainsubx=train_all;
% end
sub_all = intersect(find(limiting~=1),find(limiting~=-1));
Y0s=repmat(Y(train_all),1,size(y_accepted,2));
Y0s=y_accepted;
thinned_ts=[];
thets=thinned_thets(:,nburn:end);
trainsubx=train_all;



% calculate_posteriors
refyear=1950;
if ~exist('Seed','var')
    Seed = 1;
end
if run_type == 2 % validation tests
    if ~exist('y_accepted')
        ys0=Y(trainsub);
        ys=[];
        thetas=thets;
    elseif length(y_accepted)==0
        ys0=Y(trainsub);
        ys=[];
        thetas=thets;
    else
        ys=y_accepted;
        thetas=thets;
    end
else
    if ~exist('fs_all')
        if ~exist('y_accepted')
            ys0=Y;
            ys=[];
            thetas=thets;
            y_accepted = ys0;
        elseif isempty(y_accepted)
            ys0=Y;
            ys=[];
            thetas=thets;
            y_accepted = ys0;
        else
            ys0=y_accepted;
            thetas=thets;
        end
        ys = y_accepted;
    else
        ys = y_accepted;
    end
end

%% now do a regression
clear Loc cnt Locid oldest youngest testsitedef;
distfrom = dDist(datasets{1}.sitecoords,datasets{1}.sitecoords);
for ii=1:length(datasets{1}.sitenames)
    s0=0;
    Loc{ii}=datasets{1}.sitenames{ii};
    sub=strfind(Loc{ii},' ');
    Loc{ii}=Loc{ii}(setdiff(1:length(Loc{ii}),sub));
    sub=strfind(Loc{ii},'/');
    Loc{ii}=Loc{ii}(setdiff(1:length(Loc{ii}),sub));
    sub=find(datasets{1}.datid==datasets{1}.siteid(ii));
    cnt(ii)=length(sub);
end

clear testsitedef;
testsitedef.sites=[];
testsitedef.names={};
testsitedef.names2={};
testsitedef.firstage=[];

for ii=1:length(Loc)
        si=find(datasets{1}.datid==datasets{1}.siteid(ii)); si=si(1);
        testsitedef.sites(end+1,:)=[datasets{1}.datid(si) datasets{1}.lat(si) datasets{1}.long(si)];
        testsitedef.names2={testsitedef.names2{:}, datasets{1}.sitenames{ii}};
        testsitedef.names={testsitedef.names{:}, Loc{ii}};
        testsitedef.firstage = [testsitedef.firstage -12000];
end
trainspecs=1;
trainlabels={'default'};

regresssets=[1];
regressparams=[1];
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end

iii=1;
ii=regresssets(iii);
clear wdataset;
noiseMasks = ones(1,length(theta0));
noiseMasks(1,end)=0;
defaultmask=1;
wdataset=datasets{1,1};
labls{iii}=['_' regresslabels{iii}];
disp(labls{iii});
trainsub = find(wdataset.limiting==0);
subtimes=find(testt>=-4000+min(wdataset.time1));
collinear=[];
if ~exist('mod_eval','var')
    [f0,s0,V0,testlocs,logpp]=RegressHoloceneDataSets(wdataset,testsitedef,modelspec(1),thetas(:,1)',trainsub,noiseMasks(:,:),testt,refyear,collinear);
    tsdef=testsitedef;
else
    u=unique(wdataset.datid);
    clear fp sdp;
    clear testsitedefp;
    subps=[];
    for pp=1:length(u)
        subp=find(wdataset.datid==u(pp));
        subq=find(wdataset.siteid==u(pp));
        subq=subq(1);
        if length(subp)>0
            testtsp{pp}=wdataset.meantime(subp);
            testsitedefp.sites(pp,:)=[wdataset.siteid(subq) ...
                                wdataset.sitecoords(subq,:)];
            testsitedefp.names(pp)=wdataset.sitenames(subq);
            testsitedefp.names2(pp)=wdataset.sitenames(subq);
            testsitedefp.firstage(pp)=min(wdataset.meantime(subp));
        end
        subps=[subps ; subp];
    end
    wdataset.istg = 0*wdataset.Y;
    wdataset.compactcorr = 0*wdataset.Y;
    [f0,s0,V0,testlocs,logpp]=RegressHoloceneDataSets(wdataset,testsitedef,modelspec(1),thetas(:,1)',trainsub,noiseMasks(:,:),testt,refyear,collinear);
    tsdef=testsitedef;
end
collinear =[]; NpostSamps=10; V2s=[]; fs=[]; sds=[];
tic
clf;
if size(thetas,2)<2
    fs=zeros(length(testlocs.reg),size(ys,2));
    sds=zeros(length(testlocs.reg),size(ys,2));
    fs_all=[];
    fs_all0=zeros(NpostSamps*size(thetas,2),length(testlocs.reg));
    for iii=1:size(ys,2)
        if ~isempty(thinned_ts)
                [fs(:,iii),sds(:,iii),V2s,testlocs,logp_test]=RegressHoloceneDataSets_sampled(wdataset,tsdef,modelspec(1),thetas(:,iii)',trainsubx,noiseMasks(:,:),testt,refyear,collinear,[],[],ys(:,iii),thinned_ts(:,iii),trainsubx);
        else
                [fs(:,iii),sds(:,iii),V2s,testlocs,logp_test]=RegressHoloceneDataSets_sampled(wdataset,tsdef,modelspec(1),thetas(:,iii)',trainsubx,noiseMasks(:,:),testt,refyear,collinear,[],[],ys(:,iii),[],trainsubx);
        end
        try
            L=chol(V2s,'lower');
        catch
            logp_test=-1e12;
        end
        if logp_test < -1e10
        else
            rng(Seed*iii);
            fs_post0=mvnrnd(fs(:,iii),V2s,NpostSamps);
            hold on;
            fs_all0((iii-1)*NpostSamps+1:iii*NpostSamps,:)=fs_post0;
        end
    end    
else
    fs=zeros(length(testlocs.reg),size(thetas,2));%*size(ys,2));
    sds=zeros(length(testlocs.reg),size(thetas,2));%*size(ys,2));
    fs_all0=zeros(NpostSamps*size(thetas,2),length(testlocs.reg));
    for iii=1:size(thetas,2)
        if ~isempty(thinned_ts)
            [fs(:,iii),sds(:,iii),V2s,testlocs,logp_test]=RegressHoloceneDataSets_sampled(wdataset,tsdef,modelspec(1),thetas(:,iii)',trainsubx,noiseMasks(:,:),testt,refyear,collinear,[],[],ys(:,iii),thinned_ts(:,iii),trainsubx);
        elseif ~isempty(ys) && size(ys,2)>1
                [fs(:,iii),sds(:,iii),V2s,testlocs,logp_test]=RegressHoloceneDataSets_sampled(wdataset,tsdef,modelspec(1),thetas(:,iii)',trainsubx,noiseMasks(:,:),testt,refyear,collinear,[],[],ys(:,iii),[],trainsubx);
        else
            %% for models that don't sample the data.  Only the limiting==0 trainsub should be used.
            [fs(:,iii),sds(:,iii),V2s,testlocs,logp_test]=RegressHoloceneDataSets_sampled(wdataset,tsdef,modelspec(1),thetas(:,iii)',trainsub,noiseMasks(:,:),testt,refyear,collinear,[],[],[],[],trainsub);
        end
        if logp_test < -1e11
        else
            try
                L=chol(V2s,'lower');
            catch
                logp_test=-1e10;
            end
            if logp_test < -1e9
                rng(Seed*iii);
                fs_post0=mvnrnd(fs(:,iii),sds(:,iii)',NpostSamps);
                hold on;
                fs_all0((iii-1)*NpostSamps+1:iii*NpostSamps,:)=fs_post0;
            else
              rng(Seed*iii);
              fs_post0=mvnrnd(fs(:,iii),V2s,NpostSamps);
                hold on;
                fs_all0((iii-1)*NpostSamps+1:iii*NpostSamps,:)=fs_post0;
            end
        end
    end
end
toc
sub_bad=find(fs_all0(:,1)==0);
sub_good=find(fs_all0(:,1)~=0);
fs_all=fs_all0(sub_good,:);
toc

maxf=max(fs_all)';
minf=min(fs_all)';
if ~exist('meanSL')
    meanSL = 0;
end
medf=(median(fs_all)+meanSL)';
meanf=(mean(fs_all)+meanSL)';
f95u=(quantile(fs_all,.975)+meanSL)';
f95l=(quantile(fs_all,.025)+meanSL)';
f67u=(quantile(fs_all,.833)+meanSL)';
f67l=(quantile(fs_all,.167)+meanSL)';













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

plotdat.x=1950-testt';
plotdat.y=medf;
plotdat.dy=[f95u-medf medf-f95l] ;
clf; 
PlotWithShadedErrors(plotdat,[0 0 0],'none');

hold on;
labl = ['Posterior RSL Curve ' dateField];
title=labl;

MakePlotsSampling(PXholo,fs(:,1),sds(:,1),V2s,[],labl,[],1000,[],[],meanSL,1950-testt',coral_params);
