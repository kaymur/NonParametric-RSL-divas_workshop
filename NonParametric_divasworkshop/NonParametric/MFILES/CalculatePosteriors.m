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
if run_type == 1 % sensitivity tests
    if truth_flag==1
        truths = -20000*cos((1950-testt+6280)/2000)-20000;  
        truerate=[1.2277;    1.7209; 2.2098;    2.6931; 3.1697; 3.6384;    4.0980 ;4.5474;4.9853; 5.4109;5.8229;6.2203; 6.6022; 6.9676; 7.3156 ;   7.6453;    7.9559;    8.2466;    8.5166;    
                8.7654;    8.9923;    9.1968;    9.3782;    9.5362;    9.6703;    9.7803;    9.8658;    9.9267;    9.9628;    9.9739;    9.9602;    9.9215;    9.8581;    9.7700;    9.6574;    
                9.5208;    9.3603;    9.1765;    8.9697;    8.7405;    8.4894;    8.2171;    7.9243;    7.6117;    7.2801;    6.9302;    6.5631;    6.1795;    5.7805;    5.3670;    4.9402;    
                4.5009;    4.0505;    3.5899;    3.1203;    2.6429;    2.1590;    1.6696;    1.1760;    0.6796;    0.1814;   -0.3173;   -0.8151;   -1.3109;   -1.8034;   -2.2915;   -2.7738;   
                -3.2491;   -3.7164;   -4.1743;   -4.6218;   -5.0578;   -5.4811;   -5.8907;   -6.2856;   -6.6648;   -7.0273;   -7.3723;   -7.6989;   -8.0061;   -8.2934;   -8.5600;   -8.8051;   
                -9.0283;   -9.2289;   -9.4064;   -9.5604;   -9.6905;   -9.7964;   -9.8778;   -9.9345;   -9.9664;   -9.9734;   -9.9554;   -9.9126;    -9.8450;   -9.7527;   -9.6361;   -9.4955; 
                -9.3310;   -9.1433;   -8.9327;   -8.6998;   -8.4451;   -8.1693;   -7.8731;   -7.5573;   -7.2225;   -6.8697;   -6.4997;   -6.1135;   -5.7120;   -5.2962;   -4.8671;   -4.4259;  -3.9737];    
    elseif truth_flag==2
        sub0=find(1950-testt>7000);
        sub1=find(1950-testt<=7000 & 1950-testt>6000);
        sub2=find(1950-testt<=6000 & 1950-testt>5200);
        sub3=find(1950-testt<=5200 & 1950-testt>3800);
        sub4=find(1950-testt<=3800 & 1950-testt>3000);
        sub5=find(1950-testt<=3000);
        Y0=testt*0; Y1=Y0; Y2=Y0; Y3=Y0; Y4=Y0; Y5=Y0;

        Y5(sub5)=(1950-testt(sub5))*-20;
        Y4(sub4)=-60e3+(1950-testt(sub4)-3000)*-2.5;
        Y3(sub3)=-62e3+(1950-testt(sub3)-3800)*-18;
        Y2(sub2)=-87200+(1950-testt(sub2)-5200)*-1.2;
        Y1(sub1)=-88160+(1950-testt(sub1)-6000)*-40;
        Y0(sub0)=-128160+(1950-testt(sub0)-7000)*-5;
        truths =Y0+Y1+Y2+Y3+Y4+Y5;
        truerate=[20.0000;  20.0000;   20.0000;   20.0000;   20.0000;  20.0000; 20.0000;   20.0000;   20.0000;   20.0000;   20.0000;   20.0000;   20.0000;   20.0000;   20.0000;   20.0000;  20.0000;    
                20.0000;   20.0000;   20.0000;   20.0000;  20.0000; 20.0000;   20.0000;   20.0000;   20.0000;  16.5000;    13.0000;    9.5000;    6.0000;    2.5000;    2.5000;    2.5000;    2.5000;   
                5.6000;    8.7000;   11.8000;   14.9000;   18.0000;   18.0000;   18.0000;   18.0000;   18.0000;   18.0000;   18.0000;   18.0000;   18.0000;   18.0000;   14.6400;   11.2800;    7.9200;    
                4.5600;    1.2000;    1.2000;    1.2000;    1.2000;    8.9600;   16.7200;   24.4800;   32.2400;   40.0000;   40.0000;   40.0000;   40.0000;  40.0000;   40.0000;   33.0000;   26.0000;   
                19.0000;   12.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    
                5.0000;    5.0000;   5.0000;    5.0000;  5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    
                5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000;    5.0000];
    end
else
    truths=interp3(RSLs.Lat,RSLs.Lon,RSLs.Time,RSLs.Y0,testlocs.X(:,1),mod(testlocs.X(:,2),360),testlocs.X(:,3),'linear')';
end

% coverage at testtimes
cov_95=(length(truths')-(sum(truths'>f95u)+sum(truths'<f95l)))/length(truths');
cov_67=(length(truths')-(sum(truths'>f67u)+sum(truths'<f67l)))/length(truths');

ave_95_wi=mean(f95u-f95l);
ave_67_wi=mean(f67u-f67l);
