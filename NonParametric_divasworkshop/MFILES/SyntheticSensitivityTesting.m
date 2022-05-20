%% n is the number of data points per 1000 year period (n=1,5, or 10)
%% t_unc is the median amount of age uncertainty (75 or 225 years)
%% dat_type defines whether we are using the sedimentary (3), acropora(1), or orbicella(2) distribution

if dat_type == 9
    dat_type_0 = dat_type;
    dat_type = 1;
elseif dat_type == 10
    dat_type_0 = dat_type;
    dat_type = 2;
else
    dat_type_0 = 0;
end

%% create data
if jjj>48
    jj0=jjj;
    jjj=jjj-48;
end

[ages,dt,age_true,Y_final,Y_synth,Y_true,dY_synth,limiting,cspecies,meas_err,offset]=CreateSensitivityData(t_unc,n,[0 12000],truth_flag,distKern0,limiting,cspecies,dat_type,jjj,dateField,Seed);

if exist('jj0','var')
    jjj=jj0;
end
% PlotSensitivityData;

datid = datid(1:N); 
Y0 = Y0(1:N); 
dY0 = dY0(1:N); 
dY_init = dY_init(1:N); 
lat = lat(1:N); 
long = long(1:N); 
meanYs = Y0(1:N); 
meantime = 1950-ages;
time1 = meantime-2*dt;
time2 = meantime+2*dt;
%limiting = limiting(1:N);
Y = Y_final;
dY = dY_synth;
testt0 = [0:100:12000];

Ycv = sparse(diag(dY.^2));

if truth_flag==2
    meanSL = 0;
else
    meanSL=0;
end
Y = Y - meanSL;
meanYs = meanYs - meanSL;

PXholo.datid=round(datid);
PXholo.region=region;
PXholo.time1=time1;
PXholo.time2=time2;
PXholo.limiting=limiting;
PXholo.Y=Y;
PXholo.Y0=Y0;
PXholo.dY = dY;
PXholo.dY0=dY0;
PXholo.cspecies=cspecies;
PXholo.lat=lat;
PXholo.long=long;
PXholo.Ycv=sparse(diag(dY.^2));
PXholo.siteid=round(siteid);
PXholo.sitenames=sitenames;
PXholo.meantime=meantime;
PXholo.sitecoords = sitecoords;
PXholo.sitelen = N;
PXholo.dt=dt;

clear datasets;

datasets{1}=PXholo;
datasets{1}.label='PXholo';

for ii=1:length(datasets)
    datasets{ii}.Y0=datasets{ii}.Y0;
    datasets{ii}.Y=datasets{ii}.Y;
    t1=datasets{ii}.time1; t2=datasets{ii}.time2;
    datasets{ii}.long = mod(datasets{ii}.long,360); sub=find(datasets{ii}.long>180); datasets{ii}.long(sub)=datasets{ii}.long(sub)-360;
    datasets{ii}.meantime=mean([t1 t2],2);
    datasets{ii}.dt = abs(t1-t2)/4;
end

DefineCovarianceFunctions;

trainsub = find(Y); 
if dat_type==3 
    trainsubz=find(limiting~=20); sampnorm=1;
else
    trainsubz = find(limiting~=0&limiting~=20);
end
trainsublim = find(limiting==3|limiting==-1|limiting==1);
trainsubNotNormal = find(limiting==3|limiting==-1|limiting==1|limiting==2);
angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));
dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);
dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
obsGISfp=ones(size(t1));
fp1fp1=bsxfun(@times,obsGISfp-1,obsGISfp'-1);
thetTGG=[modelspec(1).thet0];
subfixed = modelspec(1).subfixed;

subnotfixed=setdiff(1:length(thetTGG),subfixed);

%% find optimized thetas with sampled data.
dolb = [modelspec(ii).lb]; %1e-6]; lower bound
doub = [modelspec(ii).ub]; % 1e-5]; upper bound

%starting values
initialm = thetTGG;
%%%%%

%thetas for amplitude, age, and white noise
s =  [850   130  315]; 
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

t1=X1(:,3);
t2=testXs(:,3);

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

trainsub0=find(limiting==0);
thetPXholo=modelspec(1).thet0;
theta0=thetPXholo;
% 
t1=X1(:,3);
dt12=dYears(t1(trainsub0),t1(trainsubz));
dy12=dDist(X1(trainsub0,1:2),X1(trainsubz,1:2));
fp12=bsxfun(@times,obsGISfp(trainsub0)-1,obsGISfp(trainsubz)'-1)';
dt22=dYears(t1(trainsubz),t1(trainsubz));
dy22=dDist(X1(trainsubz,1:2),X1(trainsubz,1:2));
fp22=bsxfun(@times,obsGISfp(trainsubz)-1,obsGISfp(trainsubz)'-1)';
[Ys0,~,~] = GaussianProcessRegression(meantime(trainsub0),Y(trainsub0),meantime(trainsubz),...
        modelspec.traincv(t1(trainsub0),t1(trainsub0),dt1t1(trainsub0,trainsub0),thetPXholo,Ycv(trainsub0,trainsub0),dy1y1(trainsub0,trainsub0),fp1fp1(trainsub0,trainsub0)),...
        modelspec.cvfunc(t1(trainsub0),meantime(trainsubz),dt12,thetPXholo,dy12,fp12)',...
        modelspec.cvfunc(meantime(trainsubz),meantime(trainsubz),dt22,thetPXholo,dy22,fp22),[]);
first_flag=1;
steps0=dY;

%% set up the sampling
tic
step_size=dY; 
if synthflag ==0
else
    Ys0 = meanSL + zeros(size(trainsubz,1),1);
end

tic
train_all=find(limiting~=20);
dataset=datasets{1};
% if dat_type_0 == 9 || dat_type_0 == 10
% %    modno=6; nsamps = 25000;
%     trainsubx=train_all;
%     
%     sub_unif=find(cspecies==1);
%     delev(sub_unif)=dY(sub_unif);
%     elev(sub_unif)=Y(sub_unif);
%     dY(sub_unif)= sqrt((dY(sub_unif)).^2+(2500)^2);  
%     Y(sub_unif)=Y(sub_unif)+2500;
%     limiting(sub_unif)=0;
%     
%     sub_orb=find(cspecies==9);
%     delev(sub_orb)=dY(sub_orb);
%     elev(sub_orb)=Y(sub_orb);
%     dY(sub_orb)= sqrt((dY(sub_orb)).^2+(10000)^2);  
%     Y(sub_orb)=Y(sub_orb)+10000;
%     limiting(sub_orb)=0;
%     
%     sampnorm=0;
%     nsamps = 100000;
%     
%     %CompareNormal;
if max(cspecies) == 0 && isempty(trainsublim) && sampnorm==0
    wdataset = dataset;
    clear testsitedef;
    testsitedef.sites=[];
    testsitedef.names={};
    testsitedef.names2={};
    testsitedef.firstage=[];
    ii = 1;
    si=find(datasets{1}.datid==datasets{1}.siteid(ii)); si=si(1);
    testsitedef.sites(end+1,:)=[datasets{1}.datid(si) datasets{1}.lat(si) datasets{1}.long(si)];
    testsitedef.names2={testsitedef.names2{:}, datasets{1}.sitenames{ii}};
    testsitedef.names={datasets{1}.sitenames{ii}};
    testsitedef.firstage = [testsitedef.firstage -12000];
    noiseMasks = ones(1,length(theta0));
    [f2s,sd2s,V2s,testloc]=RegressHoloceneDataSets(wdataset,testsitedef,modelspec(1),theta0,trainsub,noiseMasks,testt,refyear,[]);

    %% only sample the hyperparameter values for the normal model
    sample_Theta;
elseif dat_type_0 == 9 || dat_type_0 == 10
    sample_Ys_Thetas_Norm;
else
    %% sample sea levels and hyperparameters
    sample_Ys_Thetas;
end

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

%% post processing (can activate any of these to find further

ModelEvaluation;

Y_orig = Y+meanSL;

