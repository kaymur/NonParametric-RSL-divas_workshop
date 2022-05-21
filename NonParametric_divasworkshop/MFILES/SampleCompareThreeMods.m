    if modno==6
        dY_init=dY;
        excl=intersect(find(cspecies>1),find(cspecies<10));
        limiting(excl)=20;
        excl=find(limiting==-1);
        limiting(excl)=20;
        excl=find(limiting==1);
        limiting(excl)=20;
        sub_unif=find(cspecies==1);
        dY(sub_unif)= (2/(sqrt(3))) * dY_init(sub_unif);
        limiting(sub_unif)=0;
    end
        Ycv = sparse(diag(dY.^2));
        PXholo.datid=round(datid);
        PXholo.region=region;
        PXholo.time1=time1;
        PXholo.time2=time2;
        PXholo.limiting=limiting;
        PXholo.indic=indic;
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
        %PXholo.meantime=(PXholo.time1+PXholo.time2)/2;
        PXholo.meantime=meantime;
        PXholo.sitecoords = sitecoords;
        PXholo.sitelen = sitelen;
        PXholo.dY_init=dY_init;
        PXholo.indic=indic;
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
        
%%%%%

if spatial==1 && (modno==6||modno==3)
    DefCov_spatial;
else
    DefineCovarianceFunctions;
end

trainsub = find(Y); 
if sampnorm ==1
    trainsubz=find(limiting~=20);
else
    trainsubz = find(limiting~=0&limiting~=20);
end

trainsublim = find(limiting==3|limiting==-1|limiting==1);
trainsubNotNormal = find(limiting==3|limiting==-1|limiting==1|limiting==2);

obsGISfp=ones(size(time1));
fp1fp1=bsxfun(@times,obsGISfp-1,obsGISfp'-1);

%defines starting theta values
thetTGG=[modelspec(1).thet0];
subfixed = modelspec(1).subfixed;

subnotfixed=setdiff(1:length(thetTGG),subfixed);

dolb = [modelspec(1).lb]; %1e-6]; lower bound
doub = [modelspec(1).ub]; % 1e-5]; upper bound

%starting values
initialm = thetTGG;

%%%%%
%thetas for amplitude, age, and white noise
%can use sqrt of range of upper and lower bound to determine starting
%parameters for s
%standard deviation of the jumps in hyperparameters; 
%determines acceptance ratio-ideal is ~40%
%this is something we would likely manipulate a lot
%%%%%

s =  [850   130  315]; 

%times over which you are making predictions of SL
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

t1=X1(:,3);
dt12=dYears(t1(trainsub0),t1(trainsubz));
dy12=dDist(X1(trainsub0,1:2),X1(trainsubz,1:2));
fp12=bsxfun(@times,obsGISfp(trainsub0)-1,obsGISfp(trainsubz)'-1)';
dt22=dYears(t1(trainsubz),t1(trainsubz));
dy22=dDist(X1(trainsubz,1:2),X1(trainsubz,1:2));
fp22=bsxfun(@times,obsGISfp(trainsubz)-1,obsGISfp(trainsubz)'-1)';

if modno==6
    clear nnn;
    sample_Theta;
else
    [Ys0,~,~] = GaussianProcessRegression(meantime(trainsub0),Y(trainsub0),meantime(trainsubz),...
        modelspec.traincv(t1(trainsub0),t1(trainsub0),dt1t1(trainsub0,trainsub0),thetPXholo,Ycv(trainsub0,trainsub0),dy1y1(trainsub0,trainsub0),fp1fp1(trainsub0,trainsub0)),...
        modelspec.cvfunc(t1(trainsub0),meantime(trainsubz),dt12,thetPXholo,dy12,fp12)',...
        modelspec.cvfunc(meantime(trainsubz),meantime(trainsubz),dt22,thetPXholo,dy22,fp22),[]);
    first_flag=1;
    tic
    step_size=dY*2.5; 
    if modno==1
        step_size(suborb)=10e3;
    end
    if synthflag ==0
    else
        Ys0 = zeros(size(trainsubz,1),1);
    end
    tic
    sample_Ys_Thetas;    time_to_sample=toc;
end

train_all=find(limiting~=20);
defval('logp',1e-5);
logps=logp;
theta_samp=thet_samples;
%% redefining posteriors (only sampled data is used to train)
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

%% when normal is not sampled, then we have to add those in to be used in the regression
y_accepted=[];
if sampnorm==1
    y_accepted=thinned_ys(:,Nburn:end);
elseif modno~=6 && spatial==0 && modno~=7
    %% includes Ys for normally distributed and sampled values for sampled thinned_ys;
    y_int=repmat(Y,1,size(thinned_ys,2));
    y_int(trainsubz,:)=thinned_ys;
    y_accepted=y_int(train_all,nburn:end);
    trainsubx=train_all;
elseif modno==6 && spatial==1
    trainsubx=trainsub;
end
sub_all = intersect(find(limiting~=1),find(limiting~=-1));
Y0s=repmat(Y(train_all),1,size(y_accepted,2));
Y0s=y_accepted;
thinned_ts=[];
if nsamps > 200
    thets=thinned_thets(:,nburn:end);
else
    thets=thinned_thets;
end
