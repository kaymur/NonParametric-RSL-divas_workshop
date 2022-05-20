
maxstep=max(8*dY(trainsubz),1000);
minstep=100;
tri_max = 50e3;       % set the limit of the triangle distribution
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

% Last updated by Erica Ashe, May 19 2022

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
logp(1,:) = logp0*8e4;

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
