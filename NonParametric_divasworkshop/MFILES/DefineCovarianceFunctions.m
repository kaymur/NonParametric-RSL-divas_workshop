% Define Covariance structures for model
refyear=1950;

dYears = @(years1,years2) abs(bsxfun(@minus,years1',years2));
dYears0 = @(years1,years2) (bsxfun(@minus,years1',years2));
angd = @(Lat0,Long0,lat,long) (180/pi)*(atan2( sqrt( (cosd(lat) .* sind(long-Long0)).^2 + (cosd(Lat0) .* sind(lat) - sind(Lat0) .* cosd(lat) .* cosd(long-Long0)).^2),(sind(Lat0) .* sind(lat) + cosd(Lat0) .* cosd(lat) .* cosd(long-Long0))));
dDist = @(x1,x2) angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))' + 1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);
kMat1 = @(dx,thetas) thetas(1).^2 .* (1).*exp(-dx/thetas(2));
kMat3 = @(dx,thetas) thetas(1).^2 .* (1 + sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
kMat5 = @(dx,thetas) thetas(1).^2 .* (1 + (sqrt(5)*dx/thetas(2)).*(1 + sqrt(5)*dx/thetas(2)/3)).*exp(-sqrt(5)*dx/thetas(2));
kSE = @(dx,thetas) thetas(1).^2 * exp(-(dx.^2)/(2*thetas(2).^2));
kDELTA = @(dx,thetas) thetas(1).^2 .* (dx==0);
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));
kCONST = @(thetas) thetas(1).^2;
kSIN = @(dt,thetas) thetas(1).^2 * exp(-2*sin(pi*dt/(thetas(2))).^2/thetas(3).^2);
kMODSIN = @(dt,thetas) kSIN(dt,thetas(1:3)) .* kSE(dt,[1 thetas(4)*thetas(2)]);
kFREQ = @(dt,thetas) thetas(1).^2 * cos(2*pi*dt/thetas(2));
kRQ = @(dx,thetas) thetas(1).^2 * (1 + dx.^2/(2*thetas(2)*thetas(3))).^-thetas(3);
kMatG = @(dx,thetas) thetas(1).^2 .* 2.^(1-thetas(3))./gamma(thetas(3)) .* (sqrt(2*thetas(3))*(dx+eps)/thetas(2)).^thetas(3) .* besselk(thetas(3),sqrt(2*thetas(3))*(dx+eps)/thetas(2));
kDELTAG = @(ad,thetas)thetas(1).^2.*(abs(ad)<1e-4).*(ad<360);
kGEOG = @(ad,thetas) kMat5(ad,thetas) .* (ad<360);
kGEOGG = @(ad,thetas)kMatG(ad,thetas).*(ad<360);
kMat3d = @(years1,years2,dx,thetas)thetas(1).^2.*(-3/(thetas(2).^2)).*dx.*exp(-sqrt(3)*dx/thetas(2)).*(-1+2*bsxfun(@ge,years1',years2));
kDPd = @(years1,years2,thetas) thetas(1).^2 * repmat((years1-refyear)',length(years2),1);
kMat3dd = @(dx,thetas)thetas(1).^2.*(3/(thetas(2).^2)).*(1-sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
kDPdd = @(years1,years2,thetas) thetas(1).^2 * ones(length(years2),length(years1));
FiniteMask = @(x1,x2) (bsxfun(@plus,sum(abs(x1),2)',sum(abs(x2),2)))<1e12;


% defines covariance functions
kMat5d = @(years1,years2,dx,thetas)thetas(1).^2.*(-5*dx.*exp(-sqrt(5)*dx/thetas(2))).*(thetas(2)+sqrt(5)*dx)/(3*thetas(2).^3);
kMat5dd = @(dx,thetas)thetas(1).^2.*-(5*exp(-(sqrt(5)*dx)/thetas(2)).*(thetas(2).^2+sqrt(5)*thetas(2).*dx-5*dx.^2))/(3*thetas(2).^4);
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));

%global covariance function
cvfunc.G = @(dt1t2,thetas) kMat3(dt1t2,thetas(1:2));
cvfunc.H = @(dt1t2,thetas) kMat3(dt1t2,thetas(1:2));
cvfunc.W = @(dt1t2,ad,thetas) kDELTA(dt1t2,thetas(1));% .* kDELTAG(ad,1);

%defines temporal derivative of covariance function
dcvfunc.G = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2));
dcvfunc.H = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2));
dcvfunc.W = @(dt1t2,ad,thetas) 0;

%second derivative
ddcvfunc.G = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2));
ddcvfunc.H = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2));
ddcvfunc.W = @(dt1t2,ad,thetas) 0;

clear modelspec;

%define model specifications
modelspec(1).label = 'Full';

%model is sum of functions we defined above
modelspec(1).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.W(dt1t2,ad,thetas(5)) + cvfunc.H(dt1t2,thetas(3:4)) ;
modelspec(1).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) +dcvfunc.H(t1,t2,dt1t2,thetas(3:4));% + dcvfunc.R(t1,t2,dt1t2,ad,thetas(3:5));
modelspec(1).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) +ddcvfunc.H(dt1t2,thetas(3:4));% + ddcvfunc.R(dt1t2,ad,thetas(3:5));

%adds defined error covariance to the covariance function
modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

%%%%%
%define starting values and lower and upper bounds of the hyperparameters
%%%%%
tluTGG = [
10e3 5e3 50e3      % global amplitude (starting, lower, and upper); SL in mm; how much values can go up or down over a particular time period  
12e3 4e3 40e3    % temporal scale in years, how much it can very over a particular time scale
1000 1 20e3      % high frequency amplitude (starting, lower, and uppper); SL in mm; how much values can go up or down over a particular time period  
1e3 500 10e3    % HF temporal scale in years, how much it can very over a particular time scale
1e3 1 10e3     % white noise    
];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';
modelspec(1).subfixed=[ ];
theta0=modelspec(1).thet0;
%modelspec(1).subfixed=[ 6 7 8]; %which parameters we're not sampling

%we won't do anything with this now
modelspec(1).sublength=[]; % these are ones that will be tuned based only on tide gauges if use Optimize level 1
modelspec(1).subamp = [1 3 ];
modelspec(1).subampnoise = [3 ];

