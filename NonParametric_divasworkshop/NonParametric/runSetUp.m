%% add the paths to run all of the scripts
clearvars -except distFile run_type Seed df pc pctext pd distKern0 distKernel distKernFlor jjj all run_start run_end counter CEFILES IFILES synthflag truth_flag Seed dateField rng datGIA sfile distNorm distKern df

if ~exist('run_type', 'var')==1
    run_type = 1; Seed = 4;
end

% if contains(pwd, 'scratch')
%     pd = '/scratch/ea289/Code/NP';
%     addpath([pd '/scripts']);
% else
%     pd = '/Users/ericaashe/Dropbox/Code/NP_publish/';
% end

addpath([pd '/MFILES']);
addpath(pd);

IFILES=[pd '/IFILES'];

if run_type == 1        % sensitivity tests
    spatial = 0; modno = 3; synthflag = 1; nsamps = 500; nburn = 100; thetachange = 29;
%    spatial = 0; modno = 3; synthflag = 1; nsamps = 20000; nburn = 100; thetachange = 29;
    stepchange = 51; nthin = 10; sampnorm = 0;
elseif run_type == 2    % validation, temporal only
    spatial = 0; modno = 3; synthflag = 2; nsamps = 1000; nburn = 10; thetachange = 15;
%    spatial = 0; modno = 3; synthflag = 2; nsamps = 500000; nburn = 100; thetachange = 15;
    stepchange = 12;  nthin = 10; sampnorm = 0;
else
    spatial = 0; modno = 3; synthflag = 0; nsamps = 1000; nburn = 10; thetachange = 15;
    stepchange = 12;  nthin = 10; sampnorm = 0;
end
 
% cd(pd);
% if synthflag == 1
%     dateField = [df  '_Sens' num2str(truth_flag)];
%     distFile = 'FL_RSLdata.csv';
% elseif synthflag == 2
%     distFile = 'Florida_data.csv';
%     dateField = [df '_Val'];
% else
%     dateField = df;
%     if ~exist("distFile",'var')
%         distFile = 'Florida_data.csv';
%     end
% end

cd(pd);
if synthflag == 1
    dateField = [df  '_Sens' num2str(truth_flag)];
    distFile = 'Subregion7_WA_RSL.csv';
elseif synthflag == 2
    distFile = 'Subregion7_WA_RSL.csv';
    dateField = [df '_Val'];
else
    dateField = df;
    if ~exist("distFile",'var')
        distFile = 'Subregion7_WA_RSL.csv';
    end
end



if exist('run_start', 'var')
        if run_start > 48
            load distNorm.mat;
            load distKern.mat;
        end
else
    load distKern.mat;
end

%distKern0=distKern;

datHolo=importdata(fullfile(IFILES,distFile));
lb=-300000; ub=200000; 

if exist('Seed', 'var')
    WORKDIR = ['Results_' dateField '_seed_' num2str(Seed)];
else
    WORKDIR = ['Results_' dateField];
end
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);
WORKDIR = pwd;

HoloRegions=unique(datHolo.data(:,15));
Holositelat=[]; Holositelong=[]; Holositename={};
ndat = size(datHolo.data,1); % + length(HoloRegions);

region = zeros(ndat,1);
datid = zeros(ndat,1);
time1 = zeros(ndat,1);
time2 = zeros(ndat,1);
dt = zeros(ndat,1);
limiting = zeros(ndat,1);
indic = zeros(ndat,1);
cspecies = zeros(ndat,1);
Y0 = zeros(ndat,1);
elev = zeros(ndat,1);
delev = zeros(ndat,1);
dY0 = zeros(ndat,1);
lat0 = zeros(ndat,1);
long0 = zeros(ndat,1);
dY_init = zeros(ndat,1);
meanYs = zeros(ndat,1);
meantime = zeros(ndat,1);

time1 = 1950-(datHolo.data(:,9)+datHolo.data(:,10)) + ndat/1e5;
time2 = 1950-(datHolo.data(:,9)-datHolo.data(:,11)) + ndat/1e5;
dt = abs(time1-time2)/4;
limiting = (datHolo.data(:,12));
indic = (datHolo.data(:,14));        
cspecies = (datHolo.data(:,16));
Y0 = datHolo.data(:,3)*1000;
elev = datHolo.data(:,6)*1000;
delev = datHolo.data(:,7)/4*1000; 
dY0 = datHolo.data(:,4)/2*1000;
lat0 = datHolo.data(:,1);
long0 = datHolo.data(:,2);
lat0 = lat0+1e-8*rand;

sitelen=[]; 
siteid=[]; 
sitenames={}; 
sitecoords=zeros(length(HoloRegions), 2); 

idHolo = 3e4;
for curreg = 1:length(HoloRegions)
    sub = find((datHolo.data(:,15)==HoloRegions(curreg)));
    count = [1:length(sub)]';
    datid(sub) = ones(length(sub),1)*(idHolo+HoloRegions(curreg)*1e3);
    region(sub) = ones(length(sub),1)*(HoloRegions(curreg));
    datid = [datid ; datid(sub(1))];
    region = [region ; region(sub(1))];
    time1 = [time1 ; 1925];
    time2 = [time2 ; 1975];
    dt = [dt ; 25]; 
    limiting = [limiting ; 0];
    indic = [indic ; 0 ];
    Y0 = [Y0 ; 0 ];
    dY0 = [dY0 ; 100];
    cspecies = [cspecies ; 0];
    lat0 = [lat0 ; lat0(sub(1)) ];
    long0 = [long0 ; long0(sub(1)) ];
    sitelen(curreg,1) = length(sub) + 1;
    sitecoords(curreg,:) = [mean(lat0(sub)) mean(long0(sub))];
    sitenames{curreg} = ['Holo-' datHolo.textdata{sub(1)+1,1}];
    siteid(curreg,1)=datid(sub(1));
end
        
dY_init = dY0;
meanYs = Y0;
meantime = mean([time1 time2],2);
dY = dY0;

suborb=intersect(find(limiting==0),find(cspecies==9));
trainsub0=find(limiting==0); %normally distributed
subpeat=intersect(find(limiting==0),find(cspecies==0));
trainsub2=find(limiting==2); %nonparametric
trainsub3=find(limiting==3); %uniform
trainsubn1=find(limiting==-1); %marine limiting (RSL above)
trainsub1=find(limiting==1); % terrestrial limiting (RSL below)
Y=Y0;

Y(trainsub2)=elev(trainsub2);
dY(trainsub2)=delev(trainsub2);
    
coral_params = [
    0.5356483	0.9061358	1.7     5.4     7.3     12.3
    1.6573675	0.9572511	5.2     17.5	24.3	41.9
    2.4202311	0.8501376	11.2	32.7	43.8	71.2
    2.146351	1.194733	8.5     38.3	57.8	38.3
    1.648682	1.238584	5.2     24.6	37.7	76.5
    2.3539188	0.8721877	10.5	31.5	42.5	69.9
    2.197737	1.232576	8.9     42.3	64.6	42.3
    2.5358296	0.9388956	12.5	41      56.7	96.9
    10.03894	4.612317	10      15.9	17.5	20.1
    2.043093	1.138957	7.7     32.2	47.7	91.4
    2.3446465	0.9594187	10.4	34.8	48.4	83.7 
    -0.292      0.149       0   0   0   0
    -0.3145     0.15975     0   0   0   0
    -0.317      0.161       0   0   0   0
    -0.414      0.208       0   0   0   0
    -0.4175     0.20925     0   0   0   0];

if spatial ==0
    lat = repmat(lat0(1), size(lat0));
    long = repmat(long0(1), size(long0));
    sitecoords(1,:) = [lat(1) long(1)];
else
    lat = lat0;
    long = long0;
end
if exist('savefile') %length(savefile)>0
    if exist('Y_fin')
    else
    end
end

if exist('run_start', 'var')
else
    run_start = 1; run_end = 2;
end  

sfile=['init_dat'];
save(sfile);

if exist('distNorm','var')
    distKern=distNorm;
end
