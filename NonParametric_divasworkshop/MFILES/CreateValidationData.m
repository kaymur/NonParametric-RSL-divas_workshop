% Synthetic data for validations 
%%                      meas_unc    temp_unc
%%   all non-normal     0.85        310.44
%%   a. palmata         0.72        233.68
%%   sed                0.40        291.78
    
% 1.  Perturb true times (from original data)
% 2.  Perturb true RSLs (from GIA model interpolation)
% 3.  Use depth distributions to create elevations from perturbed RSLs

N=length(Y); 

GIALat=unique(datGIA(:,1));
GIALon=mod(unique(datGIA(:,2)),360);
GIAs = 1000.*datGIA(:,3:end);
GIAtime = 1950+1000.*[-17:1:0]';
GIA0 = reshape(GIAs,21,20,[]);
%GIA0 = reshape(GIAs,21,23,[]);

%% subtract GIA for each lat,long according to where the data was
for jj = 1:length(Y)
        scoord(1) = lat(jj);                
        scoord(2) = mod(long(jj),360);
        age(jj,1) = meantime(jj);
        trueYs(jj,1)=interp3(GIALat,GIALon,GIAtime,GIA0,scoord(1),scoord(2),age(jj));
end

age(1) = 1950; trueYs(1) = 0; limiting(1) = 0; cspecies(1) = 0;
ages = zeros(N,1);
Y_fin = ages;
Y_synth = trueYs;
dt=ages;
t_err = ages;
meas_err = ages;
offsets = ages;
%% assign temporal error
%%                      meas_unc (m)   temp_unc (yrs)
%%   all non-normal     0.85            310.44
%%   a. palmata         0.72            233.68
%%   sed                0.40            291.78
 
%% draw offset and errors from the empirical distributions and normal measurement errors

%% sedimentary data
sed=intersect(find(cspecies==0), find(limiting == 0)); ns=length(sed);
rng(Seed*2);
dt_s= max(250+100*randn([ns,1]),23);
dt_s(1)=10;
dt(sed)=dt_s;
rng(Seed*8);
temp_err_s = dt_s.*randn([ns,1]);
t_err(sed) = temp_err_s;
ages(sed)=age(sed)+temp_err_s;
rng(Seed*3);
dY_s=max(400+100*randn([ns,1]),15);
dY_s(1) = 12;
dY(sed) = dY_s;
rng(Seed*4);
meas_err_s = dY_s.*randn([ns,1]);
Y_synth(sed) = trueYs(sed) + meas_err_s;
Y_fin = Y_synth;
meas_err(sed) = meas_err_s;
rng(Seed*5);
    offset_s=1500*randn(ns,1); 
    offset_s(1) = 0;
    Y_fin(sed)=Y_synth(sed)+offset_s;
    limiting(sed)=0;
    cspecies(sed)=0;
offsets(sed) = offset_s;

%% orbicella
orb=find(cspecies==9); no=length(orb);
rng(Seed*6);
dt_o= max(375+100*randn([no,1]),100);
dt(orb)=dt_o;
rng(Seed*5);
temp_err_o = dt_o.*randn([no,1]);
t_err(orb) = temp_err_o;
ages_o=age(orb)+temp_err_o;
ages(orb) = ages_o;
rng(Seed*9);
dY_o=max(850+100*randn([no,1]),50);
rng(Seed*2);
meas_err_o = dY_o/2.*randn([no,1]);
meas_err(orb) = meas_err_o;
Y_synth(orb) = trueYs(orb) + meas_err_o;
rng(Seed*8);
offset_o=-1000*random(distKern{9},no,1); 
Y_fin(orb)=Y_synth(orb)+offset_o;
limiting(orb)=2;
cspecies(orb)=9;
offsets(orb) = offset_o;
    
%% Acropora Palmata    
acr=find(cspecies==1); na=length(acr);
rng(Seed*7);
dt_a= max(225+100*randn([na,1]),10);
dt(acr)=dt_a;
rng(Seed*1);
temp_err_a = dt_a.*randn([na,1]);
t_err(acr) = temp_err_a;
ages(acr)=age(acr)+temp_err_a;
rng(Seed*6);
dY_a=max(700+100*randn([na,1]),50);
dY(acr) = dY_a;
rng(Seed*3);
meas_err_a = dY_a.*randn([na,1]);
Y_synth(acr) = trueYs(acr) + meas_err_a;
meas_err(acr) = meas_err_a;
rng(Seed*4);
offset_a=-1000*random(distKern{1},na,1); 
Y_fin(acr)=Y_synth(acr)+offset_a;
limiting(acr)=2;
cspecies(acr)=1;
offsets(acr) = offset_a;

%% Limiting Data
sublims = union(find(limiting==1),find(limiting==-1));
nl = length(sublims);
rng(Seed*11);
dt_l= max(373+100*randn([nl,1]),10);
dt(sublims)=dt_l;
rng(Seed*5);
temp_err_l = dt_l.*randn([nl,1]);
t_err(sublims) = temp_err_l;
ages(sublims)=age(sublims)+temp_err_l;
rng(Seed*9);
dY_l=max(850+100*randn([nl,1]),50);
dY(sublims) = dY_l;
rng(Seed*2);
meas_err_l = dY_l.*randn([nl,1]);
Y_synth(sublims) = trueYs(sublims) + meas_err_l;
meas_err(sublims) = meas_err_l;
rng(Seed*8);
offset_l=70000*rand(nl,1)-35000;
Y_fin(sublims)=Y_synth(sublims)+offset_l;
b1=find(offset_l<0);
b2=find(offset_l>=0);
limiting(sublims(b1))=-1;
limiting(sublims(b2))=1;
cspecies(sublims(b1))=0;
cspecies(sublims(b2))=0;

oth=intersect(find(cspecies>1),find(cspecies~=9));
limiting(oth)=20;
cspecies(oth)=20;

cspec=cspecies';
Y=Y_fin;
