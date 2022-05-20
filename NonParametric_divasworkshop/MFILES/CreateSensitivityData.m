function [ages,dt,age,Y_fin,Y_synth,Y_syn,dY,limiting,cspec,meas_err,offset]=CreateSensitivityData(t_uncert,n,age_range,truth_flag,distKernel,limiting,cspecies,dat_type,jj,df,Seed)
% synthetic data creation for sensitivity tests
%%                      meas_unc    temp_unc
%%   all non-normal     0.85        310.44
%%   a. palmata         0.72        233.68
%%   sed                0.40        291.78
rng(Seed+jj);
N=n*12;
limiting = limiting(1:N);
cspecies = cspecies(1:N);
t_unc = repelem(t_uncert,N)';
%% data not evenly spaced, so as to account for missingness
%% randomly assign dates within the age range
a2=age_range(2);
a1=age_range(1);
age=[];     % true age
ages=[];    % perturbed ages
for ii=1:2000:a2
    wage = unifrnd(ii-1,ii+1999,[n*2,1]);
    age=[age wage'];
end
age=sort(age);
%% assign temporal error (temporal error in database is centered around 200 years with a max of 852 & min of 22
rng(Seed+2*jj);
dt=max(t_unc+t_unc*0.2.*randn([N,1]),0);
dt(1)=10;
age(1)=0;
rng(Seed);
temp_err = dt.*randn([N,1]);
ages=age'+temp_err;

%% one function with similar rates to Holocene
if truth_flag==1
    Y_syn(:,1) = -20000*cos((age+6280)/2000)-20000;
elseif truth_flag==2
%% other function, similar to post-glacial (18-8ka)
        sub0=find(age>7000);
        sub1=find(age<=7000 & age>6000);
        sub2=find(age<=6000 & age>5200);
        sub3=find(age<=5200 & age>3800);
        sub4=find(age<=3800 & age>3000);
        sub5=find(age<=3000);
        
        Y0=age*0; Y1=Y0; Y2=Y0; Y3=Y0; Y4=Y0; Y5=Y0;
        
        Y5(sub5)=(age(sub5))*-20;
        Y4(sub4)=-60e3+(age(sub4)-3000)*-2.5;
        Y3(sub3)=-62e3+(age(sub3)-3800)*-18;
        Y2(sub2)=-87200+(age(sub2)-5200)*-1.2;
        Y1(sub1)=-88160+(age(sub1)-6000)*-40;
        Y0(sub0)=-128160+(age(sub0)-7000)*-5;
        
        Y_syn(:,1) =Y0+Y1+Y2+Y3+Y4+Y5;
end

rng(Seed+jj);
dY=max(400 + 100*randn([N,1]),15);
dY(1)=50;
rng(Seed+jj);
meas_err = dY.*randn([N,1]);
Y_synth = Y_syn + meas_err;

%% draw offset from the empirical distributions
%% dat_type defines whether we are using the sedimentary (3), acropora(1), or orbicella(2) distribution
Y_fin=Y_synth;
if dat_type==1                  % a. palmata
    offset=-1000*random(distKernel{1},N,1); 
    Y_fin(2:end)=Y_synth(2:end)+offset(2:end);
    limiting(2:end)=2;
    cspecies(2:end)=1;
%     dY_a=max(250*randn([N-1,1]),250);
%     dY(2:end) = dY_a;
elseif dat_type==2              % orbicella
rng(Seed+3*jj);
    offset=-1000*random(distKernel{9},N,1); 
    Y_fin(2:end)=Y_synth(2:end)+offset(2:end);
    limiting(2:end)=2;
    cspecies(2:end)=9;
    rng(Seed*2);
%     dY_o=max(750*randn([N-1,1]),500);
%     dY(2:end) = dY_o;
elseif dat_type==3              % sedimentary
rng(Seed+4*jj);
    offset=1500*randn(N,1);
%     dY_s=max(200*randn([N-1,1]),100);
%     dY_s(1) = 12;
%     dY(2:end) = dY_s;
    Y_fin(2:end)=Y_synth(2:end)+offset(2:end);
    limiting(2:end)=0;
    cspecies(2:end)=0;
elseif dat_type==4              % normal data & A. palmata
% select half of the data to be normal and half coral
    rng(Seed+5*jj);
    all=find(limiting>-20);
    limiting(all)=20;
    nn=N/2;
    rng(Seed+6*jj);
    norms = round(unifrnd(1,N,[nn-1,1]),0);
    norms=[1; norms];
    norms = unique(norms);
    while length(norms)<nn
        norms = [norms; round(unifrnd(1,N,[1,1]),0)];
        norms = unique(norms);
    end
    rng(Seed+7*jj);
    offset=60*randn(N/2,1);
%     dY_s=max(200*randn([N/2,1]),100);
%     dY_s(1) = 12;
%     dY(norms) = dY_s;
%     dY(norms)=(dY(norms).^2+60^2).^(0.5);
    Y_fin(norms)=Y_synth(norms)+offset;
    limiting(norms)=0;
    cspecies(norms)=0;
    oth=find(limiting==20);
    rng(Seed+8*jj);
    offset=-1000*random(distKernel{1},N/2,1);
%     dY_o=max(750*randn([N/2,1]),500);
%     dY(oth)=dY_o;
    Y_fin(oth)=Y_synth(oth)+offset;
    limiting(oth)=2;
    cspecies(oth)=1;
elseif dat_type==5
    all=find(limiting>-20);
    limiting(all)=20;
    nn=N/2;
    rng(Seed+9*jj);
    norms = round(unifrnd(1,N,[nn-1,1]),0);
    norms=[1; norms];
    norms = unique(norms);
    while length(norms)<nn
        norms = [norms; round(unifrnd(1,N,[1,1]),0)];
        norms = unique(norms);
    end
    rng(Seed+10*jj);
    offset=60*randn(N/2,1);
%     dY_s=max(150*randn([nn,1]),15);
%     dY_s(1) = 12;
%     dY(norms) = dY_s;
%     dY(norms)=(dY(norms).^2+60^2).^(0.5);
    Y_fin(norms)=Y_synth(norms)+offset;
    limiting(norms)=0;
    cspecies(norms)=0;
    oth=find(limiting==20);
    rng(Seed+11*jj);
    offset=-1000*random(distKernel{9},N/2,1); 
%     dY_o=max(750*randn([N/2,1]),500);
%     dY(oth)=dY_o;
    Y_fin(oth)=Y_synth(oth)+offset;
    limiting(oth)=2;
    cspecies(oth)=9;
elseif dat_type==6                  % Limiting data
    offset=70000*rand(N,1)-35000;
%    offset=70000*rand(N,1)-35000;
    Y_fin(2:end)=Y_synth(2:end)+offset(2:end);
%     dY_o=max(750*randn([N-1,1]),500);
%     dY(2:end)=dY_o;
    b1=find(offset<0);
    b2=find(offset>=0);
    limiting(b1)=-1;
    limiting(b2)=1;
    limiting(1)=0;
    cspecies(b1)=0;
    cspecies(b2)=0;
    
elseif dat_type==7                  % limiting & orbicella
%% select half of the data to be orbicella and half limiting
    all=find(limiting>-20);
    limiting(all)=20;
    nn=N/2;
    rng(Seed+12*jj);
    cor = round(unifrnd(1,N,[nn,1]),0);
    cor = unique(cor);
    cor = cor(cor~=1);
    while length(cor)<nn
        cor = [cor; round(unifrnd(1,N,[1,1]),0)];
        cor = unique(cor);
        cor = cor(cor~=1);
    end
    rng(Seed+13*jj);
    offset=-1000*random(distKernel{9},N/2,1); 
    Y_fin(cor)=Y_synth(cor)+offset;
    dY_c=max(750*randn([N/2,1]),500);
    dY(cor)=dY_c;
    limiting(cor)=2;
    cspecies(cor)=9;
%     Y_fin(cor)=Y_synth(cor)+offset;
%     limiting(cor)=2;
%     cspecies(cor)=9;
    oth=find(limiting==20);
%     dY(oth) = max(1000*randn([N/2,1]),650);
    rng(Seed+14*jj);
dY_l=max(850+100*randn([N/2,1]),500);
dY(oth) = dY_l;
rng(Seed*2);
meas_err(oth) = dY_l.*randn([N/2,1]);
Y_synth(oth) = Y_synth(oth) + meas_err(oth);
    offset=70000*rand(N/2,1)-35000;
    Y_fin(oth)=Y_synth(oth)+offset;
    b1 = find(offset<0);b0=find(limiting==20);b2=find(offset>=0);
    bl=b0(b1);bu=b0(b2);
    limiting(bl)=-1;
    limiting(bu)=1;
    cspecies(oth)=0;
elseif dat_type==8                  % limiting & normal
%% select half of the data to be normal and limiting
    all=find(limiting>-20);
    limiting(all)=20;
    nn=N/2;
    rng(Seed+15*jj);
    norms = round(unifrnd(1,N,[nn-1,1]),0);
    norms=[1; norms];
    norms = unique(norms);
    while length(norms)<nn
        norms = [norms; round(unifrnd(1,N,[1,1]),0)];
        norms = unique(norms);
    end
    rng(Seed+16*jj);
    offset=60*randn(N/2,1);
%     dY_s=max(200*randn([nn,1]),100);
%     dY_s(1) = 12;
%     dY(norms) = dY_s;
%     dY(norms)=(dY(norms).^2+60^2).^(0.5);
    Y_fin(norms)=Y_synth(norms)+offset;
    limiting(norms)=0;
    cspecies(norms)=0;
    oth=find(limiting==20);
    rng(Seed+17*jj);
    offset=70000*rand(round(N/2),1)-35000;
    Y_fin(oth)=Y_synth(oth)+offset;
%     dY_o=max(1000*randn([N/2,1]),650);
%     dY(oth)=dY_o;
    b1 = find(offset<0);b0=find(limiting==20);b2=find(offset>=0);
    bl=b0(b1);bu=b0(b2);
    limiting(bl)=-1;
    limiting(bu)=1;
    cspecies(oth)=0;
end
Y_fin(1)=0;
limiting(1)=0;
cspecies(1)=0;
clf;
plot(ages,Y_fin,'.'); hold on;
plot(ages,Y_fin+2*dY,'.');
plot(ages,Y_fin-2*dY,'.');
if truth_flag==2
    agex=[0:100:13000]';
    sub0=find(agex>7000);
    sub1=find(agex<=7000 & agex>6000);
    sub2=find(agex<=6000 & agex>5200);
    sub3=find(agex<=5200 & agex>3800);
    sub4=find(agex<=3800 & agex>3000);
    sub5=find(agex<=3000);
    
    Y0=agex*0; Y1=Y0; Y2=Y0; Y3=Y0; Y4=Y0; Y5=Y0;
        
    Y5(sub5)=(agex(sub5))*-20;
    Y4(sub4)=-60e3+(agex(sub4)-3000)*-2.5;
    Y3(sub3)=-62e3+(agex(sub3)-3800)*-18;
    Y2(sub2)=-87200+(agex(sub2)-5200)*-1.2;
    Y1(sub1)=-88160+(agex(sub1)-6000)*-40;
    Y0(sub0)=-128160+(agex(sub0)-7000)*-5;
        
        
    truths =Y0+Y1+Y2+Y3+Y4+Y5;
    plot(agex,truths);
elseif truth_flag==1
    agesamp=[0:100:13000];
    truths=-20000*cos((agesamp+6280)/2000)-20000;
    plot(agesamp,truths);
end
if min(ages)<0
    mage = find(ages<0);
    ages(mage) = 0;
end

plot(ages,Y_fin,'*');
if dat_type ==7
       plot(ages(oth),Y_fin(oth),'r*');
       plot(ages(cor),Y_fin(cor),'k*');
end       

cspec=cspecies';
pdfwrite(['True_plus_Corrup_' df '_' num2str(jj)]);

end

