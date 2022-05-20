%% Distribution fitting for use in the runNonParametricGP script
%% must run prior to running that script, unless you have saved the distributions object
%carib=1; posSupport=1;% 1 for caribbean and 0 for florida

% specify the file that you are using to fit the distributions:
distFile = 'CaribbeanCorals.csv'; distFile2 = 'CaribbeanCorals2.csv';
%distFile = 'CaribbeanCorals.csv'; distFile2 = [];

% if you are using two files, specify the second file here 
%distFile2 = 'CaribbeanCorals2.csv';

% specify the date field
dateField = '050522';
% specify whether you are using positive support or not (0 or 1)
posSupport = 1;

if posSupport == 1
    bandWidth = .150;
else
    bandWidth = .800:.100:1.300;
end

% to change the directory you are using, substitute the appropriate folders
if ~exist('pd','var')
    pd = pwd;
end
cd(pd);
addpath([pd '/MFILES']);
ipd=pd;

% if (contains(pwd, '/Users/ericaashe'))
%     pd='/Users/ericaashe/Dropbox/';
% elseif (contains(pwd, '/projects/kopp/eashe/'))
%     pd='/projects/kopp/eashe/';
% end
% cd ([pd '/Code/NP']);
% addpath([pd '/Code/NP']);
% addpath([pd '/Code/NP/scripts']);
% addpath([pd '/Code/NP/MFILES']);
% ipd=[pd 'Code/NP/'];
IFILES=[ipd '/IFILES'];
WORKDIR=['Distributions' dateField]; 

if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end

cd(WORKDIR);
datCoral=importdata(fullfile(IFILES,distFile));

depths=zeros(length(datCoral.data(:,2)),1); taxa=depths; Cname={}; taxaVar=[];
coralTaxa=unique(datCoral.data(:,2));
for curreg=1:length(coralTaxa)
    sub=find(datCoral.data(:,2)==coralTaxa(curreg));
    count=[1:length(sub)]';
    wdepths=datCoral.data(sub,1);
    wtaxa=datCoral.data(sub,2);
    wCname=datCoral.textdata{sub(2),2};
    depths=[depths; wdepths];
    taxa=[taxa; wtaxa];
    Cname{coralTaxa(curreg)}=wCname;
    wvar = var(wdepths);
    taxaVar=[taxaVar; wvar];
end
   
if ~isempty(distFile2)
    datCoral2=importdata(fullfile(IFILES,distFile2));
    coralTaxa2=unique(datCoral2.data(:,2));
    for curreg = 1:length(coralTaxa2)
        %curreg=length(coralTaxa):length(coralTaxa)+length(coralTaxa2)
        sub=find(datCoral2.data(:,2)==coralTaxa2(curreg));
        count=[1:length(sub)]';
        wdepths=datCoral2.data(sub,1);
        wtaxa=datCoral2.data(sub,2);
        wCname=datCoral2.textdata{sub(2),2};
        depths=[depths; wdepths];
        taxa=[taxa; wtaxa];
        Cname{coralTaxa2(curreg)}=wCname;
        wvar = var(wdepths);
        taxaVar=[taxaVar; wvar];
    end
    cTaxa = [coralTaxa; coralTaxa2];
    Ntaxa=zeros(length(coralTaxa)+length(coralTaxa2),1);
else
    Ntaxa=zeros(length(coralTaxa),1);
    cTaxa = coralTaxa;
end

nlls=[];aic =[]; bic=[];
for jjj=1:length(Ntaxa)
    jj=cTaxa(jjj);
    sub=find(taxa==jj);
    Ntaxa(jjj) = length(sub);
    clf; clear title;
    edges=[0:1:min(max(depths(sub)),70)];
    h=histogram(depths(sub),edges,'normalization','probability'); hold on;
    depth_vals=min(depths(sub))-1:.200:min(max(depths(sub)),70);
    if posSupport == 1
        % this fitted distribution has a user-assigned band width
        distKern{jj}=fitdist(depths(sub)+.1,'Kernel','Support','positive','Width',bandWidth);
        % this distribution has an optimized bandwidth
        distKern0{jj}=fitdist(depths(sub)+.1,'Kernel','Support','positive');
        title(Cname{jj}(1:min(length(Cname{jj}),11)));
        pdfwrite(['Fitted_Kernels_Pos_only_hist' Cname{jj}]);
        plot(depth_vals,pdf(distKern{jj},depth_vals')); 
        nlls(jj)=negloglik(distKern{jj});
        [aic(jjj) bic(jjj)]=aicbic(-nlls(jj),1,length(sub));
        if max(depths(sub))>70
            xlim=([-5 70]);
        end
        legend(['BW ' num2str(distKern{jj}.BandWidth)]);
        %legend('BW 0.15');
        clf;clear title;
        %figure;
        h=histogram(depths(sub),edges,'normalization','probability'); hold on;
        plot(depth_vals,pdf(distKern0{jj},depth_vals')); 
        if max(depths(sub))>70
            xlim=([-5 70]);
        end
        legend(['BW ' num2str(distKern0{jj}.BandWidth)]);
        title(Cname{jj});
        pdfwrite(['Fitted_Kernels_Pos_only_BW' Cname{jj}]);
    else
        [f,xi]=ksdensity(depths(sub),'bandwidth',bandWidth); %,'function','pdf');
        plot(xi,f,'LineWidth',2);
        legend('BW 0.15');
        pdfwrite(['Fitted_Kernels_hist' Cname{jj}]);
    end
    if max(depths(sub))>70
        xlim=([-5 70]);
    end
end

distVar=[];
for ii=1:length(distKern)
     distVar=[distVar; var(distKern{ii})];
end

cd(pd);
save('distKern.mat', 'distKern', 'Cname', 'distVar');