function MakePlotsSampling(dataset,f2s_shifted,sd2s,V2s,testloc,labl,doplots,difftimestep,xlim0,maxdistfrom,meanSL,testt,coral_params)

% Last updated by Erica Ashe, April 2017
% plotting all sites on one figure

defval('testloc',0);
defval('testt',0);
defval('labl','');
defval('doplots',[]);
defval('difftimestep',100);
defval('xlim0',[0 12000]);
defval('maxdistfrom',.1);
numrows=1+(difftimestep>0);
if min(dataset.meantime) < -3000
    dataset.meantime = 1950 - dataset.meantime;
    dataset.time1 = 1950 - dataset.time1;
    dataset.time2 = 1950 - dataset.time2;
end
%indic=dataset.indic;

%% set colors for the data plots
% colors=[ 1 .55 0.05
%         1 .85 .1
%         .68 .20 .70
%         0 0 0
%         0 1 0
%         1 .22 .70
%         0.2 0.27 .9
%         %0.90 0.12 0.13
%         1 0 0
%         0 .96 .86
%         0 .45 .1
%         .18 .09 .34
%         ];
colors=[ 1 .55 0.05
        1 0 0
        .68 .20 .70
        0 0 0
        0 1 0
        1 .22 .70
        0.2 0.27 .9
        %0.90 0.12 0.13
        1 .85 .1
        0 .96 .86
        0 .45 .1
        .18 .09 .34
        ];
        
angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

if ~isstruct(testloc)
    if testt==0
        testlocs.X=dataset.meantime;
        testlocs.reg=dataset.datid(1)*ones(length(dataset.meantime),1);
        testlocs.fp=ones(length(dataset.meantime),1);
        testlocs.GIAproju=0;
        testlocs.GIAproj=zeros(length(dataset.meantime),1);
    else
        testlocs.X=testt;
        testlocs.reg=dataset.datid(1)*ones(length(testt),1);
        testlocs.fp=ones(length(testt),1);
        testlocs.GIAproju=0;
        testlocs.GIAproj=zeros(length(testt),1);
    end
    testreg=testlocs.reg;
    testsites=[dataset.siteid mean(dataset.lat) mean(dataset.long)];
    testnames=dataset.sitenames;
    testnames2=dataset.sitenames;
    testX(:,1:3)=[mean(dataset.lat)*ones(length(testt),1) mean(dataset.long)*ones(length(testt),1) testt];
    %testX(:,3)=1950-dataset.meantime;
    %testX(:,3)=1950-dataset.meantime;
else
testreg=testloc.reg;
testsites=testloc.sites;
testnames=testloc.names;
testnames2=testloc.names2;
testX(:,1:2)=testloc.X(:,1:2);
testX(:,3)=testloc.X(:,3);
end
% testreg=testlocs{1}.reg;
% testsites=testlocs{1}.sites;
% testnames=testlocs{1}.names;
% testnames2=testlocs{1}.names2;
% testX=testlocs{1}.X;

% change everything back to ages instead of common era dates
meantime=dataset.meantime;
%meantime=dataset.meantime;
lat=dataset.lat;
long=dataset.long;
Y=dataset.Y+meanSL;
Ycv=dataset.Ycv;
limiting=dataset.limiting;
%indicator=dataset.indic;
% intercalated=dataset.intercalated;
% compactcorr=dataset.compactcorr;
cspecies=dataset.cspecies;
time1=dataset.time1;
time2=dataset.time2;
%time1=dataset.time1;
%time2=dataset.time2;
dY=dataset.dY;
datid=dataset.datid;
%dY_init=dataset.dY_init;
% shift all sea levels back by the original mean of the data
f2s=f2s_shifted+meanSL;
offsetA=0;
%figure;
% 
% for i=1:size(testsites,1) % = 24 sites
     wxlim=[];
% 	figure;
    
%     hp(i)=subplot(numrows,length(js),k+length(js));
%     box on;
% 
     hold on

%             col{i} = colors(10*i,:);
%             subB=find(datid==testsites(i,1));
            subA = find(testreg == testsites(1,1));        
%             clf; clear hp;	
                js=1;        j=js;
                k=1;
%		hp(k)=subplot(numrows,length(js),k);
         
%         [tX Xis]=sort(testX(:,3));
%         plotdat.x=(tX);
%         plotdat.y=f2s(Xis);
%         plotdat.dy=[f2s(Xis)-2*sd2s(Xis)];
% 
%         PlotWithShadedErrors(plotdat,[0 0 0],[.9 .9 .9]);
%         hold on;
% 
%         plot(testX(:,3),f2s(subA)+offsetA,'Color',[.5, .5, .5]);
%         hold on;
% 		plot(testX(subA,3),f2s(subA)+sd2s(subA)+offsetA,'--','Color',[.5, .5, .5]);
% 		plot(testX(subA,3),f2s(subA)-sd2s(subA)+offsetA,'--','Color',[.5, .5, .5]);
% 		plot(testX(subA,3),f2s(subA)+2*sd2s(subA)+offsetA,':','Color',[.5, .5, .5]);
% 		plot(testX(subA,3),f2s(subA)-2*sd2s(subA)+offsetA,':','Color',[.5, .5, .5]);
means = coral_params(:,3)';
    for nn = 1:length(datid)
        if limiting(nn)==0&&cspecies(nn)==0
            %indic(nn)==1
            rectangle('Position',[min(time1(nn),time2(nn)),Y(nn)-2*dY(nn),dYears(time2(nn),time1(nn)),4*dY(nn)],'Edgecolor',[0,0,0]);  % peat index points basal, black
        elseif limiting(nn)==0&&cspecies(nn)>=10
            rectangle('Position',[min(time1(nn),time2(nn)),Y(nn)-2*dY(nn),dYears(time2(nn),time1(nn)),4*dY(nn)],'Edgecolor',[0,0,0]);  % peat index points basal, black
        elseif limiting(nn)==0&&cspecies(nn)==9
                rectangle('Position',[min(time1(nn),time2(nn)),Y(nn)+coral_params(cspecies(nn),1)*1000-2000*coral_params(cspecies(nn),2),...
                dYears(time2(nn),time1(nn)),4000*coral_params(cspecies(nn),2)],'Edgecolor',colors(cspecies(nn),:));  % peat index points basal, turquoise
%             else
%                 plot([1 1]*meantime(nn),[Y(nn)-2*dY(nn) Y(nn)+4*means(cspecies(nn))*1000+2*dY(nn)],'Color',colors(cspecies(nn),:));
%                 plot([time1(nn) time2(nn)],[1 1]*Y(nn),'Color',colors(cspecies(nn),:));
        elseif cspecies(nn)==9
                plot([1 1]*meantime(nn),[Y(nn)-2*dY(nn) Y(nn)+2*means(cspecies(nn))*1000+2*dY(nn)],'Color',colors(cspecies(nn),:));
                plot([time1(nn) time2(nn)],[1 1]*Y(nn),'Color',colors(cspecies(nn),:));
        elseif limiting(nn)==0&&cspecies(nn)==1
            rectangle('Position',[min(time1(nn),time2(nn)),Y(nn)+coral_params(cspecies(nn),1)*1000-2000*coral_params(cspecies(nn),2),...
                dYears(time2(nn),time1(nn)),4000*coral_params(cspecies(nn),2)],'Edgecolor',colors(cspecies(nn),:));  % peat index points basal, turquoise
        elseif limiting(nn)==2
            plot([1 1]*meantime(nn),[Y(nn)-2*dY(nn) Y(nn)+means(cspecies(nn))*1000+2*dY(nn)],'Color',colors(cspecies(nn),:));
            plot([time1(nn) time2(nn)],[1 1]*Y(nn),'Color',colors(cspecies(nn),:));
        elseif limiting(nn)==3 % uniform
            rectangle('Position',[min(time1(nn),time2(nn)),Y(nn)-2*dY(nn),dYears(time2(nn),time1(nn)),4*dY(nn)],'Edgecolor',[.8,.8,.8]);  % uniform, grey
            hold on
        end
    end
    for nn = 1:length(datid)
        if limiting(nn) == -1 % marine limiting
            plot([1 1]*meantime(nn),Y(nn)+[-2 2]*dY(nn),'Color',[.1,.35,.7]);
            plot([time1(nn) time2(nn)],[1 1]*(Y(nn)+2*dY(nn)),'Color',[.1,.35,.7]);
        elseif limiting(nn)==1 % terrestrial limiting
            plot([1 1]*meantime(nn),Y(nn)+[-2 2]*dY(nn),'Color',[.05,.7,.5]);
            plot([time1(nn) time2(nn)],[1 1]*(Y(nn)-2*dY(nn)),'Color',[.05,.7,.5]); 
        end
    end
%     for ii = 1:length(sub_none)
%         nn=sub_none(ii);
%         plot([1 1]*meantime(nn),[Y(nn)-2*dY(nn) Y(nn)+means(cspecies(nn))],'r');
%         plot([time1(nn) time2(nn)],[1 1]*Y(nn),'r');
%     end
        
%         if cspecies(nn)==0||cspecies(nn)==4||cspecies(nn)==2 %
%             if limiting(nn) == -1 % marine limiting
%                 plot([1 1]*meantime(nn),Y(nn)+[-2 2]*dY(nn),'Color',[.1,.35,.7]);
%                 plot([time1(nn) time2(nn)],[1 1]*(Y(nn)+2*dY(nn)),'Color',[.1,.35,.7]);
%             elseif limiting(nn)==1 % terrestrial limiting
%                 plot([1 1]*meantime(nn),Y(nn)+[-2 2]*dY(nn),'Color',[.05,.7,.5]);
%                 plot([time1(nn) time2(nn)],[1 1]*(Y(nn)-2*dY(nn)),'Color',[.05,.7,.5]); 
%             elseif limiting(nn)==0||limiting(nn)==2
%                 rectangle('Position',[min(time1(nn),time2(nn)),Y(nn)-2*dY(nn),dYears(time2(nn),time1(nn)),4*dY(nn)],'Edgecolor',[.42,.77,1]);  % peat index points basal, turquoise
%             end
%         elseif cspecies(nn)==1||cspecies(nn)==3 % log-normal
%             %% deciding how to plot the log-normal
%             plot([1 1]*meantime(nn),[Y(nn)+2*dY(nn) Y(nn)-means(cspecies(nn))],'red');
%             plot([time1(nn) time2(nn)],[1 1]*Y(nn),'red');
%             %plot([1 1]*meantime(nn),[Y(nn)-[0 means(cspecies(nn))]*dY(nn)],'red');
%         elseif cspecies(nn)==5 % uniform
%             rectangle('Position',[min(time1(nn),time2(nn)),Y(nn)-2*dY(nn),dYears(time2(nn),time1(nn)),4*dY(nn)],'Edgecolor','y');  % peat index points basal, turquoise
%         end
%         ps{i} = plot(nan,nan,'s','MarkerEdgeColor',col{i});
%         col{i}
            
            hold on
%         plot(testX(subA,3),f2s(subA)+offsetA,'Color',[.5, .5, .5]);
% 		plot(testX(subA,3),f2s(subA)+sd2s(subA)+offsetA,'--','Color',[.5, .5, .5]);
% 		plot(testX(subA,3),f2s(subA)-sd2s(subA)+offsetA,'--','Color',[.5, .5, .5]);
% 		plot(testX(subA,3),f2s(subA)+2*sd2s(subA)+offsetA,':','Color',[.5, .5, .5]);
% 		plot(testX(subA,3),f2s(subA)-2*sd2s(subA)+offsetA,':','Color',[.5, .5, .5]);
% 
% 		if j==3
%                     %			title('Regional + Local');
% 		elseif j==4
%                     %	title('Regional + Local non-linear');
% 		elseif j==6
%                     %	title('Greenland');
% 		end
% 
		ylabel('m','Color','k');
		if length(wxlim)==0
    		wxlim = xlim0;
%     		wxlim(1) = max(floor(min(testX(subA,3))/100)*100-1000,xlim0(1));
            %   if length(subB)>0
            %       wxlim(1) = min(wxlim(1),floor(min(time1(subB))/100)*100);
            %   end
        end
		xlim(wxlim);
% %        ylim([-35e3 20e3]);
%         set(gca,'YTick',[-30000:10000:20000]);
%         %set(gca,'YTick',[-40000:10000:10000]);
%         set(gca,'YTickLabel',{'-30','-20','-10','0','10','20'});
%         %set(gca,'YTickLabel',{'-35','-30','-25','-20','-15','-10','-5','0','5','10',})
    
       %legend([ps{1} ps{2} ps{3} ps{4} ],{'SE Florida';'Biscayne N.P.';'Florida Keys';'Dry Tortugas'});
       %legend([ps{1} ps{2} ps{3} ps{4} ps{5}],{'SE Florida';'Biscayne N.P.';'Florida Keys';'Dry Tortugas';'Everglades'});
	
%     if difftimestep>0
% 
%         Mdiff = bsxfun(@eq,testX(:,3),testX(:,3)')-bsxfun(@eq,testX(:,3),testX(:,3)'+ difftimestep);
%         Mdiff = Mdiff .* bsxfun(@eq,testreg,testreg');
%         sub=find(sum(Mdiff,2)==0);
%         Mdiff=Mdiff(sub,:);
%         difftimes=bsxfun(@rdivide,abs(Mdiff)*testX(:,3),sum(abs(Mdiff),2));;
%         diffreg=bsxfun(@rdivide,abs(Mdiff)*testreg,sum(abs(Mdiff),2));;
%         Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testX(:,3));
% 
%         clear df2s dV2s dsd2s;
%         for n=1:size(f2s,2)
%             df2s(:,n)=Mdiff*f2s(:,n);
%             dV2s(:,:,n)=Mdiff*V2s(:,:,n)*Mdiff';
% %             dsd2s(:,n)=sqrt(diag(dV2s(:,:,n)));
%         end
% 
%         subA2 = find(diffreg == testsites(1,1));
% 
%         for k=1:length(js)
%             offsetA=0;
%             hp(k+length(js))=subplot(numrows,length(js),k+length(js));
%             box on;
% 
%             hold on
% %             j=js(1,k);
% %             plot(difftimes(subA2),df2s(subA2),'Color',[.5, .5, .5]);
% %             plot(difftimes(subA2),df2s(subA2)+dsd2s(subA2),'--','Color',[.5, .5, .5]);
% %             plot(difftimes(subA2),df2s(subA2)-dsd2s(subA2),'--','Color',[.5, .5, .5]);
% %             plot(difftimes(subA2),df2s(subA2)+2*dsd2s(subA2),':','Color',[.5, .5, .5]);
% %             plot(difftimes(subA2),df2s(subA2)-2*dsd2s(subA2),':','Color',[.5, .5, .5]);
% 
%             ylabel(['- mm/y (= - m/ka)']);
%             xlim(wxlim);
% %            ylim([-20 10]);
%         end
%     end


    title(labl)
    %longticks(hp);
%    try; [bh,th]=label(hp,'ul',12,[],0,1,1,1.5,1.5); end

%     pdfwrite(['sl_' labl])    
% sl_rate{i}=[difftimes(subA2),df2s(subA2),dsd2s(subA2)];
end

