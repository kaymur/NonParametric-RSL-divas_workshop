colors=[ 1 .55 0.05
        1 0 0
        .68 .20 .70
        0 0 0
        0 1 0
        1 .22 .70
        0.2 0.27 .9
        1 .85 .1
        0 .96 .86
        0 .45 .1
        .18 .09 .34
        ];

for jj=1:size(sitecoords,1) 
    if spatial == 0
        clf;
        hold on;
        subB=[1:N];
    else
        subB=find(datid==siteid(jj));
    end
    means = coral_params(:,3)';
    for nn = 1:length(subB)
        if limiting(subB(nn))== 0 && cspecies(subB(nn))==0
            rectangle('Position',[ages(subB(nn))-2*dt(subB(nn)),Y_final(subB(nn))-2*dY(subB(nn)),4*dt(subB(nn)),4*dY(subB(nn))],'Edgecolor',[0,0,0]);  % peat index points basal, black
        elseif limiting(subB(nn))==0 && cspecies(subB(nn))>=10
            rectangle('Position',[ages(subB(nn))-2*dt(subB(nn)),Y_final(subB(nn))-2*dY(subB(nn)),4*dt(subB(nn)),4*dY(subB(nn))],'Edgecolor',[0,0,0]);  % peat index points basal, black
        elseif cspecies(subB(nn))==9 %&&limiting(subB(nn))==0
                plot([1 1]*ages(subB(nn)),[Y_final(subB(nn))-2*dY(subB(nn)) Y_final(subB(nn))+2*means(cspecies(subB(nn)))*1000+2*dY(subB(nn))],'Color',colors(cspecies(subB(nn)),:));
                plot([ages(subB(nn))-2*dt(subB(nn)) ages(subB(nn))+2*dt(subB(nn))],[1 1]*Y_final(subB(nn)),'Color',colors(cspecies(subB(nn)),:));
        elseif cspecies(subB(nn))==1 && limiting(subB(nn))==0
            rectangle('Position',[ages(subB(nn))-2*dt(subB(nn)),Y_final(subB(nn))+coral_params(cspecies(subB(nn)),1)*1000-2000*coral_params(cspecies(subB(nn)),2),...
                4*dt(subB(nn)),4000*coral_params(cspecies(subB(nn)),2)],'Edgecolor',colors(cspecies(subB(nn)),:));  % peat index points basal, turquoise
        elseif limiting(subB(nn))==2 % all nonparametric
            plot([1 1]*ages(subB(nn)),[Y_final(subB(nn))-2*dY(subB(nn)) Y_final(subB(nn)) + means(cspecies(subB(nn)))*1000+2*dY(subB(nn))],'Color',colors(cspecies(subB(nn)),:));
            plot([ages(subB(nn))-2*dt(subB(nn)) ages(subB(nn))+2*dt(subB(nn))],[1 1]*Y_final(subB(nn)),'Color',colors(cspecies(subB(nn)),:));
        elseif limiting(subB(nn))==3 % uniform
            rectangle('Position',[ages(subB(nn))+2*dt(subB(nn)),Y_final(subB(nn))-2*dY(subB(nn)),4*dt(subB(nn)),4*dY(subB(nn))],'Edgecolor',[.8,.8,.8]);  % uniform, grey
            hold on
        end
    end
    for nn = 1:length(subB)
        if limiting(subB(nn)) == -1 % marine limiting
            plot([1 1]*ages(subB(nn)),Y_final(subB(nn))+[-2 2]*dY(subB(nn)),'Color',[.1,.35,.7]);
            plot([ages(subB(nn))-2*dt(subB(nn)) ages(subB(nn))+2*dt(subB(nn))],[1 1]*(Y_final(subB(nn))+2*dY(subB(nn))),'Color',[.1,.35,.7]);
        elseif limiting(subB(nn))==1 % terrestrial limiting
            plot([1 1]*ages(subB(nn)),Y_final(subB(nn))+[-2 2]*dY(subB(nn)),'Color',[.05,.7,.5]);
            plot([ages(subB(nn))-2*dt(subB(nn)) ages(subB(nn))+2*dt(subB(nn))],[1 1]*(Y_final(subB(nn))-2*dY(subB(nn))),'Color',[.05,.7,.5]); 
        end
    end
    [tT, xi] = sort(age_true(subB));
    tY = Y_true(subB(xi));
    plot(tT,tY,'r','LineWidth',1);
    ylim([-45000 25000]);
    pdfwrite(['Data_v_Truth_' sitenames{jj}]);
end
