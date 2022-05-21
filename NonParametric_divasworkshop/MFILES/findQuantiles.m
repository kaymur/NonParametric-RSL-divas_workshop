index=[0:.010:120];    quants={};
    figure;
for iii=1:length(distKernel)
    probs=cdf(distKernel{iii},index');
    quants995=min(index(find(probs>=.995)));
    %quants995=mean(index(find(probs>=.995)));
    %quants995=mean(index(intersect(find(probs>=.990),find(probs<=.998))));
    %quants995=mean(index(intersect(find(probs>=.990),find(probs<=.998))));
    quants990=min(index(find(probs>=.99)));
    quants975=min(index(find(probs>=.975)));
    quants950=min(index(find(probs>=.950)));
    quants900=min(index(find(probs>=.900)));
    quants850=min(index(find(probs>=.850)));
    quants833=min(index(find(probs>=.8333)));
    quants750=min(index(find(probs>=.750)));
    quants500=min(index(find(probs>=.500)));
    quants250=min(index(find(probs>=.250)));
    quants167=min(index(find(probs>=.1667)));
    quants150=min(index(find(probs>=.15)));
    quants100=min(index(find(probs>=.10)));
    quants050=min(index(find(probs>=.05)));
    %quants025=mean(index(intersect(find(probs<=.026),find(probs>=.024))));
    quants025=min(index(find(probs>=.025)));
    quants010=min(index(find(probs>=.010)));
    quants005=min(index(find(probs>=.005)));
    %mean(index(intersect(find(probs<=.01)
%     quants99=[quantile(probs,.005) quantile(probs,.995)];
%     quants98=[quantile(probs,.01) quantile(probs,.99)];
%     quants95=[quantile(probs,.025) quantile(probs,.975)];
%     quants90=[quantile(probs,.05) quantile(probs,.95)];
%     quants80=[quantile(probs,.10) quantile(probs,.9)];
%     quants70=[quantile(probs,.15) quantile(probs,.85)];
%     quants67=[quantile(probs,.133) quantile(probs,.867)];
%     quants50=[quantile(probs,.25) quantile(probs,.75)];
    quants{iii}=[quants005 quants995  quants025 quants975  quants050 quants950 quants167 quants833 quants250 quants750].*1000;
    %     quants{iii}=[quants010 quants990 quants050 quants950 quants100 quants900 quants133 quants867 quants250 quants750];
    %quants{iii}=[quants99 quants98 quants95 quants90 quants80];
%     quants{iii}=[quants99 quants95 quants90 quants80 quants67];
    hold on;
    plot([iii iii],quants{iii}(1:2),'Color',[0.75 .88 .98],'LineWidth',5);
    plot([iii iii],quants{iii}(3:4),'Color',[.9 1 0],'LineWidth',5); 
    plot([iii iii],quants{iii}(5:6),'Color',[.9 .45 0.15],'LineWidth',5);
    plot([iii iii],quants{iii}(7:8),'Color',[.8 0.05 0.05],'LineWidth',5);
    plot([iii iii],quants{iii}(9:10),'Color',[.5 0.5 0.5],'LineWidth',5);
end
iii=iii+1;
index=[-10000:10:10000];    
    probs=normcdf(index',0,550);
    %probs=normcdf(index',0,550);
    quants995=min(index(find(probs>=.995)));
    quants990=min(index(find(probs>=.99)));
    quants975=min(index(find(probs>=.975)));
    quants950=min(index(find(probs>=.950)));
    quants900=min(index(find(probs>=.900)));
    quants850=min(index(find(probs>=.850)));
    quants833=min(index(find(probs>=.8333)));
    quants750=min(index(find(probs>=.750)));
    quants500=min(index(find(probs>=.500)));
    quants250=min(index(find(probs>=.250)));
    quants167=min(index(find(probs>=.1667)));
    quants150=min(index(find(probs>=.15)));
    quants100=min(index(find(probs>=.10)));
    quants050=min(index(find(probs>=.05)));
    quants025=min(index(find(probs>=.025)));
    quants010=min(index(find(probs>=.010)));
    quants005=min(index(find(probs>=.005)));
    quants{iii}=[quants005 quants995  quants025 quants975  quants050 quants950 quants167 quants833 quants250 quants750] ;
    hold on;
    plot([iii iii],quants{iii}(1:2),'Color',[0.75 .88 .98],'LineWidth',5);
    plot([iii iii],quants{iii}(3:4),'Color',[.9 1 0],'LineWidth',5); 
    plot([iii iii],quants{iii}(5:6),'Color',[.9 .45 0.15],'LineWidth',5);
    plot([iii iii],quants{iii}(7:8),'Color',[.8 0.05 0.05],'LineWidth',5);
    plot([iii iii],quants{iii}(9:10),'Color',[.5 0.5 0.5],'LineWidth',5);


    legend('99% CI range','95% CI range','90% CI range','66% CI range','50% CI range');

    pdfwrite('quantiles all data');
