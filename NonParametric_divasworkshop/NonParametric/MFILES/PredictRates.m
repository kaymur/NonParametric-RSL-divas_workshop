%% find the rates and summary statistics for validation and sensitivity test results

max_thets=max(thets');
min_thets=min(thets');
med_thets=median(thets');

fid=fopen(['PredictionsAndThetas.tsv'],'a');
fprintf(fid,'%0.2f\t',medf(end));
fprintf(fid,'%0.2f\t',medf(end-10));
fprintf(fid,'%0.2f\t',medf(end-20));
fprintf(fid,'%0.2f\t',medf(end-30));
fprintf(fid,'%0.2f\t',medf(end-40));
fprintf(fid,'%0.2f\t',medf(end-50));
if size(medf,2)>60
    fprintf(fid,'%0.2f\t',medf(end-60));
    fprintf(fid,'%0.2f\t',medf(end-70));
    fprintf(fid,'%0.2f\t',medf(end-80));
    fprintf(fid,'%0.2f\t',medf(end-90));
    fprintf(fid,'%0.2f\t',medf(end-100));
    fprintf(fid,'%0.2f\t',medf(end-110));
    fprintf(fid,'%0.2f\t',medf(end-120));
end
fprintf(fid,'%0.2f\t',f95u(end)-medf(end));
fprintf(fid,'%0.2f\t',f95u(end-10)-medf(end-10));
fprintf(fid,'%0.2f\t',f95u(end-20)-medf(end-20));
fprintf(fid,'%0.2f\t',f95u(end-30)-medf(end-30));
fprintf(fid,'%0.2f\t',f95u(end-40)-medf(end-40));
fprintf(fid,'%0.2f\t',f95u(end-50)-medf(end-50));
if size(medf,2)>60
    fprintf(fid,'%0.2f\t',f95u(end-60)-medf(end-60));
    fprintf(fid,'%0.2f\t',f95u(end-70)-medf(end-70));
    fprintf(fid,'%0.2f\t',f95u(end-80)-medf(end-80));
    fprintf(fid,'%0.2f\t',f95u(end-90)-medf(end-90));
    fprintf(fid,'%0.2f\t',f95u(end-100)-medf(end-100));
    fprintf(fid,'%0.2f\t',f95u(end-110)-medf(end-110));
    fprintf(fid,'%0.2f\t',f95u(end-120)-medf(end-120));
end
fprintf(fid,'%0.2f\t',medf(end)-f95l(end));
fprintf(fid,'%0.2f\t',medf(end-10)-f95l(end-10));
fprintf(fid,'%0.2f\t',medf(end-20)-f95l(end-20));
fprintf(fid,'%0.2f\t',medf(end-30)-f95l(end-30));
fprintf(fid,'%0.2f\t',medf(end-40)-f95l(end-40));
fprintf(fid,'%0.2f\t',medf(end-50)-f95l(end-50));
if size(medf,2)>60

fprintf(fid,'%0.2f\t',medf(end-60)-f95l(end-60));
fprintf(fid,'%0.2f\t',medf(end-70)-f95l(end-70));
fprintf(fid,'%0.2f\t',medf(end-80)-f95l(end-80));
fprintf(fid,'%0.2f\t',medf(end-90)-f95l(end-90));
fprintf(fid,'%0.2f\t',medf(end-100)-f95l(end-100));
fprintf(fid,'%0.2f\t',medf(end-110)-f95l(end-110));
fprintf(fid,'%0.2f\t',medf(end-120)-f95l(end-120));
end
fprintf(fid,'%0.5f\t',med_thets(1));
fprintf(fid,'%0.5f\t',med_thets(2));
fprintf(fid,'%0.5f\t',med_thets(end));
fprintf(fid,'%0.5f\t',max_thets(1));
fprintf(fid,'%0.5f\t',max_thets(2));
fprintf(fid,'%0.5f\t',max_thets(5));
fprintf(fid,'%0.5f\t',min_thets(1));
fprintf(fid,'%0.5f\t',min_thets(2));
fprintf(fid,'%0.5f\n',min_thets(end));

fclose(fid);

ave_95_wi=mean(f95u-f95l);
ave_90_wi=mean(maxf-minf);

%% for sensitivity tests, look at rates over 500 year periods
if synthflag==1
    difftimestep=500;
elseif barbflag==1
    difftimestep = [200 400 1000];
else
    difftimestep=1000;
end

%% Plot rates
testX = [mean(lat)*ones(length(testt),1) mean(long)*ones(length(testt),1) testt'];
for kk=1:length(difftimestep)   % = [200 400 1000];
    Mdiff = bsxfun(@eq,testX(:,3),testX(:,3)')-bsxfun(@eq,testX(:,3),testX(:,3)'+ difftimestep(kk));
    testreg=datid(1)*ones(length(testt),1);
    Mdiff = Mdiff .* bsxfun(@eq,testreg,testreg');
    sub=find(sum(Mdiff,2)==0);
    Mdiff=Mdiff(sub,:);
    difftimes=bsxfun(@rdivide,abs(Mdiff)*testX(:,3),sum(abs(Mdiff),2));
    diffreg=bsxfun(@rdivide,abs(Mdiff)*testreg,sum(abs(Mdiff),2));
    Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testX(:,3));
    clear dfs df2s dV2s dsd2s;

    fp=fs_all';
    for n=1:size(fp,2)    
        dfs(:,n)=Mdiff*fp(:,n);  
    end
    max_rate=max(dfs');
    r95u=quantile(dfs',.975);
    r95l=quantile(dfs',.025);
    r67u=quantile(dfs',.833);
    r67l=quantile(dfs',.167);
    medr=median(dfs');
    meanr=mean(dfs');
    maxrates=max(dfs');
    mu67r=quantile(dfs',.833);
    ml67r=quantile(dfs',.133);
    ru975=quantile(dfs',.975);
    rl025=quantile(dfs',.025);
    medrates=median(dfs');
    
    clf; 
    plotdat.x=1950-testt';
    plotdat.y=medf';
    plotdat.dy=[f95u'-medf' medf'-f95l'] ;
    PlotWithShadedErrors(plotdat,[0 0 0],'none');

    hold on;

    MakePlotsSampling(PXholo,fs(:,1),sds(:,1),V2s,[],'Posterior',[],1000,[],[],meanSL,testt',coral_params);
    plot(1950-testt,truths,'r','LineWidth',1);

    if truth_flag==2
        ylim([-180e3 15e3]);
    else
        ylim([-70e3 10e3]);
    end
    clear title;
    title(['Posterior RSL ' num2str(jjj ) ' Seed ' num2str(Seed) ' SL ' num2str(truth_flag)]);
    pdfwrite(['Synth_post_RSL_' df '_' num2str(jjj ) '_Seed_' num2str(Seed) '_SL_' num2str(truth_flag)]);        

    clf;
    plotr.x=1950-difftimes;
    plotr.y=medr';
    plotr.dy=[medr'-r95l' r95u'-medr'] ;
    PlotWithShadedErrors(plotr,[0 0 0],[.9 .9 .9]);
    hold on;
    plotr.dy=[medr'-r67l' r67u'-medr'] ;
    PlotWithShadedErrors(plotr,[0 0 0],[.9 .9 .9],'-','-.');
    plot(1950-difftimes,truerate,'r','LineWidth',1);
    if truth_flag==2
        ylim([-80 100]);
    else
        ylim([-100 100]);
    end
    if truth_flag==2
        ylim([-20 50]);
    else
        ylim([-30 30]);
    end
    title(['Rates Run ' num2str(jjj ) ' Seed ' num2str(Seed) ' SL ' num2str(truth_flag)]);
    pdfwrite(['Synth_post_Rates_' df '_' num2str(jjj) '_Seed_' num2str(Seed) '_SL_' num2str(truth_flag)]);        

    %% true rate for synthetic tests
    hold on;
    if difftimestep(kk)==1000
        ylim([-20 30]);
        pdfwrite(['Post_rates_Mod_' num2str(modno) '_1000-yr-av']);
    elseif difftimestep(kk)==100
        ylim([-40 50]);
        pdfwrite(['Post_rates_Mod_' num2str(modno) '_100-yr-av']);
    elseif difftimestep(kk)==200
        ylim([-40 50]);
        pdfwrite(['Post_rates_Mod_' num2str(modno) '_200-yr-av']);
    elseif difftimestep(kk)==400
        ylim([-40 50]);
        pdfwrite(['Post_rates_Mod_' num2str(modno) '_400-yr-av']);
    elseif difftimestep(kk)==500
        ylim([-30 40]);
        pdfwrite(['Post_rates_Mod_' num2str(modno) '_500-yr-av']);
    end
     if exist('Rates.tsv')==0
        fid=fopen(['Rates.tsv'],'a');
        fprintf(fid,['Rates (mm/y)\n']);
        fprintf(fid,['Model\t']);
        fprintf(fid,['Run No\t']);
        for pp=1:(length(testX)-1)/10
            fprintf(fid,'Rate (avg, %0.0f-%0.0f)\t',[1950-testX((pp-1)*10+1,3) 1950-testX(pp*10+1,3)]);
        end
        for pp=1:(length(testX)-1)/10
            fprintf(fid,'UErr (avg, %0.0f-%0.0f)\t',[1950-testX((pp-1)*10+1,3) 1950-testX(pp*10+1,3)]);
        end
        for pp=1:(length(testX)-1)/10
            fprintf(fid,'LErr (avg, %0.0f-%0.0f)\t',[1950-testX((pp-1)*10+1,3) 1950-testX(pp*10+1,3)]);
        end
        for pp=1:(length(testX)-1)/10
            fprintf(fid,'Upper67percentile (avg, %0.0f-%0.0f)\t',[1950-testX((pp-1)*10+1,3) 1950-testX(pp*10+1,3)]);
        end
        for pp=1:(length(testX)-1)/10
            fprintf(fid,'Lower67percentile (avg, %0.0f-%0.0f)\t',[1950-testX((pp-1)*10+1,3) 1950-testX(pp*10+1,3)]);
        end
    else
        fid=fopen(['Rates.tsv'],'a');
    end
    fprintf(fid,'\n');
    fprintf(fid,'%0.0f\t',modno);
    for x=1:(length(medrates)-1)/10+1
        fprintf(fid,'%0.2f\t',medrates((x-1)*10+1));
    end
    for x=1:(length(medrates)-1)/10+1
        fprintf(fid,'%0.2f\t',ru975((x-1)*10+1)-medrates((x-1)*10+1));
    end
    for x=1:(length(medrates)-1)/10+1
        fprintf(fid,'%0.2f\t',medrates((x-1)*10+1)-rl025((x-1)*10+1));
    end
    for x=1:(length(medrates)-1)/10+1
        fprintf(fid,'%0.2f\t',mu67r((x-1)*10+1)-medrates((x-1)*10+1));
    end
    for x=1:(length(medrates)-1)/10+1
        fprintf(fid,'%0.2f\t',medrates((x-1)*10+1)-ml67r((x-1)*10+1));
    end
    fclose(fid);
end

covr_95=(length(truerate')-(sum(truerate'>r95u)+sum(truerate'<r95l)))/length(truerate');
covr_67=(length(truerate')-(sum(truerate'>r67u)+sum(truerate'<r67l)))/length(truerate');
