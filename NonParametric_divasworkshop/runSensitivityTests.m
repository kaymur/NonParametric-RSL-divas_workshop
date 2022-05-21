%% runSensitivityTests: sensitivity tests using synthetic data

datGIA_temp = importdata(fullfile(IFILES,'GIAs.csv'));
datGIA = datGIA_temp.data;

if run_start==1
    fid=fopen(['Sensitivity_Results_' df '.tsv'],'a');
    fprintf(fid,['Truth\tDate\tSeed\tRun\tN\tCI 95 Width\tMeanErr\tMeanAbsErr\tMeanRateErr\tMax Med Rate\t67 CI Width Rate\tCI 67 Width\tMed Err\tCoverage 95\tCoverage 67\tCov Rate 95\tCov Rate 67\tAbsol Rate Err\tRate Width 95\tMax 975 Rate\tMax 025 Rate\tMax 833 Rate\tMax 167 Rate\tMSE']);
    fprintf(fid,'\n');    
    fclose(fid); 
    fid2=fopen(['Sensitivity_Output.tsv'],'a');
    fprintf(fid2,['Truth\tDate\tSeed\tRun\tN\tupper_95_f']);
    fprintf(fid2,'\n');    
    fclose(fid2);
end

counter =run_start-1;
all=[];

fprintf('Truth Flag %0.0f\t Run Start %0.0f\n', [truth_flag run_start]);

for jjj = run_start:run_end
    clearvars -except distKern0 distKernel distKernFlor jjj all run_start run_end counter CEFILES IFILES synthflag truth_flag Seed dateField rng datGIA sfile distNorm distKern df pd
    runSetUp;
    %    load(sfile);
    modno=3;
    if mod(jjj ,3)==1
        n=1; 
    elseif mod(jjj ,3)==2
        n=5;
    else
        n=10;
    end
    N=n*12;
    
    if mod(jjj -1,6)<3
        t_unc=75;
    else
        t_unc=250;
    end

    if jjj < 7
        dat_type=1; %a_palmata
    elseif jjj < 13
        dat_type=2; %orbicella
    elseif jjj < 19
        dat_type=3; %sedimentary 
    elseif jjj < 25
        dat_type=4; %a_palmata and sedimentary 
    elseif jjj < 31
        dat_type=5; %limiting
    elseif jjj < 37
        dat_type=6; %limiting and orbicella 
    elseif jjj < 43
        dat_type=7; %limiting and sedimentary
    elseif jjj < 49
        dat_type=8; 
    elseif jjj < 55
        dat_type=9; %limiting and sedimentary
    elseif jjj < 61
        dat_type=10; 
    end

    counter = counter + 1;
        
    all=[all; modno n t_unc dat_type];

    SyntheticSensitivityTesting;
    
    truth =truths';
    ave_95_wi=mean(f95u-f95l);
    ave_67_wi=mean(f67u-f67l);
    ave_mean_err=mean(meanf-truth);
    mean_abs_err=mean(abs(medf-truth));
    ave_med_err=mean(medf-truth);
    ave_mean_err_rate=mean(medr'-truerate);
    ave_abs_mean_err_rate=mean(abs(medr'-truerate));
    max_thets=max(thets(:,[1 2 5]));
    min_thets=min(thets(:,[1 2 5]));
    upp95_thets=quantile(thets(:,[1 2 5]),.975);
    low95_thets=quantile(thets(:,[1 2 5]),.025);
    med_thets=median(thets(:,[1 2 5]));
    rate_wi_95=mean(r95u-r95l);
    rate_wi_67=mean(r67u-r67l);

    fid=fopen(['../Sensitivity_Results_' df '.tsv'],'a');
    fprintf(fid,'%0.f\t', truth_flag);
    fprintf(fid,'%0s\t', df);
    fprintf(fid,'%0.f\t', Seed);
    fprintf(fid,'%0.f\t', jjj );
    fprintf(fid,'%0.3f\t', N);
    fprintf(fid,'%0.3f\t', ave_95_wi);
    fprintf(fid,'%0.3f\t', ave_mean_err);
    fprintf(fid,'%0.3f\t', mean_abs_err);
    fprintf(fid,'%0.3f\t', ave_mean_err_rate);
    fprintf(fid,'%0.3f\t', max(medr));
    fprintf(fid,'%0.3f\t', rate_wi_67);
    fprintf(fid,'%0.3f\t', ave_67_wi);    
    fprintf(fid,'%0.3f\t', ave_med_err);
    fprintf(fid,'%0.3f\t', cov_95);
    fprintf(fid,'%0.3f\t', cov_67);
    fprintf(fid,'%0.3f\t', covr_95);
    fprintf(fid,'%0.3f\t', covr_67);
    fprintf(fid,'%0.3f\t', ave_abs_mean_err_rate);
    fprintf(fid,'%0.3f\t', rate_wi_95);
    fprintf(fid,'%0.3f\t', max(r95u));
    fprintf(fid,'%0.3f\t', max(r95l));
    fprintf(fid,'%0.3f\t', max(r67u));
    fprintf(fid,'%0.3f\t', max(r67l));
    fprintf(fid,'%0.3f\t', srmse);
    %srmse = sqrt(sum(resids.^2)/length(resids))

    fprintf(fid,'\n');    
    fclose(fid);
    
fid2=fopen(['Sensitivity_Output.tsv'],'a');
    fprintf(fid2,'%0.1f\t', truth_flag);
    fprintf(fid2,'%0s\t', dateField);
    fprintf(fid2,'%0.1f\t', Seed);
    fprintf(fid2,'%0.3f\t', jjj );
    fprintf(fid2,'%0.3f\t', N);
    for ii=1:length(maxf)
        fprintf(fid2,'%0.3f\t', maxf(ii));
    end
    fprintf(fid2,'\t');
    for ii=1:length(medf)
        fprintf(fid2,'%0.3f\t', medf(ii));
    end
    fprintf(fid2,'\t');
    for ii=1:length(minf)
        fprintf(fid2,'%0.3f\t', minf(ii));
    end
    fprintf(fid2,'\t');
    for ii=1:length(max_thets)
        fprintf(fid2,'%0.3f\t', max_thets(ii));
    end
    fprintf(fid2,'\t');
    for ii=1:length(med_thets)
        fprintf(fid2,'%0.3f\t', med_thets(ii));
    end
    fprintf(fid2,'\t');
    for ii=1:length(min_thets)
        fprintf(fid2,'%0.3f\t', min_thets(ii));
    end
    fprintf(fid2,'\t');
    for ii=1:length(medr)
        fprintf(fid2,'%0.3f\t', medr(ii));
    end   
    fprintf(fid2,'\t');
    for ii=1:length(r95u)
        fprintf(fid2,'%0.3f\t', r95u(ii));
    end
    fprintf(fid2,'\t');
    for ii=1:length(r95l)
        fprintf(fid2,'%0.3f\t', r95l(ii));
    end

    fprintf(fid2,'\n');
    fclose(fid2);
    save('matlab.mat','-v7.3');
end
