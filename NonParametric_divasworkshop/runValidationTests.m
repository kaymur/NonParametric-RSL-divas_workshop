%% compare results of full non-parametric, Gaussian, and uncorrelated models, 
% all with same synthetic datasets
distKernel=distKern;

first_flag=0;
jjj=1;

% create synthetic data for validation experiments
datGIA_temp = importdata(fullfile(IFILES,'GIAs.csv'));
datGIA = datGIA_temp.data;

%****
CreateValidationData;

    synthflag = 2;
    cor=find(cspecies>0);
    noncor=find(cspecies==0);
    indic(noncor)=1;
    
    meantime = ages;
    time1 = meantime-2*dt;
    time2 = meantime+2*dt;
    Y = Y_fin;
    testt0 = 1950-[0:100:12000];
    
%% check to make sure data looks correct
%PlotValidationData;
    
SampleCompareThreeMods;
trainsubx=train_all; collinear=[];
wdataset=datasets{1};
if isempty(y_accepted)
    ys0=Y(trainsub);
    ys=[];
    thetas=thets;
else
    if size(thets,2)>1
        ss=min(size(thets,2),size(y_accepted,2));
        thetas=thets(:,1:ss);
        ys=y_accepted(:,1:ss);
    else
        ys=y_accepted;
        thetas=thets;
    end
end
if ~exist('meanSL','var')
    meanSL = 0;
end

ModelEvaluation;

% plot the data as in Hibbert et al. & evaluate as such
testloc=testlocp;
PlotUncorrelatedModel;

    fid=fopen(['Val_Stats.tsv'],'a');
    fprintf(fid,['Model\tDate\tSeed\tFinalLogLike\t95 CI\t67 CI\tMeanResid\tMedianResids\tAbsMeanResids\tMSE\tCoverage(95)\tCoverage(67)\n']);
    fprintf(fid,'%0.f\t', modno);
    fprintf(fid,'%0s\t', dateField);
    fprintf(fid,'%0.f\t', Seed);
    fprintf(fid,'%0.3f\t', fin_ll );
    fprintf(fid,'%0.3f\t', ave_95_wi);
    fprintf(fid,'%0.3f\t', ave_67_wi);
    fprintf(fid,'%0.3f\t', m_resid);
    fprintf(fid,'%0.3f\t', med_resid);
    fprintf(fid,'%0.3f\t', abs_mean_resid);
    fprintf(fid,'%0.3f\t', srmse);
    fprintf(fid,'%0.3f\t', cov_95);
    fprintf(fid,'%0.3f\t', cov_67);
    fprintf(fid,'\n');

    fprintf(fid,'Uncorrelated\t');
    fprintf(fid,'%0s\t', dateField);
    fprintf(fid,'%0.f\t', Seed);
    fprintf(fid,'%0.3f\t', logp_hibbert );
    fprintf(fid,'%0.3f\t', CI95);
    fprintf(fid,'%0.3f\t', CI67);
    fprintf(fid,'%0.3f\t', hib_m_resid);
    fprintf(fid,'%0.3f\t', hib_med_resid);
    fprintf(fid,'%0.3f\t', hib_abs_mean_resid);
    fprintf(fid,'%0.3f\t', hib_srmse);
    fprintf(fid,'%0.3f\t', hib_cov_95);
    fprintf(fid,'%0.3f\t', hib_cov_67);
    fprintf(fid,'\n');
    fclose(fid);

Nsamps=Nsamples;
save(['mat3_' num2str(Seed)]);

% reload the data, and run the Gaussian model
modno=6; nsamps = 25000;
trainsubx=train_all;
sub_unif=find(cspecies==1);
delev(sub_unif)=dY(sub_unif);
elev(sub_unif)=Y(sub_unif);
dY(sub_unif)= sqrt((dY(sub_unif)).^2+(2500)^2);  
Y(sub_unif)=Y(sub_unif)+2500;
limiting(sub_unif)=0;

sampnorm=0;

SampleCompareThreeMods;

trainsubx=train_all; collinear=[];
wdataset=datasets{1};
if isempty(y_accepted)
    ys=[];
    thetas=thets;
else
    if size(thets,2)>1
        ss=min(size(thets,2),size(y_accepted,2));
        thetas=thets(:,1:ss);
        ys=y_accepted(:,1:ss);
    else
        ys=y_accepted;
        thetas=thets;
    end
end
ModelEvaluation;
    fid=fopen(['Val_Stats.tsv'],'a');
    fprintf(fid,['Model\tDate\tSeed\tFinalLogLike\t95 CI\t67 CI\tMeanResid\tMedianResids\tAbsMeanResids\tMSE\n']);
    fprintf(fid,'%0.f\t', modno);
    fprintf(fid,'%0s\t', dateField);
    fprintf(fid,'%0.f\t', Seed);
    fprintf(fid,'%0.3f\t', fin_ll );
    fprintf(fid,'%0.3f\t', ave_95_wi);
    fprintf(fid,'%0.3f\t', ave_67_wi);
    fprintf(fid,'%0.3f\t', m_resid);
    fprintf(fid,'%0.3f\t', med_resid);
    fprintf(fid,'%0.3f\t', abs_mean_resid);
    fprintf(fid,'%0.3f\t', srmse);
    fprintf(fid,'%0.3f\t', cov_95);
    fprintf(fid,'%0.3f\t', cov_67);
    fprintf(fid,'\n');
    fclose(fid);
save(['mat6_' num2str(Seed)]);

