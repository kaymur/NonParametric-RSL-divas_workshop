all_ys = repmat(Y,1,nsamps);
trainsubz = find(limiting~=0&limiting~=20);
if spatial==1 && (modno==6||modno==3)
    DefCov_spatial;
else
    DefineCovarianceFunctions;
end
PXholo.datid=round(datid);
PXholo.region=region;
PXholo.time1=time1;
PXholo.time2=time2;
PXholo.limiting=limiting;
PXholo.indic=indic;
PXholo.Y=Y;
PXholo.Y0=Y0;
PXholo.dY = dY;
PXholo.dY0=dY0;
PXholo.cspecies=cspecies;
PXholo.lat=lat;
PXholo.long=long;
PXholo.Ycv=sparse(diag(dY.^2));
PXholo.siteid=round(siteid);
PXholo.sitenames=sitenames;
%PXholo.meantime=(PXholo.time1+PXholo.time2)/2;
PXholo.meantime=meantime;
PXholo.sitecoords = sitecoords;
PXholo.sitelen = sitelen;
PXholo.dY_init=dY_init;
PXholo.indic=indic;
PXholo.dt=dt;

clear datasets;
dataset=PXholo;
datasets{1}=PXholo;
datasets{1}.label='PXholo';


t1=datasets{1}.time1; t2=datasets{1}.time2;
datasets{1}.long = mod(datasets{1}.long,360); 
sub=find(datasets{1}.long>180); datasets{1}.long(sub)=datasets{1}.long(sub)-360;
datasets{1}.meantime=mean([t1 t2],2);
datasets{1}.dt = abs(t1-t2)/4;

step_size=dY*2.5; 
