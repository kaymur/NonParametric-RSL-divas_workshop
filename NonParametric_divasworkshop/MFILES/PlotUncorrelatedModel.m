%% make uncorrelated plots
findQuantiles;
flo=find(region==20);
dataset1.datid=datid(flo);
dataset1.Y=Y(flo);
dataset1.limiting=limiting(flo);
dataset1.lat=lat(flo);
dataset1.long=long(flo);
dataset1.cspecies=cspecies(flo);
dataset1.indic=indic(flo);
dataset1.dt=dt(flo);
dataset1.meantime=meantime(flo);

dataset1.sitecoords=[25.03 -80.68];
clf; 
if exist('testlocs')
    testloc =testlocs;
elseif exist('testlocp')
    testloc=testlocp;
end
makeplots_percentile_colors;
