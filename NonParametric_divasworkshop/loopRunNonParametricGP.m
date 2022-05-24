% Loop through all data files in a given folder
% And run the runNonParametricGP.m script for each of these data files
% make sure you are in the 'NonParametric' working directory
clear
clc

%% User inputs here 

% input name of folder containing files to loop through
% folder should be contained within the 'IFILES' folder
folder = '/Pacific RSL data'; 

%% Loop to run the runNonParametricGP.m script for all files in a folder 

% get list of files within this folder 
pd=pwd; % to make sure each time it goes back to NonParametric main folder before running
IFILES=[pd '/IFILES' folder];
filelist = {dir(IFILES).name};
filelist={filelist{3:length(filelist)}};

for i = 1:length(filelist)
    
    % user inputs required for runNonParametricGP.m:
        % distFile
        % df
        % run type
    
    run_type = 3;  % We are using real data, not running tests on the model.
    distFile=filelist{i}; 
    name = distFile(1:length(distFile)-4); 
    df = [datestr(date,'ddmmyyyy') ' ' name]; 

    % run the runNonParametricGP.m for file i within filelist
    runNonParametricGP_loop


end

