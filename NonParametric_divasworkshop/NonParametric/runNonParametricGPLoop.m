%% For Pacific database
%% Non-parametric methods
% Loops through all data files in a given folder
% And run the model for each of these data files
% Make sure you are in the 'NonParametric' working directory
clear
clc

%% USER INPUTS HERE 

% input name of folder containing files to loop through
% folder should be contained within the 'IFILES' folder
folder = '/Pacific RSL data'; 

%% DO NOT EDIT THIS SECTION

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

            clearvars -except pd df IFILES filelist IFILES distFile run_type;
            
            if exist('Seed', 'Var')
            elseif exist('run_type','Var')
            else
                Seed = 3; 
                truth_flag = 2;
                run_type = 1;
                run_start = 36;
                run_end = 36;
            end
            
            if ~exist('pd','var')
                pd = pwd;
            end
            
            NCT = maxNumCompThreads;
            fprintf('%0.0f computational threads\n',[maxNumCompThreads]);        
            LASTN = maxNumCompThreads(NCT);
            fprintf('%0.0f same number computational threads\n',[LASTN]);        
            newNCT = maxNumCompThreads(2);
            fprintf('%0.0f new comp threads\n',[newNCT]);  
            
            %% define which test we are doing here:
            % 1 for time-series sensitivity tests
            % 2 for time-series synthetic validation 
            
            %% run set-up scripts (define the directories and import data)
            addpath([pd '/MFILES']);
            addpath(pd);
            runSetUp_loop;
                
            if run_type == 1           % sensitivity tests
                synthflag = 1;
                runSensitivityTests;
            elseif run_type == 2       % validation runs
                synthflag = 2;
                runValidationTests;
            else
                prepData;
                sample_Ys_Thetas;
                finishAnalysis;
                save %"finalsamps" nthin thetachange nburn stepchange thinned_ys steps_keep logp_keep Nsamples trainsubz all_ys cspecies
            end

end



