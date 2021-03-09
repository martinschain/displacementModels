function[s] = applyDisplacementModels()
% This function creates a structure s that is needed to run the 
% displacement models, and then executes all four models. The user 
% needs to edit the fields so that his/her data are processed, and 
% save the output from each model according to his/her preference. 
% Some parameter values are suggested, but most of these settings 
% can be modified without affecting the functionality of the code.
%_____________________________________________________________________
%                      Martin Schain, Neurobiology Research Unit, 2021 

% Generic parameters
s.nbrSteps = 3;  %The number of steps used in the multistep models
s.stepSize = 1/30; %The resolution of the time grid (1/30 min = 2 seconds)

% Parameters for the specific data, to be filled by user
s.scanDur = []; % scan duration (in minutes)
s.T = []; % time of intervention (in minutes).
 
s.TAC = []; % displacement PET time activity curve [px1] vector
s.tPET = []; % Midtimes of the PET frames [px1] vector
s.dur = []; % frame durations [px1] vector
s.weights = []; % Weights for the fit [px1] vector. We recommend sqrt(s.dur).

s.t =[]; % time grid for the model, could be 0:stepSize:s.scanDur. [nx1] vector
s.inFcn = []; % metabolite corrected arterial input funtion, [nx1] vector 
s.wb = []; %whole blood curve, [nx1] vector
% NB! Input function and whole blood should already be in model time,

% Optimization parameters
s.fitParams.teStart = (s.scanDur - s.T)/2; % time AFTER T (i.e., real te is s.T+te.
s.fitParams.options = optimset;
s.fitParams.options.MaxFunEvals = 10^4;
s.fitParams.options.MaxIter = 10^4;
s.fitParams.options.TolFun = .5*10^-3;
s.fitParams.option.TolX = .5*10^-3;

%% Execute 2TCM single step
s.fitParams.startParams = [.1 .1 .1 .1 .05 .5 s.T+2]; %K1 k2 k3 k4 vB occ tb
s.fitParams.lowBound = [0 0 0 0 0 0 s.T];
s.fitParams.upBound = [inf inf inf inf inf 1 s.scanDur-5]; %Hard upper bound on tb. 
[Cnd,Cs,model,occCurve,Ks,exitFlag] = singleStepModel_2TCM(s);
% save the output somehow

%% Execute 2TCM multistep
s.fitParams.startParams = [.1 .1 .1 .1 .05 s.fitParams.teStart .5]; %K1 k2 k3 k4 vB te occ
s.fitParams.lowBound = [0 0 0 0 0 0 0]; 
s.fitParams.upBound = [inf inf inf inf inf s.scanDur-5-s.T 1]; %Hard upper bound on te.
[Cnd,Cs,model,occCurve,discOcc,Ks,exitFlag] = multiStepModel_2TCM(s);
% save the output somehow

%% Execute 1TCM single step
s.fitParams.startParams = [.1 .1 .05 1 .5 s.T]; %K1 k2 vB BP occ T
s.fitParams.lowBound = [0 0 0 0 0 0]; 
s.fitParams.upBound = [inf inf inf inf 1 s.scanDur-5];
[model,occCurve,Ks,exitFlag] = singleStepModel_1TCM(s);
% save the output somehow

%% Execute 1TCM multistep
s.fitParams.startParams = [.1 .1 .05 1 s.fitParams.teStart .5]; % K1 k2 vB BP te occ
s.fitParams.lowBound = [0 0 0 0 0 0]; 
s.fitParams.upBound = [inf inf inf inf s.scanDur-5-s.T 1];
[model,occCurve,discOcc,Ks,exitFlag] = multiStepModel_1TCM(s);
% save the output somehow 