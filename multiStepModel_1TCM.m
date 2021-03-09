function[model,occCurve,discOcc,Ks,exitFlag] = multiStepModel_1TCM(s)
% Applies the 1TC multistep displacement model
% Input data should be organized in structure s (see ReadMe.txt)

% Output arguments:
% Model: displacement model curve
% occCurve: the time course for the occupancy
% discOcc: the discretized occupancy curve used in the model
% Ks: Vector of parameters: K1 k2 vB BP te occ (te is the end-of-growth time) 
% Note: To get the actual te, do te = s.T+Ks(5);
% exitFlag: the exitFlag from the fit, see lsqnonlin help section
%_____________________________________________________________________
%                      Martin Schain, Neurobiology Research Unit, 2021  

p0 = s.fitParams.startParams; 
lb = s.fitParams.lowBound;
ub = s.fitParams.upBound;

[K,~,~,exitFlag] = lsqnonlin(@fitMultiStepModel_1TCM,p0,lb,ub,...
        s.fitParams.options,s);
    
[model,occCurve,discOcc] = multiStepModel_1TCM_createCurves(K,s);
Ks = K(:);


function[err] = fitMultiStepModel_1TCM(p,s)
    modelCurve = multiStepModel_1TCM_createCurves(p,s);
    err = (interp1(s.t,modelCurve,s.tPET) - s.TAC).*s.weights;

function[modelCurve,o_curve,fullCurve] = multiStepModel_1TCM_createCurves(p,s)


    inFcn = s.inFcn;
    stepSize = s.stepSize;
    vB = p(3);
    wb = s.wb;
    te = s.T + p(5);
    tm = s.T + .5*p(5); % This enforces symmetric growth
    o = p(6);
    t = s.t;
    tb = s.T;
    o_curve = occCurve(t,tb,tm,te,o);

    try
        [fullCurve,T_nodes,delta] = discretizeOccCurve(s,o_curve,te,0);
    catch
        error('Error in discretizing the sigmoid')
    end
    delta = 1-delta; %holds 1-occupancy at each node (i.e. the term to be multiplied to k3) 

    tt = length(T_nodes);
    for j = 1:tt 
        [~,k] = min(abs(s.t-T_nodes(j)));
        Tid(j) = k; % Holds the time indices for the nodes
    end
    Tid = unique(Tid);

    %% Pre intervention
    t_cur = s.t(1:Tid(1)); %time array from start until first jump
    irfPre = p(1)*exp(-p(2)*t_cur / (1+p(4)) );
    CtPre = stepSize*filter(inFcn(1:Tid),1,irfPre);
    CtOut = CtPre;

    %% Post intervention
    CtT = CtPre(end);

    for m = 1:length(Tid) % loop through the nodes
        if m < length(Tid) % If we are at intermediate time intervals
            timeIDs = Tid(m)+1:Tid(m+1); % timeIDs will hold all time indices for the current interval
        elseif m == length(Tid) % if we re at the last interval
            timeIDs = Tid(m)+1:length(s.t);
        end
        t_cur = s.t(timeIDs);
        tau = t_cur-T_nodes(m);
        a =  1+delta(m)*p(4); % 1+(1-occ(t))*BP 

        irf = p(1)*exp(-p(2)*tau /a);
        Ct = stepSize*filter(inFcn(Tid+1:end),1,irf) + CtT*exp(-p(2)*tau/a);
        CtT = Ct(end);
        CtOut = [CtOut(:);Ct(:)];
    end

    modelCurve = (1-vB)*CtOut + vB*wb(:);


