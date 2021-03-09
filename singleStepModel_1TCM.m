function[model,occCurve,Ks,exitFlag] = singleStepModel_1TCM(s)
% Applies the 1TC single step displacement model
% Input data should be organized in structure s (see ReadMe.txt)

% Output arguments:
% Model: displacement model curve
% occCurve: the time course for the occupancy (here its just a jump)
% Ks: Vector of parameters: K1 k2 vB BP occ T (T is when the jump occurs)
% exitFlag: the exitFlag from the fit, see lsqnonlin help section
%_____________________________________________________________________
%                      Martin Schain, Neurobiology Research Unit, 2021                                

p0 = s.fitParams.startParams;
ub = s.fitParams.upBound;
lb = s.fitParams.lowBound;

[K,~,~,exitFlag] = lsqnonlin(@fitSingleStep_1TCM,p0,lb,ub,s.fitParams.options,s);
model = singleStep_1TCM_createCurves(K,s);
occCurve = stepfun(s.t,K(6))*occ;
Ks = K(:);


function[err] = fitSingleStep_1TCM(K,s)

    model = singleStep_1TCM_createCurves(K,s);
    err   = ( interp1(s.t,model(:),s.tPET,'pchip') - s.TAC(:)).*s.weights(:);
    
    
function[model] = singleStep_1TCM_createCurves(K,s)

    inFcn = s.inFcn;
    t = s.t;
    stepSize = s.stepSize;
    wb = s.wb;
    T = K(end);
    [~,Tid] = min(abs(t-T));
    t_pre = t(1:Tid);
    t_post = t(Tid+1:end);
    vB = K(3);
    
    %% Pre T

    irfPre = K(1)*exp(-K(2)*t_pre / (1+K(4)) );
    CtPre = stepSize*filter(inFcn(1:Tid),1,irfPre);
  
    
    %% Post T
    
    tau = t_post-T;
    CtT = CtPre(end);
    a =  1+(1-K(5))*K(4);
    irfPost = K(1)*exp(-K(2)*tau /a);
    CtPost = stepSize*filter(inFcn(Tid+1:end),1,irfPost) + CtT*exp(-K(2)*tau/a);
    
    Ct = [CtPre(:);CtPost(:)];
    model = (1-vB)*Ct + vB*wb(:);



    

