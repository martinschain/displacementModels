function[Cnd,Cs,model,occCurve,Ks,exitFlag] = singleStepModel_2TCM(s)
% Applies the 2TC singlestep displacement model
% Input data should be organized in structure s (see ReadMe.txt)

% Output arguments:
% Cnd, Cs and Model: displacement model curves for ND, S and total binding 
% occCurve: the time course for the occupancy (here its just a jump)
% Ks: Vector of parameters: K1 k2 k3 k4 vB occ T (T is when the jump occurs)
% exitFlag: the exitFlag from the fit, see lsqnonlin help section
%_____________________________________________________________________
%                      Martin Schain, Neurobiology Research Unit, 2021


p0 = s.fitParams.startParams;
ub = s.fitParams.upBound;
lb = s.fitParams.lowBound;

[K,~,~,exitFlag] = lsqnonlin(@fitSingleStepModel_2TCM,p0,lb,ub,s.fitParams.options,s);
[Cnd,Cs,model] = singleStepModel_2TCM_createCurves(K,s);
occ = 1-K(6);
Ts = K(7);
occCurve = stepfun(s.t,Ts)*occ;
Ks = K;

function[err] = fitSingleStepModel_2TCM(K,s)

    [~,~,model] = singleStepModel_2TCM_createCurves(K,s);
    err  = ( interp1(s.t,model(:),s.tPET,'pchip') - s.TAC(:) ).*s.weights(:);
    
    
function[Cnd,Cs,model] = singleStepModel_2TCM_createCurves(K,s)
    
    inFcn = s.inFcn;
    t = s.t;
    stepSize = s.stepSize;
    wb = s.wb;    
    T = K(end);
    [~,Tid] = min(abs(t-T));
    t_pre = t(1:Tid);
    t_post = t(Tid+1:end);
    vB = K(5);
    
    %% Pre T
    [K1,k2,k3,k4,d,t1,t2,p1,p2] = setParameters(K(1:end-1),0); %set parameters before challenge
    Hnd = K1/d * ( (t1-k4)*exp(-t1*t_pre) - (t2-k4)*exp(-t2*t_pre));
    Hs  = K1*k3/d *(-exp(-t1*t_pre) + exp(-t2*t_pre));
    CndPre = stepSize*filter(inFcn(1:Tid),1,Hnd);
    CsPre = stepSize*filter(inFcn(1:Tid),1,Hs); 
    
    %% Post T
    [K1,k2,k3,k4,d,t1,t2,p1,p2] = setParameters(K(1:end-1),1); %set parameters after challenge
    tau = t_post-T;
    CndT = CndPre(end);
    CsT = CsPre(end);
    H = K1*( (t1-k4)/d*exp(-t1*tau) + ((t2-k4)/(-d))*exp(-t2*tau) );
    a1 = stepSize*filter(inFcn(Tid+1:end),1,H);
    a2 = ( (t1*CndT-k4*(CndT+CsT))/d )*exp(-t1*tau);
    a3 = ( (t2*CndT-k4*(CndT+CsT))/(-d) )*exp(-t2*tau);
    CndPost = (a1+a2+a3);

    H = K1*k3*(exp(-t1*tau)-exp(-t2*tau))/-d;
    aa1 = stepSize*filter(inFcn(Tid+1:end),1,H);
    aa2 = ( (-k3*CndT+(k4-t2)*CsT)/d )*exp(-t1*tau);
    aa3 = ( (-k3*CndT+(k4-t1)*CsT)/(-d) )*exp(-t2*tau);
    CsPost = aa1+aa2+aa3;
    
    Cnd = [CndPre(:); CndPost(:)];
    Cs = [CsPre(:); CsPost(:)];
    CtTot = Cnd+Cs;
    model = (1-vB)*CtTot(:) + vB*wb(:);

    
    
function[K1,k2,k3,k4,d,t1,t2,p1,p2] = setParameters(K,afterChallenge)
    K1 = K(1);
    k2 = K(2);
    k3 = K(3);
    k4 = K(4);
    delta = K(6); %delta is actually 1-occ,
    
    if afterChallenge
        k3 = delta*k3;
    end
    
    d = sqrt((k2+k3+k4)^2 - 4*k2*k4);
    t1 = (k2+k3+k4+d)/2;
    t2 = (k2+k3+k4-d)/2;
    p1 = K1*(t1-k3-k4)/d;
    p2 = K1*(t2-k3-k4)/(-d);

