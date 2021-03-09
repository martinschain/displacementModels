function[Cnd,Cs,model,occCurve,discOcc,Ks,exitFlag] = multiStepModel_2TCM(s)
% Applies the 2TC multistep displacement model
% Input data should be organized in structure s (see ReadMe.txt)

% Output arguments:
% Cnd, Cs and model: displacement model curves for ND, S and total binding
% occCurve: the time course for the occupancy
% discOcc: the discretized occupancy curve used in the model
% Ks: Vector of parameters: K1 k2 k3 k4 vB te occ (te is the end-of-growth time) 
% Note: To get the actual te, do te = s.T+Ks(5);
% exitFlag: the exitFlag from the fit, see lsqnonlin help section
%_____________________________________________________________________
%                      Martin Schain, Neurobiology Research Unit, 2021 

p0 = s.fitParams.startParams; 
lb = s.fitParams.lowBound;
ub = s.fitParams.upBound;

[K,~,~,exitFlag] = lsqnonlin(@fitMultiStepModel_2TCM,p0,lb,ub,s.fitParams.options,s);
[Cnd,Cs,model,occCurve,discOcc] = multiStepModel_2TCM_createCurves(K,s);
Ks = K(:);

function[err] = fitMultiStepModel_2TCM(p,s)
    [~,~,modelCurve] = multiStepModel_2TCM_createCurves(p,s);
    err = (interp1(s.t,modelCurve,s.tPET) - s.TAC).*s.weights;

    
function[CndOut,CsOut,modelCurve,o_curve,fullCurve] = multiStepModel_2TCM_createCurves(p,s)

    inFcn = s.inFcn;
    stepSize = s.stepSize;
    vB = p(5);
    wb = s.wb;
    te = s.T + p(6);
    tm = s.T + .5*p(6); % tm is set equidistant between tb and te (symmetric growth)
    o = p(7);
    t = s.t;
    tb = s.T;
    o_curve = occCurve(t,tb,tm,te,o);

    try
        [fullCurve,T_nodes,delta] = discretizeOccCurve(s,o_curve,te,0);
    catch
        error('Error in discretizing the sigmoid')
    end
    delta = 1-delta; %holds 1-occupancy at each node (i.e., the term to be multiplied to k3)  

    tt = length(T_nodes);
    for j = 1:tt 
        [~,k] = min(abs(s.t-T_nodes(j)));
        Tid(j) = k; % holds time indices for the jumps
    end
    Tid = unique(Tid);

    %% Pre intervention
    [K1,k2,k3,k4,d,t1,t2,p1,p2] = setParameters(p,1);
    t_cur = s.t(1:Tid(1));
    Hnd = K1/d * ( (t1-k4)*exp(-t1*t_cur) - (t2-k4)*exp(-t2*t_cur));
    Hs  = K1*k3/d *(-exp(-t1*t_cur) + exp(-t2*t_cur));
    CndPre = stepSize*filter(inFcn(1:Tid(1)),1,Hnd);
    CsPre = stepSize*filter(inFcn(1:Tid(1)),1,Hs); 
    CndOut = CndPre;
    CsOut = CsPre;

    %% Post intervention
    CndT = CndPre(end);
    CsT = CsPre(end);

    for m = 1:length(Tid) % loop through the nodes
        if m < length(Tid) % If we are at intermediate time intervals
            try
            timeIDs = Tid(m)+1:Tid(m+1);
            catch
                disp('asd')
            end
        elseif m == length(Tid) % if we re at the last interval
            timeIDs = Tid(m)+1:length(s.t);
        end
        t_cur = s.t(timeIDs);
        [K1,k2,k3,k4,d,t1,t2,p1,p2] = setParameters(p,delta(m));
        tau = t_cur-T_nodes(m);
        H = K1*( (t1-k4)/d*exp(-t1*tau) + ((t2-k4)/(-d))*exp(-t2*tau) );
        a1 = stepSize*filter(inFcn(timeIDs),1,H);
        a2 = ( (t1*CndT-k4*(CndT+CsT))/d )*exp(-t1*tau);
        a3 = ( (t2*CndT-k4*(CndT+CsT))/(-d) )*exp(-t2*tau);
        fnd = (a1+a2+a3);

        H = K1*k3*(exp(-t1*tau)-exp(-t2*tau))/-d;
        aa1 = stepSize*filter(inFcn(timeIDs),1,H);
        aa2 = ( (-k3*CndT+(k4-t2)*CsT)/d )*exp(-t1*tau);
        aa3 = ( (-k3*CndT+(k4-t1)*CsT)/(-d) )*exp(-t2*tau);
        fs = aa1+aa2+aa3;
        try
            CndT = fnd(length(fnd));
        catch
            disp('Skipped a step (Cnd)') % this can happen if a time-interval gets too short
            break
        end

        try 
            CsT = fs(length(fs));
        catch
            disp('Skipped a step (Cs)') % this can happen if a time-interval gets too short
            break
        end

        CndOut = [CndOut(:);fnd(:)];
        CsOut = [CsOut(:);fs(:)];
    end

    modelCurve = (1-vB)*(CndOut+CsOut) + vB*wb(:);

function[K1,k2,k3,k4,d,t1,t2,p1,p2] = setParameters(p,delta)
    K1 = p(1);
    k2 = p(2);
    k3 = delta*p(3);
    k4 = p(4);

    d = sqrt((k2+k3+k4)^2 - 4*k2*k4);
    t1 = (k2+k3+k4+d)/2;
    t2 = (k2+k3+k4-d)/2;
    p1 = K1*(t1-k3-k4)/d;
    p2 = K1*(t2-k3-k4)/(-d);