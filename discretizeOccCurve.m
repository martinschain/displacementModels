function[fullCurve,tNew,oNew] = discretizeOccCurve(s,occCurve,te,doPlot)
% This function reads in an occupancy curve (generated with occCurve.m) and
% dicretizes it into n number of discrete steps. The value for n is stored
% in structure s. This routine is called inside an optimizer, thus the
% value for te me be changed in every iteration, and is therefore provided 
% as an input argument. 
%_____________________________________________________________________
%                      Martin Schain, Neurobiology Research Unit, 2021

nbrSteps = s.nbrSteps;
[~,T1id] = min(abs(s.t-s.T));
[~,teid] = min(abs(s.t-te));
if teid <= T1id
    teid = T1id+1;
end

dur = s.t(teid) - s.t(T1id);
stepDuration = dur/nbrSteps;
tNew = s.t(T1id):stepDuration:s.t(teid);
oNew = interp1(s.t,occCurve,tNew,'next');
fullCurve = interp1(tNew,oNew,s.t,'next');

if doPlot
    figure,
    plot(s.t,s.occCurve,'-'), hold on
    plot(tNew,oNew,'ro')
    plot(s.t,fc,'r--')
end

