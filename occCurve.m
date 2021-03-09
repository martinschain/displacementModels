function[o_curve] = occCurve(t,tb,tm,te,o)
% This function generates an occupancy growth model based on the input
% parameters.
% t: time vector
% tb: begin time of growth
% tm: time of max derivative
% te: end time of growth
% o: maximal occupancy reached
%_____________________________________________________________________
%                      Martin Schain, Neurobiology Research Unit, 2021 
    
f = stepfun(t,tb) .* (1+(te-t)./(te-tm)).*((t-tb)./(te-tb)).^((te-tb)/(te-tm));
y = f + stepfun(t,te).*(1-f);
[~,teid] = min(abs(t-te));
y(teid:end) = y(teid);
o_curve = o*y;
    