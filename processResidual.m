function [ur,vr,rmag,rdir]=processResidual(x,y,modelTime,timeStep)
% [UR,VR,RMAG,RDIR]=PROCESSRESIDUAL(X,Y,MODELTIME,TIMESTEP) takes the X and
% Y summed vector output from DFSURESIDUAL and calculates the last few
% steps necessary to get a magntidue and direction. Using TIMESTEP (the
% model output time step), a tide window for a spring-neap cycle is
% calculated and used to create a new residual vector. The components of
% this vector are stored in UR and VR; the magnitude and direction of those
% vectors are in RMAG and RDIR. MODELTIME is the model time series from the
% DFSURESIDUAL saved output.
%
% Pierre Cazenave 01/09/2011 v1.0.
%

% Time to do this...

% The duration of a tidal cycle (12h25m)
tideCycle=(12+(25/60))/24;
tideWindow=ceil(tideCycle*(24*60*60/timeStep));
% That duration as a number of seconds (use to calculate residual as m/s).
dt=(mean(modelTime(length(modelTime)-tideWindow:length(modelTime)),2)-mean(modelTime(1:tideWindow),2))*24*60*60;

% The vectors x and y contain the residual values (i.e. the Eulerian
% displacement). uRaw and vRaw contain the vectors, unsummed. Calculating
% the mean from a single tidal cycle gives you the midpoint of that cycle.
% Likewise for the last cycle. Then, joining the two positions up with a
% line and calculating its length is the Eurlerian distance a particle has
% travelled. Dividing that length by the duration (dt) gives the residual
% flow speed. The angle is just the inverse tangent of the two positions.
% The spring-neap tidal cycle is 14.4861 days.

% Get the mean position of the start of the tide.
uStart=mean(x(:,1:tideWindow),2);
vStart=mean(y(:,1:tideWindow),2);
uEnd=mean(x(:,size(x,2)-tideWindow:size(x,2)),2);
vEnd=mean(y(:,size(y,2)-tideWindow:size(y,2)),2);

% To get the vector components, calculate the difference betwee the start
% and the end.
ur=uEnd-uStart;
vr=vEnd-vStart;

% Get the direction from the vector components
rdir=atan2(ur,vr)*(180/pi);

% Now for the magnitude
rmag=sqrt(ur.^2+vr.^2);

% Scale the magnitude by the duration (dt) to get m/s
rmag=rmag/dt;
