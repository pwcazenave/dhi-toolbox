function Z=transAsymmetry(u,v)
% Z=TRANSASYMMETRY(U,V) takes the output from DFSURESIDUAL (the U and
% V vector components) and calculates corresponding transport asymmetry,
% where that is calculated as follows:
%
%   1. Calculate the linear regression for the eastward and northward
%   transport vectors (U and V).
%   2. Remove the linear trend to produce a gradient of zero.
%   3. Calculate the vector magnitude (maintaining sign).
%   4. Separate the magnitudes greater and less than zero into two separate
%   arrays (internally stored as U1 and U2).
%   5. Using U1 and U2, calculate the transport asymmetry (R) as:
%       R=U1/(U1+U2)
%   6. R is limited to 0.5<=R<=1 (if R>0.5, R=1-R).
%
% This approach is described in Dix, J.K., Lambkin, D.O. and Cazenave, P.W.
% (2007) Development of a Regional Sediment Mobility Model for Submerged
% Archaeological Sites. English Heritage ALSF project no. 5524. School of
% Ocean and Earth Science, University of Southampton, U.K. 156pp.
%
% Based on David Lambkin's code in dfsuasymmetrySED.m.
%
% Pierre Cazenave 2011/09/01 v1.0.
%                 2011/09/05 v1.1. Fixed polyfit warning message by
%                 outputting all three outputs, but discarding the two I'm
%                 not interested in.
%                 2011/10/18 v1.2. Added step 6 and improved speed by using
%                 logical indexing for the reciprocal and step 6.
%

% Let's do this...

% Preallocate arrays to speed this up a bit.
p_array=nan(size(u,1),2);
U1=nan(size(u,1),1);
U2=nan(size(u,1),1);
uRot=nan(size(u));
vRot=nan(size(u));

% Find the gradient of the transport axis
for n=1:size(u,1)
    if max(u(n,:)>0)
        [p,~,~]=polyfit(u(n,:),v(n,:),1);
        p_array(n,:)=p;
    else
        p_array(n,:)=[nan,0];
    end
end

p_array(:,2)=[];
p_array=atan(p_array); % converts the slope of the data axis into radians

% Rotate the data so it is parallel to the X axis
for n=1:size(u,2)
    uRot(:,n)=(cos(-p_array).*u(:,n))-(sin(-p_array).*v(:,n));
    vRot(:,n)=(sin(-p_array).*u(:,n))+(cos(-p_array).*v(:,n));
end

% Find the magnitude of the transport
mag_data=sqrt(uRot.^2+vRot.^2).*sign(uRot);

% Find the sum of transport in each direction, maintaining the sign
for n=1:size(mag_data)
    ref=find(mag_data(n,:)>0);
    U1(n)=sum(mag_data(n,ref));
    ref=find(mag_data(n,:)<0);
    U2(n)=abs(sum(mag_data(n,ref)));
end

% Calculate asymmetry of the signal
Z=U1./(U1+U2);

% for n=1:length(Z)
%     if Z(n)>1
%         Z(n)=1/Z(n);
%     end
% end
% Z(Z>1)=1/Z(Z>1);
Z(Z<0.5)=1-Z(Z<0.5);

% Make the output compatible with the X and Y arrays for use with patch().
Z=Z';

