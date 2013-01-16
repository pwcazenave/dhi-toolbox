function [elemAreaDeg,elemAreaMet,sepXYMetres,sepXYMetresMine,elemDist]=dsfuStats(filename)
% DFSUSTATS calculated the element area in degrees and metres. It also
%   obtains the distance between every node in metres. Finally, distances
%   between element centres are calculated in degrees.
%
%   [ELEMAREADEG,ELEMAREAMET,SEPXYMETRES,SEPXYMETRESMIN,ELEMDIST]=DFSUSTATS(FILENAME)
%   reads in the dfsu file FILENAME and calcualtes the element area in
%   degrees, the area in metres, the separation between nodes in metres
%   (with an alternative formulation) and the distance between element
%   centres in degrees respectively.
%
%   Pierre Cazenave 21/10/2010 v1.0
%

%  Let's get to it...

dfsu=dfsManager(filename);
tn=readElmtNodeConnectivity(dfsu);
% Node coordinates
[xn,yn,zn]=readNodes(dfsu); %#ok<NASGU>
% Element center coordinates
[xe,ye,ze]=readElmts(dfsu); %#ok<NASGU>
X(:,1)=xn(tn(:,1));
X(:,2)=xn(tn(:,2));
X(:,3)=xn(tn(:,3));
Y(:,1)=yn(tn(:,1));
Y(:,2)=yn(tn(:,2));
Y(:,3)=yn(tn(:,3));

% calculate the distance between each node
distX(:,1)=abs(X(:,1)-X(:,2));
distY(:,1)=abs(Y(:,1)-Y(:,2));
distX(:,2)=abs(X(:,2)-X(:,3));
distY(:,2)=abs(Y(:,2)-Y(:,3));
distX(:,3)=abs(X(:,3)-X(:,1));
distY(:,3)=abs(Y(:,3)-Y(:,1));

sepXY(:,1)=sqrt(distX(:,1).^2+distY(:,2).^2);
sepXY(:,2)=sqrt(distX(:,2).^2+distY(:,3).^2);
sepXY(:,3)=sqrt(distX(:,3).^2+distY(:,1).^2);

% Get (approximate) metres
sepXYMetres=distdim(sepXY(:),'deg','metres');

% Get (more accurate) metres
a=6378137; % radii of the earth in metres
b=6356752.3;
sepXYMetresMine=nan(size(sepXY));
for j=1:size(sepXY,1)
    for k=1:size(sepXY,2)
        sepXYMetresMine(j,k)=sepXY(j,k)*(((pi/180)*cosd(ye(j)))*sqrt(((a^4*cosd(ye(j)).^2)+(b^4*sind(ye(j)).^2))/((a*cosd(ye(j))).^2+(b*sind(ye(j))).^2)));
    end
end

clear a b j k Num

% Use the 3 distances between each node to calculate the area of the
% element using Heron's formula (numerically stable version) where a>=b>=c
workingDataDeg=sort(sepXY,2,'descend');
a=workingDataDeg(:,1);
b=workingDataDeg(:,2);
c=workingDataDeg(:,3);
elemAreaDeg=0.25*sqrt((a+(b+c)).*(c-(a-b)).*(c+(a-b)).*(a+(b-c)));
% Numerically unstable formula:
% s=(a+b+c)/2; % semiperimeter
% elemAreaDeg=sqrt(s*(s-a)*(s-b)*(s-c))
% clear s
clear workingDataDeg a b c
workingDataMet=sort(sepXYMetresMine,2,'descend'); % Not going to bother with distdim conversion
a=workingDataMet(:,1);
b=workingDataMet(:,2);
c=workingDataMet(:,3);
elemAreaMet=0.25*sqrt((a+(b+c)).*(c-(a-b)).*(c+(a-b)).*(a+(b-c)));
clear workingDataMet a b c

% Let's figure out the distance between element centres. Based on Adrian's
% code for determining the closest node between two meshes.
ref=nan(size(xe));
elemDist=nan(size(xe));

for n=1:length(xe)
    dist=sqrt(((xe-xe(n)).^2)+((ye-ye(n)).^2));
    % Use the first distance if there is more than one as they're
    % identical in distance.
    ref(n)=find(dist==min(dist(dist~=0)),1);
    elemDist(n)=min(dist(dist~=0));
end

