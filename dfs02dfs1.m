% function [exitstatus]=dfs02dfs1(indir,outdir,prefix)
% Function to generate a dfs1 file from a series of dfs0 files.
% Originally scripted by David Lambkin.
% Modified to a function by Pierre Cazenave 07/06/2010 pwc101@soton.ac.uk
%
% DFS02DFS1 takes a INDIR with

clear; close all; clc
base=checkos;
indir=[base,'\modelling\data\boundary_files\round_8_palaeo\uehara_boundary\dfs0\'];
outdir=('C:\Users\goosin\Desktop\test');
prefix='bta0-00ka_north';

filenames=dir([indir,'/*.dfs0']);%{'N1','N2','N3','N4','N5','N6','N7','N8','N9'};
data=[];

for n=1:length(filenames);
    filename=[indir,filenames(n).name];
    dfs0=dfsTSO(filename);
    data(:,n)=dfs0(1);

    if n==1
        t=readTimes(dfs0);
        startDate=[t(1,1),t(1,2),t(1,3)];
        t2=datenum(t)-datenum(startDate); % Take out absolute dates and have only relative dates
        timestep=floor((t2(2,1)-t2(1,1))*86400);
%         timestep=t(2,6)-t(1,6)+(60*(t(2,5)-t(1,5)))+(60*60*(t(2,4)-t(1,4)));
    end
    close (dfs0);
end

% data=fliplr(data); % Reverse order of columns

% data=(round(data*1000))/1000; % round water levels to nearest millimeter

% filename=[indir,'SETP_North_boundary_',num2str(year),'.dfs1'];
dfs=dfsManager([outdir,'\',prefix,'.dfs1'],1);
createTemporalDef(dfs,'Equidistant_Calendar',t(1,:),size(t,1),timestep);
setSpatialDef(dfs,'LONG/LAT',[0,0],0,[0,1,length(filenames)]);
addItem(dfs,'Surface height','Surface Elevation','m');
create(dfs);

h=waitbar(0,'Writing North boundary data (stage 1/4), progress...');
for i=0:size(t,1)-1
  dfs(1,i)=data(i+1,:);
  waitbar(i/size(t,1),h);
end
close(h)

saveAndClose(dfs);
fprintf('\nFile created: ''%s''\n',[path,filename]);