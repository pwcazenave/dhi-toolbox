function bodc2dfs0(infiles,indir,outdir,prefix)
% ASCII2DFS0 Convert BODC tidal data into dsf0 files.
%   ASCII2DFS0(INFILES,INDIR,OUTDIR,PREFIX) takes a file INDIR/INFILES and
%   writes surface elevation to OUTDIR/PREFIX.dfs0
%
%   INFILES must be a struct of the input filenames (most easily
%   accomplished with a dir command).
%
%   Each INFILE must be formatted:
%       line# year month day hour minute second tide residual
%
% Pierre Cazenave 14/01/2011 v1.0

% Let's get it on...

if nargin~=4
    error('Incorrect number of inputs')
end

% Read in the data

% Not preallocating because the size of the input files varies as a
% function of the year (sampling changed from hourly to every 15 minutes
% some time in the late nineties). Sucking the data off the disk is slower
% than not preallocating, so there's no massive performance loss here
% anyway.
for i=1:size(infiles,1)
    data(:,1:9,i)=load([indir,'\',infiles(i).name]); %#ok<AGROW>
end

% Start date from the data
startDate=data(1,2:7,1);

% Create the dfs0 output file in outdir with prefix as the name prefix.
outfile=[outdir,'\',prefix,'.dfs0'];
dfs0=dfsTSO(outfile,1);
% set(dfs0,'filetitle','Average of water depth and flux');
set(dfs0,'startdate',startDate);
set(dfs0,'timestep',data(2,2:7,1)-data(1,2:7,1));
set(dfs0,'deletevalue',-99);
addTimesteps(dfs0,size(data(:,1,1),1));
for j=1:size(data,3)
    site=infiles(j).name;
    addItem(dfs0,[site(5:7),': Surface elevation'],'Surface elevation','meter');
end

% Add the observed data to the dfs0 object.
for n=1:size(data,3) % each site
    dfs0(n)=single(data(:,8,n));
end

% for n=1:size(data,3) % each site
%     for m=1:size(data,1) % each timestep
%         % Convert to single to avoid warnings.
%         dfs0(m,n)=single(data(m,8,n));
%     end
% end

save(dfs0)
close(dfs0)

% fprintf('\nFile created: ''%s''\n',outfile);
