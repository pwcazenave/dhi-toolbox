function status=exportdfs1(data,filename,starttime)
% EXPORTDFS1(DATA,FILENAME) takes a 2D array (DATA) with time and a number
% of points as follows:
%
%  t   va1 va2 va3 va4 va5 ...
%  ---------------------------
%  1   0.1 0.1 0.2 0.2 0.3 ...
%  2   0.2 0.2 0.3 0.3 0.4 ...
%  3   0.2 0.2 0.3 0.3 0.4 ...
%  4   0.1 0.1 0.2 0.2 0.3 ...
%
% and exports it formatted correctly for input into MIKE21 to FILENAME.
%
% EXPORTDFS1(DATA,FILENAME,STARTTIME) takes a 2D array (DATA) as exports as
% above, but sets the start time to STARTTIME (hh:mm:ss).
%

% Pierre Cazenave pwc101@soton.ac.uk v0.1 02/06/2010
% 24/08/2010 Added filename to the progress indicator.

if nargin<3
    starttime='01:00:00';
end

timeStep=floor((data(2,1)-data(1,1))*86400);

f=fopen(filename,'w');

% Header
fprintf(f,'Title\tE:\\modelling\\data\\boundary_files\\round_8_palaeo\\uehara_boundary\\dfs1\\bta0-00ka_east.dfs1\n');
% fprintf(f,['Time\t1990-01-01\t',starttime,'\t',num2str(size(data,1)),'\t',num2str(timeStep),'\n']);
fprintf(f,['Time\t1990-01-01\t','%s','\t','%i','\t','%i','\n'],starttime,size(data,1),timeStep);
fprintf(f,['NoGridPoints\t','%i','\n'],size(data,2)-1);
% From MIKE21 Help
% "Grid Step: The grid step for the data. The unit is meters."
% Using the uehara_v2.mesh file and generating dfs1 files from that, then
% looking this value in those files gives:
% 0.0481489 0.0458233 0.0461409 0.0641 0.0455633 0.0524667
% These all seem implausibly small, but I'm going to go with the median
% for the time being (0.0471449m).
fprintf(f,'Spacing	0.0471449\n');
fprintf(f,'NoDynamicItems\t1\n');
fprintf(f,'Item	Predicted\tWater Level	meter\n');
fprintf(f,'Delete\t-1E-030\n');
fprintf(f,'\n');

% Get the output filename from FILENAME
outname=regexpi(regexprep(filename,'_',' '),'\','split');
outname=outname(end);

% All the data
h=waitbar(0,['Processing ',num2str(size(data,1)),' timesteps'],'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)
for i=1:size(data,1);
    fprintf(f,'tstep\t%i\n',i-1);
    if getappdata(h,'canceling')
        break
    end
    for k=1:size(data(i,2:end),2);
        fprintf(f,[num2str(data(i,k+1)),'\t\n']);
    end
    waitbar(i/size(data,1),h,sprintf('%s: %i of %i',char(outname),i,size(data,1)))
    fprintf(f,'\n');
end
delete(h)
status=fclose(f);
