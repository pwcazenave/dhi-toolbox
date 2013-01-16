function data=dfs02csv(filename,outdir,spec)
% DATA=DFS02CSV(FILENAME) converts a timeseries model output to CSV.
%   DFS02CSV(FILENAME) converts said dfs0 file to a csv file in the same
%   directory.
%   DFS02CSV(FILENAME,OUTDIR) converts the dfs0 and outputs the CSV file to
%   OUTDIR.
%   DFS02CSV(FILENAME,OUTDIR,SPEC) converts the dfs0 and outputs the CSV
%   file to OUTDIR with specified output format (C-style format).
%
%   For example, to export to the same directory as the input with an
%   alternative format specifier, run:
%
%       dfs02csv(filename,[],'%.2f')
%
%   No other conversions are made and all data are stored as doubles
%   internally. Default output format is ten decimal places (%.10f).
%
%   Pierre Cazenave v1.0 2011/06/02

% Let's do this...

% Get the output directory
filesplit=regexp(regexprep(filename,'\','/'),'/','split');
if nargin==1 || isempty(outdir)
    outdir=[];
    for i=1:size(filesplit,2)-1
        outdir=cat(2,outdir,char(filesplit(i)));
    end
    spec='%.10f';
elseif nargin==2
    spec='%.10f';
elseif nargin>3
    error('Please supply the correct number of arguments.')
end

% Get the basename. If the last five characters aren't '.dfs0', quit.
if strcmpi(filename(end-4:end),'.dfs0')~=1
    error('Incorrect file extension.')
end

basename=[filesplit{end}(1:end-5),'.csv'];

outFile=[outdir,'/',basename];

% Right, read in the data.

% Use the dfsTSO version.
dfs0=dfsTSO(filename);
t=readTimes(dfs0);
itemnames=get(dfs0,'itemnames');

data=nan(size(t,1),size(itemnames,1));
for i=1:size(itemnames,1)
    data(:,i)=dfs0(i);
end

% Use savecsv to write the output CSV file.
savecsv(outFile,[t,data],[],spec);
