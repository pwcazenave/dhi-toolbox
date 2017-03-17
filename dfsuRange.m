function [xn,yn,tn,zmin,zmax]=dfsuRange(path,filename,verbose)
% [XN,YN,TE,ZMIN,ZMAX]=DFSURANGE(PATH,FILENAME) takes a given dfsu
% file and extracts the minimum and maximum surface elevation for each
% node. 
%
% Intermediate files are saved in a subdirectory within the dfsu's
% directory.
%
% XN, YN and TN are the node locations and the element centres.
% 
% Optionally supply a VERBOSE argument ('yes', 'no', 'quiet' or 'noisy') to
% suppress output messages.

% Pierre Cazenave   v1.0 21/02/2012
%                       The way this is written at the moment means if your
%                       time series is long or your time step is short, you
%                       might run out of memory. Probably best to resample
%                       your output (temporally) to get it to fit in
%                       memory.

% Let's get it on...

% Set up the new MIKE interface to dfs files
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;

if nargin~=2 && nargin~=3
    error('Please supply the correct number of arguments.')
end

if exist('verbose','var')
    if strcmpi(verbose,'no') || strcmpi(verbose,'quiet')
        benoisy=0;
    else
        benoisy=1;
    end
else
    benoisy=0;
end


suffix='_TR';

% If saved files exist, we'll load them preferentially
existed=0;
if exist([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'.mat'],'file')==2
    % The ZMIN and ZMAX data have already been extracted
    load([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'.mat']);
    existed=1;
end

if existed==0
    % At least one of the steps wasn't complete, so we need to work on the
    % dfsu file directly. So, load said dfsu file.
    dfsu=DfsFileFactory.DfsuFileOpen([path,filename]);

    % Check the intermediate file directory exists, and if not, make it.
    if exist([path,'/',filename(1:end-5),'_files/'],'dir')~=7
        mkdir([path,'/',filename(1:end-5),'_files/'])
    end

    % Get some basic information

    % Read some item information
    if benoisy
        fprintf('Extracting item information from %s.\n',filename)
    end
    items=cell(dfsu.ItemInfo.Count,3);
    for i=0:dfsu.ItemInfo.Count-1
       item=dfsu.ItemInfo.Item(i);
       items{i+1,1}=char(item.Name);
       items{i+1,2}=char(item.Quantity.Unit);
       items{i+1,3}=char(item.Quantity.UnitAbbreviation);
    end
    % Number of timesteps
    num_timesteps=double(dfsu.NumberOfTimeSteps);
    % Timestep in seconds
    timestep=double(dfsu.TimeStepInSeconds);
    % Start time
    start_time=double([dfsu.StartDateTime.Year,dfsu.StartDateTime.Month,dfsu.StartDateTime.Day,dfsu.StartDateTime.Hour,dfsu.StartDateTime.Minute,dfsu.StartDateTime.Second]);
    % End time
    end_time=datevec(datenum(start_time)+((num_timesteps-1)*(timestep/(24*60*60))));
    % Model base time i.e. start of model year (hope it doesn't cross the
    % end of the year).
    base_time=[start_time(1),0,0,0,0,0];
    % Model time as decimal day of year
    time=datenum(start_time)-datenum(base_time)-1:timestep/(24*60*60):datenum(end_time)-datenum(base_time)-1; % Matlab format time

    % Identify the surface elevation item from the dfsu item list
    z_num=find(strcmp('Surface elevation',items));
    % Check whether we've got anything useful. If not, bomb out with an
    % error.
    if isempty(z_num)
        warning('MATLAB:paramAmbiguous','Unrecognised analysis type.\n')
        dfsu.Close();
        error('Check value of ''type'' or that the input file has correct output items (Surface elevation) and try again.')
    end

    % Create element table in Matlab format
    if benoisy
        fprintf('Extracting mesh from %s.\n',filename)
    end
    tn=mzNetFromElmtArray(dfsu.ElementTable);
    % Node coordinates
    xn=double(dfsu.X);
    yn=double(dfsu.Y);
    zn=double(dfsu.Z);
    % Element centres
    [xe,ye,ze]=mzCalcElmtCenterCoords(tn,xn,yn,zn); %#ok<ASGLU>
    
    % Find out how many iterations we've got, and if it's greater than a
    % fortnight's worth, clip it to a fortnight. Otherwise, it's the
    % model's duration.
    duration=timestep*num_timesteps;
    spring_neap_days=((12+(25/60))*28)/24; % 12h25m period for ~ 14 days.
    if duration>(60*60*24*spring_neap_days)
        % Analyse maximum of a spring-neap cycle (roughly a fortnight's
        % data).
        warning(['Truncating timeseries to ',num2str(spring_neap_days),' days.']) %#ok<WNTAG>
        num_proc_steps=round((60*60*24*spring_neap_days)/timestep);
    else
        num_proc_steps=num_timesteps;
    end

    % Check if we've already saved all the basic info to _raw.mat, and
    % if not, then save it.
    if exist([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'_raw.mat'],'file')~=2
        % Save just about everything
        try
            save([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'_raw.mat'],...
                'base_time','duration','spring_neap_days','end_time','items',...
                'num_proc_steps','num_timesteps',...
                'start_time','suffix','tn','time','timestep','tn',...
                'z_num','xe','xn','ye','yn','ze','zn');
        catch
            dfsu.Close();
            error('Unable to save intermediate ''raw'' MATLAB file.')
        end
    end

end

% Load the general info
if benoisy
    fprintf('Load the saved basic data (%s)... ',[filename(1:end-5),'_raw.mat'])
end
load([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'_raw.mat'])
if benoisy
    fprintf('done.\n')
end

if existed==0
    % Need to find the minimum and maximum surface elevation values for
    % each node.
    zmin=nan(dfsu.NumberOfElements,1,'single');
    zmax=nan(dfsu.NumberOfElements,1,'single');
    z=nan(dfsu.NumberOfElements,num_proc_steps,'single');

    % Make sure we omit the first three day's worth of values (on the safe
    % side).
    start_iter=345600/timestep;

    % Make an array of iterable values.
    iter=start_iter:(num_proc_steps)+start_iter;
    
    % How are we doing?
    if benoisy
        h=waitbar(0,'Extracting surface elevation data.','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)
    end

    try
        for n=1:size(iter,2)
            if benoisy
                if getappdata(h,'canceling')
                    break
                end
            end
            getVal=iter(n)-1; % DHI start from zero.
            % Get the current surface elevation
            try
                zCurr=double(dfsu.ReadItemTimeStep(z_num,getVal).Data);
            catch
                dfsu.Close();
                error('Can''t find the necessary surface elevation values. Aborting.')
            end
            % Save the current time step's surface elevation for later
            % comparison.
            z(:,n)=zCurr;
            if benoisy
                waitbar(n/size(iter,2),h);
            end
        end
        if benoisy
            delete(h)
        end
    catch
        if benoisy
            delete(h)
        end
    end
    % Now find the minimum and maximum values for each set of nodal surface
    % elevation values.
    zmin=min(z,[],2);
    zmax=max(z,[],2);

    % Save the results in the appropriate .mat file
    if benoisy
        fprintf('Saving tidal range results to disk... ')
    end
    save([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'.mat'],...
        'duration','num_proc_steps','zmin','zmax');
    if benoisy
        fprintf('done.\n')
    end
end

if existed==1
    % All the data have been loaded already, in which case, we're golden.
    if benoisy
        fprintf('All tidal range data have been loaded successfully.\n')
    end
end
