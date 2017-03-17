function [xn,yn,tn,u,v,x,y]=dfsuResidual(path,filename,type)
% [XN,YN,TE,U,V,X,Y]=DFSURESIDUAL(PATH,FILENAME,TYPE) takes a given dfsu
% file and extracts the U and V components. These are used to calculate
% the residual flow and transport directions for a given model run, stored
% in X and Y.
%
% Intermediate files are saved in a subdirectory within the dfsu's
% directory.
%
% The string TYPE can be one of:
%   'total' - Total sediment transport
%   'sus' - Suspended load transport
%   'bed' - Bedload transport
%   'flow' - Tidal flow
%
% XN, YN and TN are the node locations and the element centres.


% Based on Dave Lambkin's residual analysis code.
% Pierre Cazenave   v1.0 20/09/2010
%                   v1.1 12/10/2010 Updated with more robust item checks,
%                   and make sure we were correctly handling the
%                   subsampling of the times.
%                   v1.2 09/05/2011 Fixed the range of iter to ensure it
%                   can't go beyond the last timestep in the input filename
%                   because of the skip of the first three days' data.
%                   Added a warning if we're starting from the first
%                   timestep.
%                   v1.3 26/05/2011 Removed all the 32/64 bit checks since
%                   the latest release of the DHI MATLAB tools work in 64
%                   bit natively. Also removed the summation code when
%                   skipping data. Now, the value of timestep is increased
%                   to reflect the number of data points skipped. Thus,
%                   when the residual is calculated, the correct duration
%                   (based on the new timestep value) is used. Finally,
%                   converted everything to the new DHI MATLAB tools
%                   interface, which, incidentally, has made everything
%                   MUCH faster.
%                   v1.3.1 29/07/2011 Added try/catch statements so if
%                   reading in the dfsu file fails for some reason, the
%                   progress bars close, eliminating the requirement to
%                   kill MATLAB to exit it.
%                   v1.3.2 21/02/2012 Moved the code creating the node
%                   tables to after the check for the right type of
%                   analysis so the code bombs out if we've been given the
%                   wrong file before it slowly reads the node and element
%                   info. Also fixed an issue where the u and v vector data
%                   weren't being read from the start of the iter array but
%                   instead from the first value in the dfsu file. This is
%                   problematic if you didn't skip the first couple of
%                   tidal cycles in the creation of the output file as the
%                   initial tidal surge would contribute to the overall
%                   residual calculation. 
%
%

% Let's get it on...

% Set up the new MIKE interface to dfs files
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;

if nargin~=3
    error('Please supply the correct number of arguments.')
end

% Generate appropriate suffix for the output files.
if strcmpi(type,'total')==1 % Total sediment transport
    suffix='_TL';
elseif strcmpi(type,'sus')==1 % Suspended load transport
    suffix='_SL';
elseif strcmpi(type,'bed')==1 % Bedload transport
    suffix='_BL';
elseif strcmpi(type,'flow')==1 % Tidal residual
    suffix='_UV';
else
    warning('MATLAB:paramAmbiguous','Unrecognised analysis type: %s\n',type)
    error('Check value of ''type'' and try again.')
end

% If saved files exist, we'll load them preferentially
existed=[0,0];
if exist([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'.mat'],'file')==2
    % The U and V data have already been extracted
    load([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'.mat']);
    existed(1)=1;
end
if exist([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'res.mat'],'file')==2
    % The residuals have already been calculated and saved
    load([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'res.mat']);
    existed(2)=1;
end

if sum(existed)~=2;
    % At least one of the steps wasn't complete, so we need to work on the
    % dfsu file directly. So, load said dfsu file.
    dfsu=DfsFileFactory.DfsuFileOpen([path,filename]);

    % Check the intermediate file directory exists, and if not, make it.
    if exist([path,'/',filename(1:end-5),'_files/'],'dir')~=7
        mkdir([path,'/',filename(1:end-5),'_files/'])
    end

    % Initialise the magDir variable (used if we don't have U and V
    % vector componenets.
    magDir=0;

    % Read some item information
    fprintf('Extracting item information from %s.\n',filename)
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

    % Identify the relevant vectors from the dfsu item list based on the
    % TYPE of analysis we're doing
    if strcmpi(type,'total')==1 % Total sediment transport
        u_num=find(strcmp('Total load - x-component',items));
        v_num=find(strcmp('Total load - y-component',items));
        % Check for direction and magnitude, from which we can make U
        % and V.
        if isempty(u_num) || isempty(v_num)
            magDir=1;
            u_num=find(strcmp('Suspended load - magnitude',items));
            v_num=find(strcmp('Suspended load - direction',items));
        end
    elseif strcmpi(type,'sus')==1 % Suspended load transport
        u_num=find(strcmp('Suspended load - x-component',items));
        v_num=find(strcmp('Suspended load - y-component',items));
        % Check for direction and magnitude, from which we can make U
        % and V.
        if isempty(u_num) || isempty(v_num)
            magDir=1;
            u_num=find(strcmp('Suspended load - magnitude',items));
            v_num=find(strcmp('Suspended load - direction',items));
        end
    elseif strcmpi(type,'bed')==1 % Bedload transport
        u_num=find(strcmp('Bed load - x-component',items));
        v_num=find(strcmp('Bed load - y-component',items));
        % Check for direction and magnitude, from which we can make U
        % and V.
        if isempty(u_num) || isempty(v_num)
            magDir=1;
            u_num=find(strcmp('Suspended load - magnitude',items));
            v_num=find(strcmp('Suspended load - direction',items));
        end
    elseif strcmpi(type,'flow')==1 % Tidal residual
        u_num=find(strcmp('U velocity',items(:,1)),1);
        v_num=find(strcmp('V velocity',items(:,1)),1);
        if isempty(u_num) || isempty(v_num)
            % Try the depth averaged values.
            u_num=find(strcmp('Depth average U velocity',items(:,1)),1);
            v_num=find(strcmp('Depth average V velocity',items(:,1)),1);
                % Check for direction and magnitude, from which we can
                % make U and V.
                if isempty(u_num) || isempty(v_num)
                    magDir=1;
                    u_num=find(strcmp('Suspended load - magnitude',items));
                    v_num=find(strcmp('Suspended load - direction',items));
                end
        end
    end
    % Check whether we've got anything useful. If not, bomb out with an
    % error.
    if isempty(u_num) || isempty(v_num)
        warning('MATLAB:paramAmbiguous','Unrecognised analysis type: %s\n',type)
        dfsu.Close();
        error('Check value of ''type'' or that the input file has correct output items (U and V vector components) and try again.')
    end

    % Create element table in Matlab format
    fprintf('Extracting mesh from %s.\n',filename)
    tn=mzNetFromElmtArray(dfsu.ElementTable);
    % Node coordinates
    xn=double(dfsu.X);
    yn=double(dfsu.Y);
    zn=double(dfsu.Z);
    % Element centres
    [xe,ye,ze]=mzCalcElmtCenterCoords(tn,xn,yn,zn); %#ok<ASGLU>
    
    % Make sure we give a warning about calculating U and V rather than
    % reading them from the input file.
    if magDir==1
        warning('Calculating U and V from magnitude and direction rather than reading from the output file.') %#ok<WNTAG>
    end

    % Thin the input data to be quarter hourly, rather than every five
    % minutes. This should alleviate some of the memory problems, and
    % doesn't require using MIKE to thin and create new dfsu files. We
    % also need to update the timestep value to reflect the new
    % sampling interval as well the the time vector.
    ssample=60*15;
    if timestep<ssample
        skipfac=ssample/timestep;
        % timestep=ssample;
        time=time(1:skipfac:end); %#ok<*NASGU>
    else
        skipfac=1;
    end

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

    % Amend the number of processing steps to account for the half
    % hourly adjustment.
    if timestep<ssample
        warning('Sampling every %.1f minutes, instead of %.1f minutes.',ssample/60,timestep/60) %#ok<WNTAG>
        num_proc_steps=round(num_proc_steps/skipfac);
        % We need to adjust the value of timestep to accommodate this
        % skip when calculating the residual. Simply multiply timestep
        % by the skip factor.
        timestep=timestep*skipfac;
    end

    % Check if we've already saved all the basic info to _raw.mat, and
    % if not, then save it.
    if exist([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'_raw.mat'],'file')~=2
        % Save just about everything
        try
            save([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'_raw.mat'],...
                'base_time','duration','spring_neap_days','end_time','items',...
                'num_proc_steps','num_timesteps','skipfac','ssample',...
                'start_time','suffix','tn','time','timestep','tn',...
                'type','u_num','v_num','xe','xn','ye','yn','ze','zn');
        catch %#ok<CTCH>
            dfsu.Close();
        end
    end

end

% Load the general info
fprintf('Load the saved basic data (%s)... ',[filename(1:end-5),'_raw.mat'])
load([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'_raw.mat'])
fprintf('done.\n')

if existed(1)==0
    % Need to extract the U and V components.
    u=nan(dfsu.NumberOfElements,num_proc_steps,'single');
    v=nan(dfsu.NumberOfElements,num_proc_steps,'single');

    % Make sure we omit the first three day's worth of values (on the safe
    % side).
    start_iter=345600/timestep;

    % Adjust number of iterations include the skip factor to extract only
    % every half hour's worth of data and adjusted for the first three
    % day's of values being omitted. We need to make sure this doesn't
    % extend beyond the end of the original file.
    if (num_proc_steps*skipfac)+start_iter>num_timesteps
        warning('Can''t omit first three days'' data because original output is too short. Check you didn''t account for this in the model setup.') %#ok<WNTAG>
        iter=1:skipfac:num_proc_steps*skipfac;
    else
        iter=start_iter:skipfac:(num_proc_steps*skipfac)+start_iter;
    end

    % How are we doing?
    h=waitbar(0,'Extracting vector data.','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0)

    for n=1:size(iter,2)
        if getappdata(h,'canceling')
            break
        end
        % Integrate the values of u and v over the data omitted due to
        % subsampling, if appropriate. I'm not sure this is right: so long
        % as I amend the value of timestep to account for the skipped data,
        % I don't need to sum the intermediate values. The residual
        % calculation takes the time between observations into account.
        if n==1
            getVal=iter(n)-1;
            % If magDir is one, we need to calculate U and V vectors from
            % magnitude and direction, otherwise we can just use them as
            % is.
            if magDir==0
                try
                    % We've got U and V already, so just carry on as normal.
                    u(:,n)=double(dfsu.ReadItemTimeStep(u_num,getVal).Data);
                    v(:,n)=double(dfsu.ReadItemTimeStep(v_num,getVal).Data);
                catch %#ok<CTCH>
                    break
                end
            elseif magDir==1
                try
                    % Calculate U and V from some simple trig.
                    % u_num is magnitude, v_num is direction.
                    u(:,n)=double(dfsu.ReadItemTimeStep(u_num,getVal).Data).*sind(double(dfsu.ReadItemTimeStep(v_num,getVal).Data));
                    v(:,n)=double(dfsu.ReadItemTimeStep(u_num,getVal).Data).*cosd(double(dfsu.ReadItemTimeStep(v_num,getVal).Data));
                catch %#ok<CTCH>
                    break
                end
            else
                dfsu.Close();
                error('Can''t find any direction or magnitude components. As such, can''t calculate U and V. Aborting.')
            end
        else
            getVal=iter(n)-1;
            uSum=nan(dfsu.NumberOfElements,skipfac-1);
            vSum=nan(dfsu.NumberOfElements,skipfac-1);
            % All these calls to dfsu makes this slow as hell again. So,
            % check if we need to sum data between subsamples (shouldn't if
            % I've made my output timestep sensible).
            if timestep<ssample
                warning('You should never see this message. If you do, the processing is going wrong.') %#ok<WNTAG>
                for s=1:skipfac-1
                    getVal=iter(n-s)-1;
                    if magDir==0
                        try
                            % We've got U and V already, so just carry on as normal.
                            uSum(:,s)=double(dfsu.ReadItemTimeStep(u_num,getVal).Data);
                            vSum(:,s)=double(dfsu.ReadItemTimeStep(v_num,getVal).Data);
                        catch %#ok<CTCH>
                            break
                        end
                    elseif magDir==1
                        try
                            % Calculate U and V from some simple trig.
                            % u_num is magnitude, v_num is direction.
                            uSum(:,n)=double(dfsu.ReadItemTimeStep(u_num,getVal).Data).*sind(double(dfsu.ReadItemTimeStep(v_num,getVal).Data));
                            vSum(:,n)=double(dfsu.ReadItemTimeStep(u_num,getVal).Data).*cosd(double(dfsu.ReadItemTimeStep(v_num,getVal).Data));
                        catch %#ok<CTCH>
                            break
                        end
                    else
                        dfsu.Close();
                        error('Can''t find any direction or magnitude components. As such, can''t calculate U and V. Aborting.')
                    end
                end
                u(:,n)=sum(uSum,2);
                v(:,n)=sum(vSum,2);
                clear uSum vSum
            else
                if magDir==0
                    try
                        % We've got U and V in the model output file, so just
                        % carry on as normal.
                        u(:,n)=double(dfsu.ReadItemTimeStep(u_num,getVal).Data);
                        v(:,n)=double(dfsu.ReadItemTimeStep(v_num,getVal).Data);
                    catch %#ok<CTCH>
                        break
                    end
                elseif magDir==1
                    try
                        % Calculate U and V from some simple trig.
                        % u_num is magnitude, v_num is direction.
                        u(:,n)=double(dfsu.ReadItemTimeStep(u_num,getVal).Data).*sind(double(dfsu.ReadItemTimeStep(v_num,getVal).Data));
                        v(:,n)=double(dfsu.ReadItemTimeStep(u_num,getVal).Data).*cosd(double(dfsu.ReadItemTimeStep(v_num,getVal).Data));
                    catch %#ok<CTCH>
                        break
                    end
                else
                    dfsu.Close();
                    error('Can''t find any direction or magnitude components. As such, can''t calculate U and V. Aborting.')
                end
            end
        end
        waitbar(n/size(iter,2),h);
    end
    delete(h)

    % Save the results in the appropriate .mat file
    fprintf('Saving flow vectors results to disk... ')
    save([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'.mat'],...
        'duration','num_proc_steps','u','v');
    fprintf('done.\n')
end

if existed(2)==0
    % We need to compute the residuals. We'll assume we've got the basic
    % and U and V vectors already extracted and loaded.
    x=nan(size(u),'single');
    y=nan(size(v),'single');

    h=waitbar(0,'Analysing residual flow data.','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0)
    for t=1:size(u,2)
        if getappdata(h,'canceling')
            break
        end
        % Sum the flow vectors in time.
        if t==1 % Units are m/s/m
            try
                x(:,t)=(u(:,t)*timestep);
                y(:,t)=(v(:,t)*timestep);
            catch %#ok<CTCH>
                delete(h)
            end
        else
            try
                x(:,t)=x(:,t-1)+(u(:,t)*timestep);
                y(:,t)=y(:,t-1)+(v(:,t)*timestep);
            catch %#ok<CTCH>
                delete(h)
            end
        end
        waitbar(t/size(u,2));
    end
    delete(h)

    % Save our hard work
    fprintf('Saving residual analysis results to disk... ')
    save([path,'/',filename(1:end-5),'_files/',filename(1:end-5),suffix,'res.mat'],...
        'x','y')
    fprintf('done.\n')
end

if sum(existed)==2
    % All the data have been loaded already, in which case, we're golden.
    fprintf('All ''%s'' data have been loaded successfully.\n',type)
end

