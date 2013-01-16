function [X,Y,z,zFixed]=dfsuRep(filename,repData,wesn)
% DFSUREP(FILENAME,REPDATA) takes a given dfsu file path (FILENAME) and
% searches for all instances of the first attribute which aren't contained
% in the array REPDATA.
%
% All values which aren't those in REPDATA are replaced by the closest value
% in REPDATA, effectively eliminating the results of the interpolation
% resulting in much sharper boundaries between adjacent groups of values.
% This is particularly important when considering grain size data, where
% boundaries between different sediment types may be very sharp.
%
% DFSUREP(FILENAME,REPDATA,WESN) works as above, but limits replacements to
% the subset described by the array WESN. WESN contains coordinates of a
% bounding box of interest:
%
% WESN=[0.5 0.7 50.66 50.8];
%
% This WESN limits changes to a region in the English Channel (0.5E to 0.7E
% and 50.66N to 50.8N). Specifying a subdomain increases the time to run
% the command since each value which does not match a value in REPDATA must
% be checked to see if it lies within the subdomain coordinates in WESN.
%
% [X,Y,Z,ZFIXED]=DFSUREP(...) will output the X and Y mesh coordinates as
% well as the original Z values and the fixed ones for use in plotting. If
% these outputs are omitted, the function will merely create the new dfsu
% file.
%
% Due to the way the element values are represented by the face as the
% average of the three adjacent nodal values, it is impossible to get a
% completely clean break, but this approach should improve the situation.
%
% Pierre Cazenave 2011/09/07 v1.0
%                 2011/09/12 v1.1 Added [X,Y,Z,ZFIXED] to outputs for use
%                 in plotting from within a script.

% Right then, onwards...

if nargin==2
    wesn=nan;
elseif nargin~=3
    error('Please supply the correct number of arguments.')
end

%%% TESTING INPUTS
% The functionality for working within a subset severely increases run time
% because we have to search through each result to make sure we only change
% those values which fall within that subset.
% wesn=[-4.5,-2.5,50.8,51.9]; % This is Culver sands
% wesn=[0.5 0.7 50.66 50.8]; % This is Hastings Shingle Bank
% wesn=[-2.5 2.5 48 52]; % Central English Channel

% filename='Z:\modelling\data\sediment_distribution\round_10_culver\bank_transport\csm_culver_v4_seds.dfsu';
% Numbers are a combination of those in convert_bgs_seds.m,
% convert_planpain_seds.m and extracted_grain_size_information.xlsx for the
% Culver Sands.
% Shelf-wide grain sizes
% repData=[0.03,0.04194,0.04567,0.05,0.07108,0.0737,0.08,0.09649,0.12452,...
%     0.12825,0.15628,0.2,0.28,0.2871,0.31886,0.33,0.42463,0.66233,...
%     0.70883,0.74491,1.14,1.25215,3.39517,6.41,8.93,19.375,20,76];
% Culver sands grain sizes (BGS and backscatter derived averages)
% repData=[0.013910,0.045670,0.096490,0.124520,0.128250,0.156280,0.28,...
%     0.287100,0.318860,0.662330,0.708830,0.744910,1.252150,3.395170,...
%     19.375,20,76];

%%% END OF TESTING INPUTS

% Check our coordinates make sense.
if max(isnan(wesn))~=1
    if wesn(1)>wesn(2) || wesn(3)>wesn(4)
        error('Supplied bounding box coordinates do not make sense. Make sure order is correct (WESN)')
    end
end

debugMe=0; % Debugging: 1=on, 0=off.

% Find the midpoints between each group in REPDATA so we can find the
% elements within those ranges and set all those values to the value
% specified in REPDATA, making for a much more discrete dfsu file. Remove
% duplicates in REPDATA too.
valuesRep=unique(repData);
midValues=nan(size(valuesRep,2)-1,1);
for i=1:size(valuesRep,2)
    if i<size(valuesRep,2)
        midValues(i)=valuesRep(i)+((valuesRep(i+1)-valuesRep(i))/2);
    end
end

% Set up the new MIKE interface to dfs files.
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;

% Read in the dfsu file.
dfsu=DfsFileFactory.DfsuFileOpen(filename);
% Have to store it as a single otherwise we get very small differences in
% the actual values from the dfsu file.
z=single(dfsu.ReadItemTimeStep(1,0).Data);

% Close the dfsu file handle.
dfsu.Close()

% Create element table in Matlab format.
tn=mzNetFromElmtArray(dfsu.ElementTable);
% Node coordinates.
xn=double(dfsu.X);
yn=double(dfsu.Y);
% zn=double(dfsu.Z); (unused here)
% Element centres (unused here).
% [xe,ye,ze]=mzCalcElmtCenterCoords(tn,xn,yn,zn);

% General position data for the patch command.
X(:,1)=xn(tn(:,1));
X(:,2)=xn(tn(:,2));
X(:,3)=xn(tn(:,3));
Y(:,1)=yn(tn(:,1));
Y(:,2)=yn(tn(:,2));
Y(:,3)=yn(tn(:,3));

% We need to find all values which aren't one of the values in the REPDATA
% array.
posToChange=cell(1,size(repData,2));
zFixed=z;

for val=1:size(repData,2)
    % For the first and last data points in REPDATA, we don't have a lower
    % or upper bound, so use less than/greater than instead.
    if val==1
        if debugMe
            fprintf('First: z (%.5f) < %.5f\n',repData(val),midValues(val)) %#ok<UNRCH>
        end
        posToChange{val}=find(z<midValues(val));
    elseif val>1 && val<size(repData,2)
        if debugMe
            fprintf('Middle: %.5f < z (%.5f) < %.5f\n',midValues(val-1),repData(val),midValues(val)) %#ok<UNRCH>
        end
        posToChange{val}=find(z>midValues(val-1) & z<midValues(val));
    else
        if debugMe
            fprintf('End: z (%.5f) > %.5f\n',repData(val),midValues(val-1)) %#ok<UNRCH>
        end
        posToChange{val}=find(z>midValues(val-1));
    end

    % Using those found indices, replace all those values with the
    % appropriate value from REPDATA.
    if debugMe
        if val==1 %#ok<UNRCH>
            fprintf('Replacing %i values: z (%.5f) < %.5f\n',numel(posToChange{val}(:)),repData(val),midValues(val))
        elseif val>1 && val<size(repData,2)
            fprintf('Replacing %i values: %.5f < z (%.5f) < %.5f\n',numel(posToChange{val}(:)),midValues(val-1),repData(val),midValues(val))
        else
            fprintf('Replacing %i values: z (%.5f) > %.5f\n',numel(posToChange{val}(:)),repData(val),midValues(val-1))
        end
    end

    if max(isnan(wesn))==1
        % We're not constrained to a particular location, so we can just to
        % the lot all at once.
        zFixed(posToChange{val})=repData(val);
    else
        % We have to check every result against the corresponding latitude
        % and longitude and replace only if the data fall in that area.
        for che=1:size(posToChange{val},2)
            posX=X(posToChange{val}(che));
            posY=Y(posToChange{val}(che));
            if (posX>wesn(1) && posX<wesn(2) && posY>wesn(3) && posY<wesn(4))
                if debugMe
                    fprintf('Replaced %.5f at %.5f/%.5f with %.5f.\n',zFixed(posToChange{val}(che)),posX,posY,repData(val)) %#ok<UNRCH>
                end
                zFixed(posToChange{val}(che))=repData(val);
            end
        end
    end
end

% Some debug info.
if debugMe
    fprintf('Replaced %i of %i total values.\n',cellSum,size(z,2)) %#ok<UNRCH>
end

% We shouldn't have replaced more data points than existed in the input...
cellSum=0; for i=1:size(posToChange,2); cellSum=cellSum+numel(posToChange{i}(:)); end
if cellSum>size(z,2)
    warning('Replaced more data than existed in the input file... ') %#ok<WNTAG>
end

% Create a new output file and write the modified values.
if exist([filename(1:end-5),'_fixed.dfsu'],'file')~=2
    % Copy the input file and use as the basis for the modified output.
    copyfile(filename,[filename(1:end-5),'_fixed.dfsu'])
    % Open the newly created file.
    dfsuOut=DfsFileFactory.DfsuFileOpenEdit([filename(1:end-5),'_fixed.dfsu']);

    % Read first time step from file.
    itemData=dfsuOut.ReadItemTimeStep(1,0);
    % Write to memory.
    dfsuOut.WriteItemTimeStep(1,0,itemData.Time,NET.convertArray(single(zFixed(:))));
    % Save and close the file.
    dfsuOut.Close()
else
    error('Will not overwrite existing file.')
end