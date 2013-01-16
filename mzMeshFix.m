function mzMeshFix(mesh,attributes)
% MZMESHFIX Use on a MESH output from mzMeshMerge to eliminate arc
% attributes of value(s) specified in ATTRIBUTES.
%
%   MZMESHFIX(MESH,ATTRIBUTES) where MESH is a given mesh file and
%   ATTRIBUTES is a single or multiple attribute values. Output mesh is
%   input base name (i.e. no extension) appended with _fixed.
%
%   ATTRIBUTES can be an array of values or a string number (i.e. '20'). If
%   it is a string, the str2double of that is used as a minimum value. Any
%   arc attribute greater than that will be set to zero. This is useful if
%   you have set all your internal arc attributes to very high values
%   relative to the 'real' arc attributes. It significantly decreases the
%   run time for MZMESHFIX.
%
% Pierre Cazenave 17/01/2012
%   v1.0 - initial version.

if nargin~=2
    error('Incorrect number of arguments specified. See ''help MZMESHFIX'' for more information.')
end

fid = fopen(mesh);
if fid == -1
    error('Could not open input file %s.',mesh)
end

% Count which line we're on.
counter = 1;
% To count where in the output array we've got to.
i = 1;

% Create output file
fidOut = fopen([mesh(1:end-5),'_fixed.mesh'],'w');

tline = fgets(fid);
while ischar(tline)
    if counter == 1
        % Print the header
        fprintf(fidOut,'%s',tline);
        counter=counter+1;
    elseif counter > 1
        test = textscan(tline,'%s');
        if numel(test{1})~=5
            % We're only fiddling around with the arc attributes, so we're
            % only interested in lines with five values. If there aren't
            % five, just print the line to the file as is.
            fprintf(fidOut,'%s',tline);
        else
            if ischar(attributes);
                % We have a minimum, so don't loop through a series
                % of values, only keep those below that minimum.
                checkAtt=str2double(attributes);
                if isnan(checkAtt)
                    error('Please specify a numeric value for ''ATTRIBUTES''')
                end
                if numel(test{1}) == 5 && str2double(char(test{1}(end))) > checkAtt
                    % We can't preallocate because we don't know how many
                    % boundary points we have. Also, this is fugly...
                    oline(i,1) = str2double(char(test{1}(1)));
                    oline(i,2) = str2double(char(test{1}(2)));
                    oline(i,3) = str2double(char(test{1}(3)));
                    oline(i,4) = str2double(char(test{1}(4)));
                    oline(i,5) = 0;
                    fprintf(fidOut,'%d %f %f %f %d\n',oline(i,:));
                    i = i+1;
                elseif numel(test{1}) == 5 && str2double(char(test{1}(end))) <= checkAtt
                    fprintf(fidOut,'%s',tline);
                end
            else
                % Here, we need to check how many arc attributes we're
                % looking for.
                if numel(attributes) == 1 % easy, only one value to check.
                    if numel(test{1}) == 5 && str2double(char(test{1}(end))) == attributes
                        % We can't preallocate because we don't know how
                        % many boundary points we have. Also, this is
                        % fugly...
                        oline(i,1) = str2double(char(test{1}(1)));
                        oline(i,2) = str2double(char(test{1}(2)));
                        oline(i,3) = str2double(char(test{1}(3)));
                        oline(i,4) = str2double(char(test{1}(4)));
                        oline(i,5) = 0;
                        fprintf(fidOut,'%d %f %f %f %d\n',oline(i,:));
                        i = i+1;
                    elseif numel(test{1}) == 5 && str2double(char(test{1}(end))) ~= attributes
                        fprintf(fidOut,'%s',tline);
                    end
                else % More complicated, but nothing a for-loop can't solve.
                    if numel(test{1}) == 5
                        % We can't preallocate because we don't know how
                        % many boundary points we have. Also, this is
                        % fugly...
                        oline(i,1) = str2double(char(test{1}(1)));
                        oline(i,2) = str2double(char(test{1}(2)));
                        oline(i,3) = str2double(char(test{1}(3)));
                        oline(i,4) = str2double(char(test{1}(4)));
                        oline(i,5) = str2double(char(test{1}(end)));
                        for checkBoundary=1:numel(attributes)
                            % Check all the attributes and if it matches
                            % any of them, overwrite the existing output
                            % with zero.
                            if str2double(char(test{1}(end))) == attributes(checkBoundary)
                                oline(i,5) = 0;
                            end
                        end
                        fprintf(fidOut,'%d %f %f %f %d\n',oline(i,:));
                        i = i+1;
                    end
                end
            end
        end
    end
    tline = fgets(fid);
    counter = counter+1;
end

fclose(fid);
fclose(fidOut);