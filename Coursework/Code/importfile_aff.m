function output = importfile_aff(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   OUTPUT = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   OUTPUT = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   output = importfile('output.out', [8,10], [8,10]);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/04/11 17:44:03

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = [8,10];
    endRow = [8,10];
end

%% Format for each line of text:
%   column1: categorical (%C)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%C%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
output = table(dataArray{1:end-1}, 'VariableNames', {'Resourcesallocatedcpupercent0'});

