function [scan_id, tag] = ahe_importfile_ms(filename, dataLines)
%IMPORTFILE Import data from a text file
%  [SCAN_ID, TAG_1, TAG_2] = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the data as column
%  vectors.
%
%  [SCAN_ID, TAG_1, TAG_2] = IMPORTFILE(FILE, DATALINES) reads data for
%  the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  [scan_id, tag_1, tag_2] = importfile("E:\Research\COBRA\COBRA CTR Data Files\SPEC File Folder\ahe_scan_data_04032021.csv", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 13-May-2021 20:25:26

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "scan_id", "Var3", "tag_1", "tag_2", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11"];
opts.SelectedVariableNames = ["scan_id", "tag_1", "tag_2"];
opts.VariableTypes = ["string", "double", "string", "double", "double", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var3", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var3", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
scan_id = tbl.scan_id;
tag = strcat(num2str(tbl.tag_1),num2str(tbl.tag_2));

end