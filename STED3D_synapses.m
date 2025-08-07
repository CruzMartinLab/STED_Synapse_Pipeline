% This script runs Fiji to threshold and record 3D puncta data in STED
% images, then analyzes the separate presynaptic and postsynaptic puncta 
% data as well as the puncta overlap data. 
% Author: Grace Kirkpatrick
% Date: July 28, 2025
%
% Instructions:
% 1. Ensure that the 3D STED image to be analyzed has only 2 channels, one for the
% presynaptic marker and one for the postsynaptic marker. Extra channels can be
% rearranged and deleted in Fiji.
% 2. Run the script by clicking "Run" or typing "STED3D_synapses" into the 
% command window.
% 3. Select the 2-channel 3D STED image to be analyzed.


% --- Start of Code---

% Prompt user to select an image file
[filename, pathname] = uigetfile({'*.*', 'All Files (*.*)'}, 'Select an image file');

if isequal(filename,0)
    disp('User canceled file selection');
    return;
end

imagePath = fullfile(pathname, filename);
imagePath = strrep(imagePath, '\', '/');

% Extract output directory and basename (without extension)
[outputDir, name, ext] = fileparts(imagePath);

% Remove additional extension .ome if present in original image name
% Any image returned by Huygens deconvolution should end in .ome so no
% other extensions need to be checked
if endsWith(name, '.ome', 'IgnoreCase', true)
    basename = extractBefore(name, '.ome');
else
    basename = name;
end

% Define results folder where Fiji macro will save outputs
resultsFolder = fullfile(outputDir, [basename '_Results']);
resultsFolder = strrep(resultsFolder, '\', '/');

% Make sure the results folder exists BEFORE running Fiji
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

% Assign channels 
ch1 = 'PSD95';
ch2 = 'RIM1';

% Combine args for command
args = sprintf('%s;%s;%s;%s', imagePath, resultsFolder, ch1, ch2);


%logFile = 'C:/Users/kirkpagr/Desktop/20250616/fiji_macro_log.txt';

% Paths to Fiji and macro
%fijiPath = 'C:/Users/kirkpagr/Documents/Fiji/fiji-windows-x64.exe';
fijiDir = 'C:/Users/kirkpagr/Downloads/Fiji';
%fijiExe = 'fiji-windows-x64.exe';
fijiExe = 'ImageJ-win64.exe';
macroPath = 'C:/Users/kirkpagr/Documents/MATLAB/STED3D_synapses_macro.ijm'; % name and location of Fiji macro

% Run Fiji macro
%cmd = sprintf('start /wait "" "%s" --headless --console -macro "%s" "%s" > "%s" 2>&1', fijiPath, macroPath, args, logFile);
% cmd = sprintf('start /wait "" "%s" --headless --console -macro "%s" "%s"', fijiPath, macroPath, args);
%cmd = sprintf('"%s" --headless --console -macro "%s" "%s"', fijiPath, macroPath, args);
% cmd = sprintf('cmd /C "%s --headless --console -macro \\"%s\\" \\"%s\\""', fijiPath, macroPath, args);

%cmd = sprintf('cmd /C "cd /d %s && %s --headless --console -macro \\"%s\\" \\"%s\\""', fijiDir, fijiExe, macroPath, args);
cmd = sprintf('cmd /S /C ""%s\\%s" --headless -macro "%s" "%s"""', fijiDir, fijiExe, macroPath, args);
disp(args)
disp(['Fiji path: ', fijiPath])
disp(['Macro path: ', macroPath])
disp(['Running command: ', cmd]);
[status, result] = system(cmd);
disp(result);

% Check that Fiji ran succcessfully
if status ~= 0
    error('Fiji macro failed to run. Exit code: %d', status);
end

% Build paths to macro output files
PSDresultsCSV = fullfile(resultsFolder, 'PSD95_3DResults.csv');
RIMresultsCSV = fullfile(resultsFolder, 'RIM1_3DResults.csv');
PSDobjectmapTIF = fullfile(resultsFolder, 'PSD95_ObjectMap.tif');
RIMobjectmapTIF = fullfile(resultsFolder, 'RIM1_ObjectMap.tif');

% Check that files exist before loading
if ~isfile(PSDobjectmapTIF)
    error('PSD95 object map not found: %s', PSDobjectmapTIF);
end
if ~isfile(RIMobjectmapTIF)
    error('RIM1 object map not found: %s', RIMobjectmapTIF);
end

% Load labeled PSD95 object map as matrix while retaining object IDs
PSDinfo = imfinfo(PSDobjectmapTIF);
PSDnumSlices = numel(PSDinfo);
PSDheight = PSDinfo(1).Height;
PSDwidth = PSDinfo(1).Width;

PSDobjectmap = zeros(PSDheight, PSDwidth, PSDnumSlices, 'uint16');

for k = 1:PSDnumSlices
    PSDobjectmap(:,:,k) = imread(PSDobjectmapTIF, k);
end

% Load labeled RIM1 object map as matrix while retaining object IDs
RIMinfo = imfinfo(RIMobjectmapTIF);
RIMnumSlices = numel(RIMinfo);
RIMheight = RIMinfo(1).Height;
RIMwidth = RIMinfo(1).Width;

RIMobjectmap = zeros(RIMheight, RIMwidth, RIMnumSlices, 'uint16');

for k = 1:RIMnumSlices
    RIMobjectmap(:,:,k) = imread(RIMobjectmapTIF, k);
end

% Load results table for PSD95
if isfile(PSDresultsCSV)
    PSDresults = readtable(PSDresultsCSV);
    disp('Loaded PSD95 results table:');
    disp(PSDresults(1:min(5,height(PSDresults)), :));

%     % Add column with PSD95 ID's
%     psdIDs = 1:height(PSDresults);
%     PSDresults = addvars(PSDresults, psdIDs, 'NewVariableNames', 'PSD95_ID');
else
    warning('Results table not found: %s', PSDresultsCSV);
end

% Load results table for RIM1
if isfile(RIMresultsCSV)
    RIMresults = readtable(RIMresultsCSV);
    disp('Loaded RIM1 results table:');
    disp(RIMresults(1:min(5,height(RIMresults)), :));

    % % Add column with RIM1 ID's
    % rimIDs = 1:height(RIMresults);
    % RIMresults = addvars(RIMresults, rimIDs, 'NewVariableNames', 'RIM1_ID');
else
    warning('Results table not found: %s', RIMresultsCSV);
end

% Preallocate table to store overlap data
overTypes = ["string", "string", "string", "double", "double", "double"];
overCols = ["Synapse ID", "PSD95 ID", "RIM1 ID", "PSD95 volume", "RIM1 volume", "Overlap volume"];
overTable = table('Size', [0 6], 'VariableTypes', overTypes, 'VariableNames', overCols);

% Initialize overlap matrix 
[rows, cols, stacks] = size(PSDobjectmap);
overlapMatrix = zeros(rows, cols, stacks, 'uint16');

% Track which PSD95-RIM1 pairs we've already counted 
seenPairs = containers.Map;


% Add to overlap matrix and overlap table 
for i = 1:rows
    for j = 1:cols
        for k = 1:stacks
            % Get labels
            psdID = PSDobjectmap(i,j,k);
            rimID = RIMobjectmap(i,j,k);

            % Skip background and number labels (value = 65535)
            if psdID == 0 || rimID == 0 || psdID == 65535 || rimID == 65535
                continue
            end

            % Construct a unique key for this pair
            pairKey = sprintf('%d_%d', psdID, rimID);

            % Only process if this pair hasn't been added yet
            if ~isKey(seenPairs, pairKey)
                % Mark this pair as seen
                seenPairs(pairKey) = true;
                synapseID = length(seenPairs);
       
                % Compute overlap mask and volume
                overlapMask = (PSDobjectmap == psdID) & (RIMobjectmap == rimID);
                overlapVolume = sum(overlapMask(:)); 

                % Find corresponding PSD95 and RIM1 volumes
                psdVol = sum(PSDobjectmap == psdID);
                rimVol = sum(RIMobjectmap == rimID);

                % Assign unique value to this overlap in overlap matrix
                overlapMatrix(overlapMask) = synapseID;

                % Store overlap data in the table
                newRow = {synapseID, psdID, rimID, psdVol, rimVol, overlapVolume};
                overTable = [overTable; newRow]; 
            end
        end
    end
end
