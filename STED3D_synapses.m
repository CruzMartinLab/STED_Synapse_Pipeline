% This script runs Fiji to threshold and record 3D puncta data in STED
% images, then analyzes the separate presynaptic and postsynaptic puncta 
% data as well as the puncta overlap data. 
% Author: Grace Kirkpatrick
% Date: July 28, 2025
%
% Instructions:
% 1. Choose the correct directory that includes all MATLAB and Fiji files
% (including the macro and all Miji files). 
% 2. Ensure that the 3D STED image to be analyzed has only 2 channels, one for the
% presynaptic marker and one for the postsynaptic marker. Extra channels can be
% rearranged and deleted in Fiji.
% 3. Run the script by clicking "Run" or typing "STED3D_synapses" into the 
% command window.
% 4. Select the 2-channel 3D STED image to be analyzed.


%% Image Analysis

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

% Use existing results folder or create a new one 
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

% Assign channels 
ch1 = 'PSD95';
ch2 = 'RIM1';

% Combine args for command
args = sprintf('%s;%s;%s;%s', imagePath, resultsFolder, ch1, ch2);

% Path to macro
macroPath = 'C:/Users/kirkpagr/Documents/MATLAB/STED3D_synapses_macro.ijm'; % name and location of Fiji macro

% Run Fiji macro
Miji(false); % false = do NOT show GUI
ij.IJ.runMacroFile(macroPath, args);

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
else
    warning('Results table not found: %s', PSDresultsCSV);
end

% Load results table for RIM1
if isfile(RIMresultsCSV)
    RIMresults = readtable(RIMresultsCSV);
    disp('Loaded RIM1 results table:');
    disp(RIMresults(1:min(5,height(RIMresults)), :));
else
    warning('Results table not found: %s', RIMresultsCSV);
end

%% Compute Overlaps

voxelSizeX = 0.05; % um
voxelSizeY = 0.05;
voxelSizeZ = 0.05;
voxelVolume = voxelSizeX * voxelSizeY * voxelSizeZ; % um^3

% Initialize overlap matrix 
[rows, cols, stacks] = size(PSDobjectmap);
overlapMatrix = zeros(rows, cols, stacks, 'uint16');

% Store all PSD and RIM labels (IDs)
psdLabels = unique(PSDobjectmap(:));
psdLabels = psdLabels(psdLabels > 0  & psdLabels < 65535);

rimLabels = unique(RIMobjectmap(:));
rimLabels = rimLabels(rimLabels > 0 & rimLabels < 65535);

% Preallocate table to store overlap data
% overTypes = ["string", "string", "string", "double", "double", "double"];
overCols = ["Synapse ID", "PSD95 ID", "RIM1 ID", "PSD95 volume", "RIM1 volume", "Overlap volume"];
maxSynapses = length(psdLabels) * 10; % overestimate max number of possible overlaps
overTableData = cell(maxSynapses, 6);

% Initialize synapse ID
synapseID = 0;

% Loop through all PSD objects
for i = 1:length(psdLabels)
    psdID = psdLabels(i);
    psdMask = (PSDobjectmap == psdID);

    % Find which RIM labels overlap with this PSD
    overlappingRIM = unique(RIMobjectmap(psdMask));
    overlappingRIM = overlappingRIM(overlappingRIM > 0 & overlappingRIM < 65535);

    % Loop through all RIM that overlap with this PSD
    for j = 1:length(overlappingRIM)
        rimID = overlappingRIM(j);

        % Compute overlap mask
        overlapMask = psdMask & (RIMobjectmap == rimID);
        overlapVolume = sum(overlapMask(:));

        % Get full volumes of current PSD and RIM objects
        psdVol = sum(psdMask(:));
        rimVol = sum(RIMobjectmap(:) == rimID);

        % Assign unique synapse ID
        synapseID = synapseID + 1;
        overlapMatrix(overlapMask) = synapseID;

        % Save overlap data
        overTableData(synapseID, :) = {synapseID, psdID, rimID, psdVol, rimVol, overlapVolume};
    end
end

% Make overlap matrix binary and convert to image, then save in results folder
binOM = overlapMatrix > 0;
saveastiff(uint8(binOM * 255), fullfile(resultsFolder, 'overlap.tif'));
saveastiff(overlapMatrix, fullfile(resultsFolder, 'overlap_labeled.tif'));

% Trim overlap table to remove empty rows
overTableData = overTableData(1:synapseID, :);

% Convert overlap table data to a table format
overTable = cell2table(overTableData, 'VariableNames', overCols);

%% Data Analysis

% Save overlap table CSV to results folder
writetable(overTable, fullfile(resultsFolder, 'OverlapResults.csv'));

% Count synapses per PSD95
[counts, psd_ids] = groupcounts(overTable.("PSD95 ID"));
countTable = table(psd_ids, counts, 'VariableNames', {'PSD95_ID', 'NumSynapses'});

% Sum overlap volume per PSD95 using splitapply (more robust)
[totalVols, uniqueIDs] = findgroups(overTable.("PSD95 ID"));
summedOverlap = splitapply(@sum, overTable.("Overlap volume"), totalVols);
overlapTable = table(uniqueIDs, summedOverlap, 'VariableNames', {'PSD95_ID', 'TotalOverlapVol'});

% PSD95 volume (take first instance per PSD95_ID)
[uniqueIDs, ia] = unique(overTable.("PSD95 ID"));
psdVols = overTable.("PSD95 volume")(ia);
volumeTable = table(uniqueIDs, psdVols, 'VariableNames', {'PSD95_ID', 'PSD95_Vol'});

% Merge all PSD summary tables
PsummaryTable = outerjoin(countTable, overlapTable, 'Keys', 'PSD95_ID', 'MergeKeys', true);
PsummaryTable = outerjoin(PsummaryTable, volumeTable, 'Keys', 'PSD95_ID', 'MergeKeys', true);

PsummaryTable.TotalOverlapVol_um3 = PsummaryTable.TotalOverlapVol * voxelVolume;
PsummaryTable.PSD95_Vol_um3 = PsummaryTable.PSD95_Vol * voxelVolume;
PsummaryTable.OverlapPercent = 100 * (PsummaryTable.TotalOverlapVol ./ PsummaryTable.PSD95_Vol);

% Save PSD95 summary table CSV to results folder
writetable(PsummaryTable, fullfile(resultsFolder, 'PSD95_Summary.csv'));

% Count synapses per RIM1
[counts, psd_ids] = groupcounts(overTable.("RIM1 ID"));
countTable = table(psd_ids, counts, 'VariableNames', {'RIM1_ID', 'NumSynapses'});

% Sum overlap volume per RIM
[totalVols, uniqueIDs] = findgroups(overTable.("RIM1 ID"));
summedOverlap = splitapply(@sum, overTable.("Overlap volume"), totalVols);
overlapTable = table(uniqueIDs, summedOverlap, 'VariableNames', {'RIM1_ID', 'TotalOverlapVol'});

% RIM volume (take first instance per RIM1_ID)
[uniqueIDs, ia] = unique(overTable.("RIM1 ID"));
psdVols = overTable.("RIM1 volume")(ia);
volumeTable = table(uniqueIDs, psdVols, 'VariableNames', {'RIM1_ID', 'RIM1_Vol'});

% Merge all RIM summary tables
RsummaryTable = outerjoin(countTable, overlapTable, 'Keys', 'RIM1_ID', 'MergeKeys', true);
RsummaryTable = outerjoin(RsummaryTable, volumeTable, 'Keys', 'RIM1_ID', 'MergeKeys', true);

RsummaryTable.TotalOverlapVol_um3 = RsummaryTable.TotalOverlapVol * voxelVolume;
RsummaryTable.RIM1_Vol_um3 = RsummaryTable.RIM1_Vol * voxelVolume;
RsummaryTable.OverlapPercent = 100 * (RsummaryTable.TotalOverlapVol ./ RsummaryTable.RIM1_Vol);

% Save RIM1 summary table CSV to results folder
writetable(RsummaryTable, fullfile(resultsFolder, 'RIM1_Summary.csv'));

%% Plot Data

% PSD95 volume (in um^3) vs. number of synapses it makes
figure; 
scatter(PsummaryTable.PSD95_Vol_um3, PsummaryTable.NumSynapses, 50, 'filled');
xlabel('PSD95 Volume (µm^3)');
ylabel('Number of Synapses');
title('PSD95 Volume vs Number of Synapses');
grid on;
% Add linear fit line
hold on;
p = polyfit(PsummaryTable.PSD95_Vol_um3, PsummaryTable.NumSynapses, 1);
xFit = linspace(min(PsummaryTable.PSD95_Vol_um3), max(PsummaryTable.PSD95_Vol_um3), 100);
yFit = polyval(p, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2);
hold off;
% Compute correlation
[r1, p1] = corr(PsummaryTable.PSD95_Vol_um3, PsummaryTable.NumSynapses);
r_2 = r1*r1;
txt = sprintf('r = %.2f, r^2 = %.2f, p = %.3g', r1, r_2, p1);
% Place text in upper left corner of plot
xLimits = xlim;
yLimits = ylim;
xText = xLimits(1) + 0.05 * range(xLimits);
yText = yLimits(2) - 0.1 * range(yLimits);
text(xText, yText, txt, 'FontSize', 12, 'FontWeight', 'bold');

% Save the figure
filename = fullfile(resultsFolder, 'PSD_vs_Synapses.png');
saveas(gcf, filename); 

% % PSD95 volume (in um^3) vs. percentage of volume with overlap
% figure;
% scatter(PsummaryTable.PSD95_Vol_um3, PsummaryTable.OverlapPercent, 50, 'filled');
% xlabel('PSD95 Volume (µm³)');
% ylabel('Overlap Percentage (%)');
% title('PSD95 Volume vs Overlap Percentage');
% grid on;
% % Add linear fit line
% hold on;
% p = polyfit(PsummaryTable.PSD95_Vol_um3, PsummaryTable.OverlapPercent, 1);
% xFit = linspace(min(PsummaryTable.PSD95_Vol_um3), max(PsummaryTable.PSD95_Vol_um3), 100);
% yFit = polyval(p, xFit);
% plot(xFit, yFit, 'r-', 'LineWidth', 2);
% hold off;

% RIM1 volume (in um^3) vs. number of synapses it makes
figure; 
scatter(RsummaryTable.RIM1_Vol_um3, RsummaryTable.NumSynapses, 50, 'filled');
xlabel('RIM1 Volume (µm^3)');
ylabel('Number of Synapses');
title('RIM1 Volume vs Number of Synapses');
grid on;
% Add linear fit line
hold on;
p = polyfit(RsummaryTable.RIM1_Vol_um3, RsummaryTable.NumSynapses, 1);
xFit = linspace(min(RsummaryTable.RIM1_Vol_um3), max(RsummaryTable.RIM1_Vol_um3), 100);
yFit = polyval(p, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2);
hold off;
% Compute correlation
[r1, p1] = corr(RsummaryTable.RIM1_Vol_um3, RsummaryTable.NumSynapses);
r_2 = r1*r1;
txt = sprintf('r = %.2f, r^2 = %.2f, p = %.3g', r1, r_2, p1);
% Place text in upper left corner of plot
xLimits = xlim;
yLimits = ylim;
xText = xLimits(1) + 0.05 * range(xLimits);
yText = yLimits(2) - 0.1 * range(yLimits);
text(xText, yText, txt, 'FontSize', 12, 'FontWeight', 'bold');

% Save the figure
filename = fullfile(resultsFolder, 'RIM_vs_Synapses.png');
saveas(gcf, filename); 

% % RIM1 volume (in um^3) vs. percentage of volume with overlap
% figure;
% scatter(RsummaryTable.RIM1_Vol_um3, RsummaryTable.OverlapPercent, 50, 'filled');
% xlabel('RIM1 Volume (µm³)');
% ylabel('Overlap Percentage (%)');
% title('RIM1 Volume vs Overlap Percentage');
% grid on;
% % Add linear fit line
% hold on;
% p = polyfit(RsummaryTable.RIM1_Vol_um3, RsummaryTable.OverlapPercent, 1);
% xFit = linspace(min(RsummaryTable.RIM1_Vol_um3), max(RsummaryTable.RIM1_Vol_um3), 100);
% yFit = polyval(p, xFit);
% plot(xFit, yFit, 'r-', 'LineWidth', 2);
% hold off;

% PSD95 volume vs. RIM1 volume in each synapse
figure; 
scatter(overTable.("PSD95 volume"), overTable.("RIM1 volume"), 50, 'filled');
xlabel('PSD95 Volume (voxels)');
ylabel('RIM1 Volume (voxels)');
title('PSD95 Volume vs. RIM1 Volume per Synapse');
grid on;
% Add linear fit line 
hold on;
p = polyfit(overTable.("PSD95 volume"), overTable.("RIM1 volume"), 1);
xFit = linspace(min(overTable.("PSD95 volume")), max(overTable.("PSD95 volume")), 100);
yFit = polyval(p, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2);
hold off;
% Compute correlation
[r1, p1] = corr(overTable.("PSD95 volume"), overTable.("RIM1 volume"));
r_2 = r1*r1;
txt = sprintf('r = %.2f, r^2 = %.2f, p = %.3g', r1, r_2, p1);
% Place text in upper left corner of plot
xLimits = xlim;
yLimits = ylim;
xText = xLimits(1) + 0.05 * range(xLimits);
yText = yLimits(2) - 0.1 * range(yLimits);
text(xText, yText, txt, 'FontSize', 12, 'FontWeight', 'bold');

% Save the figure
filename = fullfile(resultsFolder, 'PSD_vs_RIM.png');
saveas(gcf, filename); 