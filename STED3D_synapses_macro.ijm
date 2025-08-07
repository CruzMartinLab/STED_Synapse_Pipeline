// STED3D_synapses_macro.ijm
// Macro for thresholding and analyzing multi-stack presynaptic and postsynaptic STED channels
// To be run through STED3D_synapses.m using Miji MATLAB integration

// Get the image path and channel assignments from MATLAB
args = getArgument();

// To test Fiji macro indpendently of MATLAB, manually input argument line (imagePath, resultsFolder, ch1, ch2)
// args = 

// This should never be a problem as long as MATLAB file selection code is not edited
if (args == "") {
    exit("No input file path provided.");
}

// Reset Fiji windows
run("Close All");
setBatchMode(true);


// Parse arguments
tokens = split(args,";");
if (lengthOf(tokens) < 4){
	exit("Missing arguments. Got: " + args);
}
oldImagePath = tokens[0];
resultsFolder = tokens[1];
ch1 = tokens[2];
ch2 = tokens[3];
ch2 = replace(ch2, "\"", "");

// Logging
logPath = tokens[1] + File.separator + "macro_log.txt";
File.append("Starting macro...\n", logPath);
File.append("Args: " + args + "\n", logPath);

// Extract base filename (no extension)
parts = split(oldImagePath, "/");
filename = parts[lengthOf(parts) - 1];
dotIndex = lastIndexOf(filename, ".");
if (dotIndex > -1){
	basename = substring(filename, 0, dotIndex);
}
else{
	basename = filename;
}

// Create new filename and copy image to results folder
shortName = "STEDinput.tif";
imagePath = resultsFolder + File.separator + shortName;
File.copy(oldImagePath, imagePath);
print("Copied original image to: " + imagePath);

// Import copied image using Bio-Formats
run("Bio-Formats Importer", "open=[" + imagePath + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT quiet");

// Update basename for all downstream naming
basename = "STEDinput";
rename(basename);

// Process PSD95 channel (Ch1)
selectWindow(basename);
run("Duplicate...", "duplicate channels=1 title=" + ch1);
threshold_and_3DOC(ch1);

// Process RIM1 channel (Ch2)
selectWindow(basename);
run("Duplicate...", "duplicate channels=2 title=" + ch2);
threshold_and_3DOC(ch2);

selectWindow(basename);
close(); // Clean up original hyperstack

run("Close All");
setBatchMode(false);
print("Macro finished");
exit();

// Function for thresholding and running 3D Objects Counter on one channel
function threshold_and_3DOC(label) {
    
    if (!isOpen(label)) {
        print("Channel window not found: " + label);
        return;
    }

    selectWindow(label);
    
    // Same threshold applied to all channels
    run("8-bit");
	run("Auto Local Threshold", "method=Bernsen radius=15 parameter_1=0 parameter_2=0 white stack");
	
	// Save thresholded channel
    saveAs("Tiff", resultsFolder + File.separator + label + "_thresholded.tif");
    
    // Run 3D Objects Counter on thresholded channel
    run("3D Objects Counter", "threshold = 128 min=40 max=99999999 objects statistics summary");
    
    // Save labeled objects map
	objectMapTitle = "Objects map of " + label + "_thresholded.tif";
    if (isOpen(objectMapTitle)) {
        selectWindow(objectMapTitle);
        saveAs("Tiff", resultsFolder + File.separator + label + "_ObjectMap.tif");
        close();
    }
    
	// Save results table

	statsTitle = "Statistics for " + label + "_thresholded.tif";
	if (isOpen(statsTitle)) {
    	selectWindow(statsTitle);
    	wait(100);
    	saveAs("Table", resultsFolder + File.separator + label + "_3DResults.csv");
    	close();
	}
	
	// Close original channel window if still open
	if (isOpen(label)) close();
}   

