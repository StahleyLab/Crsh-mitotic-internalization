// ImageJ macro to filter a 3D image by keeping only specified pixel values

#@ File (label = "CSV file containing tidy dataset with filenames, integers, mitotic phase", style = "file") csv_file
#@ File (label = "Directory of label images", style = "directory") dir_labels
#@ File (label = "Output directory of filtered label images", style="directory") dir_filtered

setBatchMode(true);

// feels a little gross to use global variables
var filenames, labelGroups;

// Output directories should have the proper delimiter for cleaner code later
if (!endsWith(dir_filtered, File.separator)) {
    dir_filtered = dir_filtered + File.separator;
}

function GetFilenamesAndIntegersFromCSV(csv_file) { 
	// Read a CSV file containing a tidy dataset
	// with headers `filename`, `label_integer`, `mitotic_phase`
	// Combine the different integers into one array
	
	// Read the CSV file
	fileContent = File.openAsString(csv_file);
	lines = split(fileContent, "\n");
	
	// Skip header line and process data
	filenames = newArray();
	labelGroups = newArray();
	
	// Process each line (starting from line 1 to skip header)
	for (i = 1; i < lines.length; i++) {
	    if (lines[i] != "") {
	        parts = split(lines[i], ",");
	        filename = parts[0];
	        labelInteger = parts[1];
	        
	        // Find if filename already exists in our arrays
	        existingIndex = -1;
	        for (j = 0; j < filenames.length; j++) {
	            if (filenames[j] == filename) {
	                existingIndex = j;
	                break;
	            }
	        }
	        
	        // If filename exists, append the label to existing group
	        if (existingIndex >= 0) {
	            labelGroups[existingIndex] = labelGroups[existingIndex] + "," + labelInteger;
	        } else {
	            // If filename is new, add it to arrays
	            filenames = Array.concat(filenames, filename);
	            labelGroups = Array.concat(labelGroups, labelInteger);
	        }
	    }
	}
	
	//return newArray(filenames, labelGroups);
}

GetFilenamesAndIntegersFromCSV(csv_file);

filenames_arr = filenames;
label_arr = labelGroups;
//file_label_array = GetFilenamesAndIntegersFromCSV(csv_file);
//filenames_arr = file_label_array[0];
//label_arr = file_label_array[1];

processFolder(filenames_arr, label_arr);

function processFolder(filenames, labels) {
	//list = getFileList(input);
	//filenames = Array.sort(filenames);
	//dir_labels_list = getFileList(dir_labels);
	//dir_labels_list = Array.sort(dir_labels_list);
	
	for (i = 0; i < filenames.length; i++) {
		processFile(dir_labels, dir_filtered, filenames[i], labels[i]);
	}
}

function processFile(input, output, file, labels) {
	print(input + File.separator + file);
	open(input + File.separator + file);
	
	if (getImageID == 0) { // If the image wasn't opened (getImageID returns 0 for no image)
	    print("Default open() failed for: " + file);
	    // Fallback to Bio-Formats
	    run("Bio-Formats Importer", "open=[" + file + "] windowless=true");
	} 
	else {
	    print("File opened successfully with default open(): " + file);
	}
	
	// Get image info
	if (nSlices > 1) {
	    Stack.getStatistics(_, _, _, max, _);
	}
	else {
	    getStatistics(_, _, _, max, _, _);
	}
	maxValue = max;
	
	print("Maximum pixel value in the image: " + maxValue);
	
	keepValuesString = labels;
	
	// Parse the comma-separated string into an array
	keepValues = split(keepValuesString, ",");
	keepArray = newArray(keepValues.length);
	
	// Convert strings to integers and store in array
	for (i = 0; i < keepValues.length; i++) {
	    keepArray[i] = parseInt(keepValues[i]);
	}
	
	// Print the values that will be kept
	print("Values to keep: " + keepValuesString);
	
	// Create Boolean array for fast lookup - determines if value i should be changed to 0
	changeToZero = newArray(maxValue + 1);
	for (i = 0; i <= maxValue; i++) {
	    changeToZero[i] = true; // Default: change to 0
	}
	
	// Mark values to keep as false (don't change to 0)
	for (i = 0; i < keepArray.length; i++) {
	    if (keepArray[i] >= 0 && keepArray[i] <= maxValue) {
	        changeToZero[keepArray[i]] = false;
	    }
	}
	
	// Get original image title for saving
	originalTitle = getTitle();
	dotIndex = indexOf(originalTitle, ".");
	if (dotIndex > 0) {
	    baseName = substring(originalTitle, 0, dotIndex);
	    extension = substring(originalTitle, dotIndex);
	} else {
	    baseName = originalTitle;
	    extension = "";
	}
	
	// Loop through all slices
	for (slice = 1; slice <= nSlices; slice++) {
	    setSlice(slice);
	    print("Processing slice " + slice + " of " + nSlices);
	    
	    getStatistics(area, mean, sliceMin, sliceMax, std, histogram);
	    
	    // Loop through possible pixel values in slice
	    for (i = sliceMin; i <= sliceMax; i++) {
	        if (changeToZero[i]) {
	            changeValues(i, i, 0);
	        }
	    }
	}
	
	// Save the modified image with suffix
	newTitle = baseName + "_filtered" + extension;
	
	run("glasbey on dark ");
	saveAs("Tiff", dir_filtered + newTitle);
	
	print("Image saved as: " + newTitle);
	run("Close All");
}

setBatchMode(false);
print("Processing complete!");