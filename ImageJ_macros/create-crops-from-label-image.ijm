// ImageJ macro to save 3D crops for individual objects from multi-channel hyperstack
// that also crops label images and applies binary masks

#@ File (label = "Directory of original fluorescence images", style = "directory") input
#@ File (label = "Directory of label images", style = "directory") dir_labels
#@ File (label = "Output directory for cropped and masked fluorescence images", style = "directory") outputDir
#@ File (label = "Output directory for cropped binary mask images", style = "directory") outputDirMasks
#@ File (label = "Output directory for cropped label images", style = "directory") outputDirLabels

#@ String (label = "Input file suffix", value = ".tif") suffix


setBatchMode(true);

// Output directories should have the proper delimiter for cleaner code later
if (!endsWith(outputDir, File.separator)) {
    outputDir = outputDir + File.separator;
}
if (!endsWith(outputDirMasks, File.separator)) {
    outputDirMasks = outputDirMasks + File.separator;
}
if (!endsWith(outputDirLabels, File.separator)) {
    outputDirLabels = outputDirLabels + File.separator;
}

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	dir_labels_list = getFileList(dir_labels);
	dir_labels_list = Array.sort(dir_labels_list);
	
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, outputDir, list[i], dir_labels_list[i]);
	}
}

function processFile(input, output, file, label_img) {
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
	originalImageTitle = getTitle();
	
	print(dir_labels + label_img);
	open(dir_labels + File.separator + label_img);
	labelImageTitle = getTitle();
	// Select and analyze label image to find bounding boxes
	selectWindow(labelImageTitle);
	
	if (nSlices > 1) {
	    Stack.getStatistics(_, _, _, labelMax, _);
	} else {
	    getStatistics(_, _, _, labelMax, _, _);
	}
	
	print("Found " + labelMax + " labeled objects");
	
	// Get dimensions of label image
	labelWidth = getWidth();
	labelHeight = getHeight();
	labelSlices = nSlices;
	
	// Create arrays to store bounding box coordinates for each object
	minX = newArray(labelMax + 1);
	maxX = newArray(labelMax + 1);
	minY = newArray(labelMax + 1);
	maxY = newArray(labelMax + 1);
	minZ = newArray(labelMax + 1);
	maxZ = newArray(labelMax + 1);
	
	// Initialize bounding box arrays
	for (i = 0; i <= labelMax; i++) {
	    minX[i] = labelWidth;
	    maxX[i] = -1;
	    minY[i] = labelHeight;
	    maxY[i] = -1;
	    minZ[i] = labelSlices;
	    maxZ[i] = -1;
	}
	
	// Find bounding boxes for each labeled object
	print("Calculating bounding boxes...");
	for (z = 1; z <= labelSlices; z++) {
	    setSlice(z);
	    for (y = 0; y < labelHeight; y++) {
	        for (x = 0; x < labelWidth; x++) {
	            pixelValue = getPixel(x, y);
	            if (pixelValue > 0 && pixelValue <= labelMax) {
	                label = pixelValue;
	                // Update bounding box
	                if (x < minX[label]) minX[label] = x;
	                if (x > maxX[label]) maxX[label] = x;
	                if (y < minY[label]) minY[label] = y;
	                if (y > maxY[label]) maxY[label] = y;
	                if (z < minZ[label]) minZ[label] = z;
	                if (z > maxZ[label]) maxZ[label] = z;
	            }
	        }
	    }
	    print("Processed slice " + z + " of " + labelSlices);
	}
	
	// Switch to original hyperstack
	selectWindow(originalImageTitle);
	
	// Get hyperstack dimensions
	Stack.getDimensions(stackWidth, stackHeight, stackChannels, stackSlices, stackFrames);
	print("Original hyperstack: " + stackWidth + "x" + stackHeight + "x" + stackChannels + "C x" + stackSlices + "Z x" + stackFrames + "T");
	
	// Get base name for output files
	baseName = substring(originalImageTitle, 0, indexOf(originalImageTitle, "."));
	if (baseName == "") baseName = originalImageTitle;
	
	// Create crops for each object
	validObjects = 0;
	for (label = 1; label <= labelMax; label++) {
	    // Check if object exists (has valid bounding box)
	    if (maxX[label] >= 0 && maxY[label] >= 0 && maxZ[label] >= 0) {
	        validObjects++;
	        
	        // Calculate crop dimensions
	        cropX = minX[label];
	        cropY = minY[label];
	        cropZ = minZ[label];
	        cropWidth = maxX[label] - minX[label] + 1;
	        cropHeight = maxY[label] - minY[label] + 1;
	        cropDepth = maxZ[label] - minZ[label] + 1;
	        
	        print("Object " + label + ": " + cropWidth + "x" + cropHeight + "x" + cropDepth + " at (" + cropX + "," + cropY + "," + cropZ + ")");
	        
	        // Create crop for all channels and time points from original hyperstack
	        selectWindow(originalImageTitle);
	        makeRectangle(cropX, cropY, cropWidth, cropHeight);
	        run("Duplicate...", "title=crop_" + label + " duplicate");
	        selectWindow("crop_"+label);       
	        
	        // Crop in Z dimension if needed
	        if (cropZ > 1 || cropZ + cropDepth - 1 < stackSlices) {
	            run("Make Substack...", "channels=1-" + stackChannels + " slices=" + cropZ + "-" + (cropZ + cropDepth - 1) + " frames=1-" + stackFrames);
	            close("crop_" + label);
	            
	            if (stackChannels < 2){
	            selectWindow("Substack (" + cropZ + "-" + (cropZ + cropDepth - 1)+")");
	            }
	            else { // TODO

	            }
	            rename("crop_" + label);
	        }
	        
	        channels = newArray(stackChannels);
	        // If needed, split by channel
	        if (stackChannels > 1) {
	        	selectWindow("crop_" + label);
	        	run("Split Channels");
		        
		        for (c = 1; c <= stackChannels; c++) {
		        	channels[c-1] = "C" + c + "-" + "crop_" + label;	
		        }
	        }
	        else { channels[0] = "crop_" + label;
	        }
	        
	        // Create crop from label image with same dimensions
	        selectWindow(labelImageTitle);
	        makeRectangle(cropX, cropY, cropWidth, cropHeight);
	        run("Duplicate...", "title=label_crop_" + label + " duplicate");
	        selectWindow("label_crop_"+label);
	        
	        // Crop label in Z dimension if needed
	        if (cropZ > 1 || cropZ + cropDepth - 1 < labelSlices) {
	            run("Make Substack...", "slices=" + cropZ + "-" + (cropZ + cropDepth - 1));
	            close("label_crop_" + label);
	            
	            selectWindow("Substack (" + cropZ + "-" + (cropZ + cropDepth - 1)+")");
	            rename("label_crop_" + label);
	        }
	        
	        // Create binary mask by thresholding the label crop
	        selectWindow("label_crop_" + label);
	        run("Duplicate...", "title=mask_" + label + " duplicate");
	        selectWindow("mask_" + label);
	        
	        // Threshold to create binary mask for this specific label
	        setThreshold(label, label);
	        run("Convert to Mask", "background=Dark black");
	        
	        // Convert binary mask to 0s and 1s by dividing by 255
	        run("Divide...", "value=255 stack");
	        
	        for (c=1; c <= stackChannels; c++){
	        	
				// Apply mask to the original crop
		        imageCalculator("Multiply create stack", channels[c-1], "mask_" + label);
		        selectWindow("Result of " + channels[c-1]);
		        masked_obj_name_c = "C"+c+"-"+baseName + "_object_" +IJ.pad(label, 3) + ".tif";
		        rename(masked_obj_name_c);
		        
		        // Save the masked crop
		        outputPath = outputDir + masked_obj_name_c;
		        saveAs("Tiff", outputPath);
		        
		        // Clean up intermediate images

		        close(masked_obj_name_c);		        
		        print("Saved: " + outputPath);
		        
	        }
	        // Save the binary mask
	        selectWindow("mask_" + label);
	        mask_name_out = baseName + "_mask_" + IJ.pad(label, 3) + ".tif";
	        maskOutputPath = outputDirMasks + mask_name_out;
	        saveAs("Tiff", maskOutputPath);
	        
	        // Save the label crop
	        selectWindow("label_crop_" + label);
	        label_name_out = baseName + "_label_" + IJ.pad(label, 3) + ".tif";
	        labelOutputPath = outputDirLabels + label_name_out;
	        saveAs("Tiff", labelOutputPath);
	        
	        print("Saved: " + maskOutputPath);
	        print("Saved: " + labelOutputPath);
	        
	        close(label_name_out);
	        close(mask_name_out);
	    }
	}
	close(labelImageTitle);
	
	print("Processing complete!");
	print("Created " + validObjects + " object crops in: " + outputDir);
		
}

setBatchMode(false);