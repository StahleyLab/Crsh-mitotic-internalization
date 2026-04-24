/*
 * 2026-04-07
 * William Giang
 * 
 * for Sarah Latario
 * 
 * Goal(s):
 *     1. Determine how many Celsr1+ puncta are Fz6+ 
 *     
 * Context:
 * 		Multi-channel z-stacks were acquired of COS7 cells expressing Vang2 (GFP) (C1) along with 
 * 		either Celsr-WT (mCherry) or Celsr-Crsh (mCherry) (C2).
 * 		
 * 		ilastik's Pixel Classification workflow was used to segment foreground from background
 * 		for the Celsr1 and Vangl2 channels
 * 		
 * 		An additional Object Classification step was used to
 * 		- split up nearby puncta
 * 		- enhance segmentation selection
 * 		 		
 * Macro plan: 
 * 	1. Load fluorescence images and their segmentation results
 * 	2. For Celsr1 and Vangl2, threshold to select correct class
 * 	3. Size filter and connected components analysis
 * 	3. Convert label images to binary images with values of 0 and 1
 * 	4. Mask fluorescence images by multiplying fluorescence images with corresponding binary masks
 * 	5. Using the label images, measure the integrated/mean/max intensities of the masked fluorescence images
 * 	6. Save CSV files
 * 	  
 * After: Data wrangle with Python
 *  - Celsr1 puncta have measurements for Celsr1 and (masked) Fz6
 *  - Fz6 puncta have measurements for Fz6
 *  - A non-zero max intensity for the masked channel suggests colocalization
 * 
 * Input: 
 * 	- 8-bit image of object predictions
 * 	  threshold to keep only pixels with a value of 2.
 * 	- hand-drawn selection on multi-channel fluorescense MIP
 * Output: csv files (and if desired, images showing intermediate steps)
 * 
 * Assumes nothing else is in the folders but the images
 * Assumes that the desired object class is second for all channels
 */


#@ File    (label = "Input fluorescence directory", style = "directory") input_dir_fluor_1
#@ File    (label = "Input segmentation CH 1", style = "directory") input_dir_seg_1
#@ File    (label = "Input segmentation CH 2", style = "directory") input_dir_seg_2
#@ String  (label = "CH1 name", value="Vangl2") CH1_name
#@ String  (label = "CH2 name", value="Celsr1")   CH2_name
#@ File    (label = "Output directory", style = "directory") output

/*
function extractObjectsFromObjectProbabilities(img, n_classes, desired_channel, prob_threshold, max_label_size) {
	
	// input ilastik files are interleaved with number of classes
	selectImage(img);
	run("Deinterleave", "how="+n_classes);
	desired_img = img + " #" + desired_channel;
	
	// threshold by probability to create a binary mask
	selectWindow(desired_img);
	setThreshold(prob_threshold, 1);
	run("Convert to Mask", "background=Dark black create");
	temp_mask = "MASK_" + desired_img;
	
	// Use 26-connected connected components analysis to create 3D objects
	selectImage(temp_mask);
	run("Connected Components Labeling", "connectivity=26 type=[16 bits]");
	tmp_name = getTitle();
	selectWindow(tmp_name);
	sanitycheck(want_sanity_checks);
	
	// Size exclusion filter to remove abnormally large objects (clumped segmentations)
	tmp_name = getTitle();
	selectWindow(tmp_name);
	run("Label Size Filtering", "operation=Lower_Than size="+max_label_size);
	sanitycheck(want_sanity_checks); // note this gets overwritten later after remapping
	
	// For convenience down the road, remap labels
	tmp_name = getTitle();
	selectImage(tmp_name);
	run("Remap Labels");
	tmp_title = getTitle();
	selectImage(tmp_title);
	img_ID = getImageID();
	sanitycheck(want_sanity_checks); // this overwrites the previous sanitycheck after size filtering
	
	return img_ID;
}

function makeMaskOfOverlappingChannels(binary_mask_A, binary_mask_B, AND_result_name, label_result_name) {
	imageCalculator("AND create stack", binary_mask_A, binary_mask_B);
	tmp_name = getTitle();
	selectWindow(tmp_name);
	rename(AND_result_name);
}

function makeLabelImageUsingOtherChannelAsValues(label_image_to_get_values_from, binary_mask, matched_label_image_name) {
	selectWindow(binary_mask);
	run("Duplicate...", "title=" + matched_label_image_name + " duplicate");
	selectWindow(matched_label_image_name);
	run("Divide...", "value=255.000 stack");
	setOption("ScaleConversions", true);
	run("16-bit");
	tmp_name = getTitle();
	imageCalculator("Multiply create stack", tmp_name, label_image_to_get_values_from);
	tmp = getTitle();
	close(tmp_name);
	selectWindow(tmp);
	rename(matched_label_image_name);
}

function getAndSetMultiColoc(img_a, img_b, main_row_to_write, table) {
	img_a = removeTifExtension(img_a);
	img_b = removeTifExtension(img_b);
	run("3D MultiColoc", "image_a="+img_a + " image_b=" + img_b);
	
	row_to_write = main_row_to_write;
	tmp_table_name = "Colocalisation";
	selectWindow(tmp_table_name);
	for (j = 0; j < Table.size; j++) {
		LabelObj = Table.get("LabelObj", j, tmp_table_name);
		O1 =       Table.get("O1",       j, tmp_table_name);
		V1 =       Table.get("V1",       j, tmp_table_name);
		P1 =       Table.get("P1",       j, tmp_table_name);
		
		Table.set("LabelObj", row_to_write, LabelObj, table);
		Table.set("O1",       row_to_write, O1,       table);
		Table.set("V1",       row_to_write, V1,       table);
		Table.set("P1",       row_to_write, P1,       table);
		
		row_to_write += 1;
	}
	Table.update;
	close(tmp_table_name);
	return row_to_write;
	
}

function thresholdChannelAsMask(img, threshold_val, mask_name) {
	// threshold a channel 
	selectWindow(img);
	setThreshold(threshold_val, 1.0);
	setOption("BlackBackground", true);
	run("Convert to Mask", "background=Dark black");
	tmp = getTitle();
	selectWindow(tmp);
	rename(mask_name);
}
*/
function removeTifExtension(img_name) {
	if (endsWith(img_name, ".tif")) {
		return substring(img_name,0,lengthOf(img_name)-4);
	}
	if (endsWith(img_name, ".tiff")) {
		return substring(img_name,0,lengthOf(img_name)-5);
	}
	else {
		return img_name
	}
}
/*
function getAndSetVoxelCountAndMean(input_img, labels_img, input_ch, labels_name, main_row_to_write, table, want_label, want_NVoxels, want_Volume, want_img_name_without_ch) {
	run("Intensity Measurements 2D/3D", "input=" + input_img + " labels=" + labels_img + " mean stddev max numberofvoxels volume");
	input_img_name_without_ext = removeTifExtension(input_img);
	temp_table_name = input_img_name_without_ext + "-intensity-measurements";
	input_img_name_without_ext_or_ch = substring(input_img_name_without_ext, 3, lengthOf(input_img_name_without_ext));

	row_to_write = main_row_to_write;

	for (z = 0; z < Table.size; z++) {
	    label_obj = Table.getString("Label",    z, temp_table_name);
	    Mean_int  = Table.get("Mean",           z, temp_table_name);
	    StdDev    = Table.get("StdDev",         z, temp_table_name);
	    Max       = Table.get("Max",            z, temp_table_name);
	    NVoxels   = Table.get("NumberOfVoxels", z, temp_table_name);
	    Volume    = Table.get("Volume",         z, temp_table_name);
	    
	    Table.set("img", row_to_write, input_img_name_without_ext, table);
	    if (want_img_name_without_ch) Table.set("img", row_to_write, input_img_name_without_ext_or_ch, table);
	    if (want_label)   Table.set("Label_"+labels_name,   row_to_write, label_obj, table);
	    if (want_NVoxels) Table.set("NVoxels_"+labels_name, row_to_write, NVoxels,   table);
	    if (want_Volume)  Table.set("Volume_"+labels_name,  row_to_write, Volume,    table);
	    Table.set("Mean_"+input_ch,   row_to_write, Mean_int, table);
	    Table.set("StdDev_"+input_ch, row_to_write, StdDev,   table);
	    Table.set("Max_"+input_ch,    row_to_write, Max,      table);

	    row_to_write += 1;
	    Table.update;
	}
	close(temp_table_name);
	return row_to_write;
}

function getIntensityForLabelImageAndSave(objects_image, signal_image, ch_suffix) {
	//print("objects_image: " + objects_image);
	//print("signal_image: " + signal_image);
	// remove file extension since the table name strips it
	//objects_image = removeTifExtension(objects_image);
	signal_image  = removeTifExtension(signal_image);
	
	run("3D Intensity Measure", "objects=" + objects_image + " signal=" + signal_image);
	if (nResults > 0){
			// rename specific columns	
		col_names = newArray("IntensityMin","IntensityMax","IntensityAvg","IntensitySum");
		for (c = 0; c < col_names.length; c++) {
			Table.renameColumn(col_names[c], col_names[c]+ch_suffix, "Results");
		}
		Table.update;
		
		tmp_save_name = objects_image + "_" + signal_image + "_" + ch_suffix + ".csv";
		
		saveAs("Results", output + File.separator + tmp_save_name);
		close(tmp_save_name);
		close("Results");
	}
}
*/
function getIntensityAndSave2D(label_image, signal_image, ch_suffix) {
	run("Intensity Measurements 2D/3D",
	"input="+signal_image+" labels="+label_image + " mean max numberofvoxels");
	
	signal_image_no_ext = removeTifExtension(signal_image);
	tmp_save_name = label_image + "_" + signal_image_no_ext + "_" + ch_suffix + ".csv";
	saveAs("Results", output + File.separator + tmp_save_name);
	close(tmp_save_name);
	close("Results");

	
}

/*
function makeLabelImageThenSizeFilterAndRemap(binary_mask_image, min_label_size, label_image_name) {
	selectWindow(binary_mask_image);
	run("Connected Components Labeling", "connectivity=26 type=[16 bits]");
	run("Label Size Filtering", "operation=Greater_Than size="+min_label_size);
	run("Remap Labels");
	remapped_orig_name = getTitle();
	remapped_orig_name_no_ext = removeTifExtension(remapped_orig_name);
	selectWindow(remapped_orig_name);
	rename(label_image_name);
	
	// clean up intermediate image
	name_without_ext = removeTifExtension(binary_mask_image);
	close(name_without_ext+"-lbl");
}
*/

function makeBinaryMaskFromLabelImage(label_image) { 
	selectImage(label_image);
	run("Duplicate...", "title="+label_image+"_MASK");
	setThreshold(1, 65535, "raw");
	run("Convert to Mask", "background=Dark black create");
	//resulting image has a "MASK_" prefix
	
	new_title = getTitle();
	// reset things just in case
	selectImage(label_image);
	resetThreshold;
	return new_title;
}

function selectObjectClassificationClassAndLabel(image, integer, connectivity) {
	selectWindow(image);
	setThreshold(integer, integer, "raw");
	run("Convert to Mask", "background=Dark black create");
	run("Connected Components Labeling", "connectivity="+connectivity+" type=[16 bits]");
	//resulting image has a "MASK_" prefix
	new_title = getTitle();
	
	// reset things just in case
	selectWindow(image);
	resetThreshold;

	return new_title;
}

function getAndSortFiles(directory) {
	filelist = getFileList(directory);
	
	return Array.sort(filelist);
}

function getMaskedImage(mask_img, img_to_mask) {
	selectImage(mask_img);
	run("Divide...", "value=255.000 stack");
	imageCalculator("Multiply create 32-bit stack", mask_img, img_to_mask);
	desired_name = "Masked_" + img_to_mask;
	rename(desired_name);
	
	return desired_name;
}

function getNameUsingImageID(image_id) { 
	selectImage(image_id);
	return getTitle();
}

function clearOutsideROI(image_id) {
	selectImage(image_id);
	roiManager("Select", 0);
	setBackgroundColor(0, 0, 0);
	run("Clear Outside", "stack");
}

setBatchMode(true);

list_input1 = getAndSortFiles(input_dir_fluor_1);

list_input_seg1 = getAndSortFiles(input_dir_seg_1);
list_input_seg2 = getAndSortFiles(input_dir_seg_2);

for (i = 0; i < list_input1.length; i++) {
//for (i = 0; i < 2; i++) { //testing
	// Fresh start by clearing Results table and ROI manager
	run("Fresh Start");
	
	// load and if needed, assign physical dimensions
	//print(list_input1[i]);
	open(input_dir_fluor_1 + File.separator + list_input1[i]);
	fluor = getTitle();
	
	getVoxelSize(width, height, depth, unit);
	
	// Add selection from MIP
	run("Add Selection..."); // adds selection to overlay
	run("To ROI Manager");   // moves overlay to ROI Manager
	
	// isolate the fluorescent channels...could split but whatever
	selectImage(fluor);
	run("Duplicate...", "duplicate channels=1");
	c1_fluor = getTitle();
	selectImage(fluor);
	run("Duplicate...", "duplicate channels=2");
	c2_fluor = getTitle();
	
	
	// open segmentation results, assign spatial calibration, clear outside using ROI
	open(input_dir_seg_1 + File.separator + list_input_seg1[i]);
	setVoxelSize(width, height, depth, unit);
	c1_seg_orig = getImageID();
	clearOutsideROI(c1_seg_orig);
	c1_seg_ID = getImageID();
	c1_seg_name = getTitle();
	
	open(input_dir_seg_2 + File.separator + list_input_seg2[i]);
	setVoxelSize(width, height, depth, unit);
	c2_seg_orig = getImageID();
	clearOutsideROI(c2_seg_orig);
	c2_seg_ID = getImageID();
	c2_seg_name = getTitle();
	
	c1_seg_label = selectObjectClassificationClassAndLabel(c1_seg_name, 2, 8);
	c2_seg_label = selectObjectClassificationClassAndLabel(c2_seg_name, 2, 8);
	

	// Get binary masks from thresholded and filtered label images	
	c1_mask_ID = makeBinaryMaskFromLabelImage(c1_seg_label);	
	c2_mask_ID = makeBinaryMaskFromLabelImage(c2_seg_label);

	// prepare image window titles without extensions
	c1_labels_name = getNameUsingImageID(c1_seg_label);
	c2_labels_name = getNameUsingImageID(c2_seg_label);
	
	// Mask fluorescence channels
	c1_fluor_masked = getMaskedImage(c1_mask_ID, c1_fluor);
	c2_fluor_masked = getMaskedImage(c2_mask_ID, c2_fluor);
	
	c1_fluor_masked_no_ext = removeTifExtension(c1_fluor_masked);
	c2_fluor_masked_no_ext = removeTifExtension(c2_fluor_masked);
	
	// Measure intensities on masked fluorescence images
	// individual puncta
	//getIntensityAndSave2D(c1_labels_name, c1_fluor_masked, CH1_name+"-"+CH1_name);
	//getIntensityAndSave2D(c2_labels_name, c2_fluor_masked, CH2_name+"-"+CH2_name);
	
	// for Celsr+ puncta, is there also Vang?
	getIntensityAndSave2D(c2_labels_name, c1_fluor_masked, CH2_name+"-"+CH1_name);

	
	run("Close All");
}

setBatchMode(false);