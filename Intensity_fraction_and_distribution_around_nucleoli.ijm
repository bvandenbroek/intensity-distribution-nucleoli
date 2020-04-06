/*
 * Macro to quantify the intensity fraction of a signal (here centromeres) inside nucleoli,
 * versus outside, and calculate the 3D intensity distribution outside the nucleoli.
 * 
 * Bram van den Broek, Netherlands Cancer Institute, 2019
 * b.vd.broek@nki.nl
 * 
 * Required input: a folder containing 3D images with 3 channels: nuclei, nucleoli and signal of interest.
 * 
 * Workflow summary:
 * 
 * 1. Cell nuclei segmentation in 2D by applying:
 *    - rolling ball background subtraction (radius dependent on estimated nucleus size)
 *    - 3D median filter (currently 2 pixels)
 *    - maximum intensity projection
 *    - Otsu hresholding
 *    - Distance transform watershed operation to separate touching nuclei.
 *    - Manual fixing of segmentation mistakes by deletion/addition/combining of ROIs.
 * 
 * 2. 3D segmentation of nucleoli using 3D Object Counter.
 *    - For every nucleus, an individual threshold for 3D segmentation is set equal to the
 *      Otsu threshold level of the (background-subtracted) 2D maximum intensity projection
 *      of the nucleoli.
 * 
 * 3. 3D segmentation of the centromeres.
 *    - A masked centromere stack (Background set to NaN) is created using the a similar strategy,
 *      except that here a single threshold is determined for the centromere signal in all valid
 *      nuclear ROIs in the image. We found this method more robust than individual thresholds.
 *      Thus, voxels with gray levels below 0.5 times the Otsu threshold are set to NaN.
 *    - The found threshold can be adjusted manually if unsatisfactory.
 * 
 * 4. For every nucleus, determine the probability density of centromere signal, as a function
 *    of the distance to the nucleoli edge.
 *    - Sets of 3D masks ("shells") around the nucleoli with certain thickness (e.g. 200 nm) are created
 *      by using the 3D distance map plugin. (Ollion et al., Bioinformatics 2013)
 *    - These distance masks are multiplied by the background-removed centromere stack of that nucleus
 *      to get the average intensity at different distances from the nucleoli edge.
 *    - The measurements were normalized by dividing by the total centromere intensity,
 *      to obtain the intensity distribution around the nucleoli.
 * 
 * 5. Display result tables, overlays and plots for visual inspection.
 */

#@ File (label = "Input folder", style = "directory") input
#@ File (label = "Output folder", style = "directory") output
#@ String (label = "Files to analyze", value = ".D3D.dv") fileExtension

#@ Integer (label = "Crop edges of the image (pixels)", value = 32) cropBorder
#@ Integer (label = "Nuclei channel", value = 1) chNuc
#@ Integer (label = "Nucleoli channel", value = 2) chNucleoli
#@ Integer (label = "Centromer channel", value = 3) chCentromers
#@ Integer (label = "Estimated nucleus size (um)", value = 8) nucSize
#@ String (label = "Nuclei threshold method", choices={"Global", "Local"}, style="radioButtonHorizontal", value = "local") thresholdMethod
#@ Boolean (label = "Automatic threshold (global only, local is always auto)", value=true) autoThreshold
#@ Boolean (label = "Exclude nuclei on edges of the image", value=true) excludeEdges
#@ Boolean (label = "Manually edit segmented nuclei?", value = true) editNuclei
#@ Boolean (label = "Automatic threshold on centromers?", value = true) autoCentromerThreshold

#@ Float (label = "Analyze 3D intensity profile up to max distance (um)", value = 5) maxDist
#@ Float (label = "with increments of (um)", value = 0,2) increment

#@ Boolean (label = "Display images while analyzing?", value=false) display_images

if(editNuclei == true) display_images=true;

threshold_nuc_bias = 0;
watershed = true;

nBins = maxDist/increment;
probDensityTable = "3D_profile";
Table.create(probDensityTable);

//Initialize
saveSettings();

print("\\Clear");
if(nImages>0) run("Close All");
run("Clear Results");
roiManager("reset");
setOption("BlackBackground", true);
setBackgroundColor(0,0,0);
run("Conversions...", " ");
run("Set Measurements...", "area mean standard integrated median redirect=None decimal=4");

var nrOfImages=0;
var current_image_nr=0;
var processtime=0;
outputSubfolder = output;	//initialize this variable


if(!File.exists(output)) {
	create = getBoolean("The specified output folder "+output+" does not exist. Create?");
	if(create==true) File.makeDirectory(output);		//create the output folder if it doesn't exist
	else exit;
}

setBatchMode(true);

scanFolder(input);
processFolder(input);

restoreSettings();





//////////  FUNCTIONS  //////////

// function to scan folders/subfolders/files to count files with correct fileExtension
function scanFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			scanFolder(input + File.separator + list[i]);
		if(endsWith(list[i], fileExtension))
			nrOfImages++;
	}
}



// function to scan folders/subfolders/files to find files with correct fileExtension
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) {
			outputFolder = output + File.separator + list[i];	
			if(!File.exists(outputSubfolder)) File.makeDirectory(outputSubfolder);	//create the output subfolder if it doesn't exist
			processFolder(input + File.separator + list[i]);
		}
		if(endsWith(list[i], fileExtension)) {
			current_image_nr++;
			showProgress(current_image_nr/nrOfImages);
			processFile(input, outputSubfolder, list[i]);
		}
	}
//	print("\\Clear");
	print("\\Update1:Finished processing "+nrOfImages+" files.");
	print("\\Update2:Average speed: "+d2s(current_image_nr/processtime,1)+" images per minute.");
	print("\\Update3:Total run time: "+d2s(processtime,1)+" minutes.");
	print("\\Update4:-------------------------------------------------------------------------");

}



function processFile(input, outputSubfolder, file) {
	if(nImages>0) run("Close All");
	roiManager("Reset");
	print("\\Clear");
	
	starttime = getTime();
	print("\\Update1:Processing file "+current_image_nr+"/"+nrOfImages+": " + input + file);
	print("\\Update2:Average speed: "+d2s((current_image_nr-1)/processtime,1)+" images per minute.");
	time_to_run = (nrOfImages-(current_image_nr-1)) * processtime/(current_image_nr-1);
	if(time_to_run<5) print("\\Update3:Projected run time: "+d2s(time_to_run*60,0)+" seconds ("+d2s(time_to_run,1)+" minutes).");
	else if(time_to_run<60) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. You'd better get some coffee.");
	else if(time_to_run<480) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes ("+d2s(time_to_run/60,1)+" hours). You'd better go and do something useful.");
	else if(time_to_run<1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. ("+d2s(time_to_run/60,1)+" hours). You'd better come back tomorrow.");
	else if(time_to_run>1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. This is never going to work. Give it up!");
	print("\\Update4:-------------------------------------------------------------------------");

	name = substring(file,0,lastIndexOf(file, "."));	//filename without extension
	name = replace(name,"\\/","-");	//replace slashes by dashes in the name
	name = replace(name," ","_");	//replace slashes by dashes in the name

	//START OF ANALYSIS-SPECIFIC CODE
	run("Bio-Formats Importer", "open=["+input + File.separator + file+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
	run("32-bit");	//Get rid of stupid Deltavision calibration
	run("16-bit");
	getDimensions(width, height, channels, slices, frames);
	if(cropBorder>0) {
		makeRectangle(cropBorder, cropBorder, width-2*cropBorder, height-2*cropBorder);
		run("Crop");
	}
	getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(unit=="microns") unit="micron";

	original = getTitle();
	selectWindow(original);
	if(display_images==true) {
		Stack.setDisplayMode("composite");
		Stack.setSlice(floor(slices/2));
		Stack.setChannel(3);
		run("Cyan Hot");
		resetMinAndMax;
		resetMinAndMax;
		resetMinAndMax;
		Stack.setChannel(2);
		resetMinAndMax;
		Stack.setChannel(1);
		resetMinAndMax;
		updateDisplay;
		setBatchMode("show");
	}

	nucleusArea = segmentNuclei(original);	//segment nuclei and return the area as an array
	nrNuclei = nucleusArea.length;

	//Determine the size of the largest nucleus. All nucleus images will be resized to this format
	roiWidthMax = 0; roiHeightMax = 0;
	for(i=0;i<nrNuclei;i++) {
		roiManager("select",i);
		getSelectionBounds(roiX, roiY, roiWidth, roiHeight);
		roiWidthMax = maxOf(roiWidthMax,roiWidth);
		roiHeightMax = maxOf(roiHeightMax,roiHeight);
	}
	concatenationString = "";


	//Create masked Centromers stack
	selectWindow(original);
	run("Select None");
	run("Duplicate...", "duplicate title=Centromers_masked channels="+chCentromers);
	run("32-bit");
	run("Subtract Background...", "rolling=10 sliding stack");	//Subtract background in Centromer channel
	run("Grays");
	run("Z Project...", "projection=[Max Intensity]");	//Determine threshold on Max projection...
	setAutoThreshold("Otsu dark");
	roiManager("Select All");
	if(nrNuclei>1) roiManager("Combine");
	setAutoThreshold("Otsu dark stack");				//...and only in the combined nuclei selection
	getThreshold(CentromerThreshold, upperCentromers);
	close();
	selectWindow("Centromers_masked");
	CentromerThreshold = CentromerThreshold/2;			//Set threshold to half the one from the max projection
	setThreshold(CentromerThreshold, upperCentromers);	
	if(autoCentromerThreshold==false) {
		setBatchMode("show");
		run("Threshold...");
		waitForUser("Adjust lower threshold for centromer signal and press OK");
		getThreshold(CentromerThreshold, upperCentromers);
		setBatchMode("hide");
	}
	print("\\Update5:Centromers threshold: "+CentromerThreshold);
	run("NaN Background", "stack");						//Masked Centromers 3D
	resetThreshold();

	//Declare data containers and initalize with zeros.
	CentromerSignalTotal = newArray(nrNuclei);
	CentromerSignalInNucleoli = newArray(nrNuclei);
	nucleoliVolume = newArray(nrNuclei);
	Array.fill(CentromerSignalTotal, 0);
	Array.fill(CentromerSignalInNucleoli, 0);
	Array.fill(nucleoliVolume, 0);

	//Create plot
	xValues=Array.getSequence(nBins);
	xValues=divideArraybyScalar(xValues, 1/increment);
	Table.reset(probDensityTable);
	Table.setColumn("distance (um)", xValues, probDensityTable);


	//Loop over all nuclei in the image
	for(i=0;i<nrNuclei;i++) {
		showStatus("Analyzing nucleus "+i+"...");
		showProgress(i/nrNuclei);
		selectWindow(original);
		roiManager("Select",i);
		run("Duplicate...", "duplicate title=Nucleoli_"+i+" channels="+chNucleoli);
		run("Clear Outside", "stack");
		//Determine threshold
			run("Z Project...", "projection=[Max Intensity]");
			run("Restore Selection");
			setAutoThreshold("Otsu dark");
			getThreshold(lowerNucleoli, upperNucleoli);
			close();
		selectWindow("Centromers_masked");
		roiManager("Select",i);

		run("Duplicate...", "duplicate title=Centromers_"+i+" channels="+chCentromers);
		run("Clear Outside", "stack");
		//Measure Centromer total signal
		run("Measure Stack...");

		for(j=0;j<nResults;j++) {
			if(!isNaN(getResult("RawIntDen",j))) CentromerSignalTotal[i] = CentromerSignalTotal[i] + getResult("RawIntDen",j);	//Add the signal of all slices
		}
		run("Clear Results");

		selectWindow("Nucleoli_"+i);
		//Segment nucleoli (and no redirection of measurements to the Masked Centromers channel any more; yields strange results)
		run("3D OC Options", "volume nb_of_surf._voxels integrated_density mean_gray_value minimum_gray_value maximum_gray_value dots_size=5 font_size=10 white_numbers");
		run("3D Objects Counter", "threshold="+lowerNucleoli+" min.=100 max.=99999 objects statistics summary");

		nucleoliVolumeArray = Table.getColumn("Volume (micron^3)", "Results");
		nucleoliVolume[i] = sumArray(nucleoliVolumeArray);
	
		close("Nucleoli_"+i);
		//Create 3D mask and multiply with Centromers image and measure in the same way as above
		for(z=1;z<=slices;z++) {
			setSlice(z);
			changeValues(1, 255, 1);	//Change all nucleoli volumes to 1
		}
		run("32-bit");
		imageCalculator("Multiply create stack", "Centromers_"+i,"Objects map of Nucleoli_"+i);
		rename("Centromers_in_nucleoli_"+i);
		run("Measure Stack...");
		//waitForUser("Centromer signal in nucleoli "+i);
		for(j=0;j<nResults;j++) {
			if(!isNaN(getResult("RawIntDen",j))) CentromerSignalInNucleoli[i] = CentromerSignalInNucleoli[i] + getResult("RawIntDen",j);
		}
		run("Clear Results");


		//Calculate probability density of centromer intensity with respect to the distance to the nucleoli
		intensity=newArray(nBins);
		volume=newArray(nBins);
		k=0;
		m=0;

		selectWindow("Centromers_"+i);
		setSlice(slices/2);
		run("16-bit");
		//measure total centromers intensity
		centromerInt = 0;
		for(p=0;p<slices;p++) {
				setSlice(p+1);
				List.setMeasurements;
				centromerInt = centromerInt + List.getValue("RawIntDen");
		}

		selectWindow("Objects map of Nucleoli_"+i);
		rename("nucleoli_"+i);
		setSlice(slices/2);
		setMinAndMax(0, 1);
		run("3D Distance Map", "map=EDT image=nucleoli_"+i+" mask=Same threshold=0.5 inverse");
		setSlice(slices/2);
		setMinAndMax(0, 3);
		run("Fire");
		
		//Discretize by multiplying by increment, converting to 16-bit, back to 32-bit, and divide by increment again.
		run("Divide...", "value="+increment+" stack");
		run("16-bit");
		run("32-bit");
		run("Multiply...", "value="+increment+" stack");
		run("Enhance Contrast", "saturated=0.35");
		rename("distance_map");

		//Create shell masks from the distance map with increments
		for(n=0; n<nBins ;n++) {
			selectWindow("distance_map");
			run("Duplicate...", "title=shell_"+n*increment+" duplicate");
			for(s=0;s<slices;s++) {
				setSlice(s+1);
				changeValues(0,n*increment-increment/2,0);
				changeValues(n*increment-increment/2,(n+1)*increment-increment/2,-1);
				changeValues((n+1)*increment-increment/2,65536,0);
				changeValues(-1,-1,1);
			}
			setSlice(slices/2);

			//measure total volume of the shell
			run("16-bit");	//Convert NaNs to zeros
			for(j=0;j<slices;j++) {
				setSlice(j+1);
				List.setMeasurements;
				volume[k] = volume[k] + List.getValue("RawIntDen");	//Only 1 and 0, so total intensity = total area in voxels
			}
			k++;
		
			imageCalculator("Multiply stack", "shell_"+n*increment,"Centromers_"+i);
			//waitForUser("masked centromers multiplied by shell "+n*increment);

			//measure total intensity in shell
			for(j=0;j<slices;j++) {
				setSlice(j+1);
				List.setMeasurements;
				intensity[m] = intensity[m] + List.getValue("RawIntDen");
			}
			m++;
			close("shell_"+n*increment);
		}
		close("distance_map");
		normInt = divideArrays(intensity, volume);	//average intensity in shells
		//Array.print(volume);
		//Array.print(intensity);
		//Array.print(normInt);
		sumInt = sumArray(intensity);	//sum of average intensity in all shells
		
		probabilityDensity = divideArraybyScalar(intensity, sumInt);
		//Array.print(probabilityDensity);
		Table.setColumn("nucleus_"+i+1, probabilityDensity, probDensityTable);
		Table.update;
		
		//Prepare merged image with nucleoli (blue), centromers (red) and centromers inside nucleoli (green -> appear as white)
		selectWindow("Centromers_"+i);
		run("32-bit");
		run("Merge Channels...", "c1=Centromers_"+i+" c2=[Centromers_in_nucleoli_"+i+"] c3=[nucleoli_"+i+"] create ignore");
		rename("Nucleus_"+i);
		close("Centromers_"+i);
		close("Centromers_in_nucleoli"+i);
		close("nucleoli"+i);
		Stack.setSlice(floor(slices/2));
		Stack.setChannel(3);
		resetMinAndMax;
		Stack.setChannel(2);
		resetMinAndMax;
		Stack.setChannel(1);
		resetMinAndMax;
		roiManager("Select", i);
		concatenationString += "image"+i+1+"=Nucleus_"+i+" ";
		//if(display_images==true) setBatchMode("show");
	}
	Table.update;
	Table.save(outputSubfolder + File.separator + name + "_3D_profile.tsv");

	//Generate line plot of probability density with random colors
	Plot.create("Probability Density", "distance (um)", "frequency");
	for(i=0;i<nrNuclei;i++) {
		data = Table.getColumn("nucleus_"+i+1);
		color1 = toHex(random*255);
		color2 = toHex(random*255);
		color3 = toHex(random*255);
		Plot.setColor("#"+color1+color2+color3);
		
		Plot.add("line", xValues, data);
	}
	Plot.show();
	Plot.setLimitsToFit();
	setBatchMode("show");
	saveAs("Tiff", outputSubfolder + File.separator + name + "_plot");

	roiManager("Select All");
	roiManager("Remove Channel Info");
	roiManager("Remove Slice Info");
	roiManager("Remove Frame Info");

	selectWindow(original);
	run("From ROI Manager");
	saveAs("Tiff", outputSubfolder + File.separator + name + "_overview");

	//Concatenate all 3D nuclei in this image file
	if(nrNuclei>1) {
		run("Concatenate...", "  title=Nuclei open "+concatenationString);
		run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices="+slices+" frames="+nrNuclei+" display=Composite");
		selectWindow("Nuclei");
		setBatchMode("show");
	}
	else {
		selectWindow("Nucleus_0");
		rename("Nuclei");
		setBatchMode("show");
	}

	//Draw boundaries of nuclei (Note: The image has to be shown first. Otherwise drawing happens on the wrong image. Then hide image again to increase speed)
	setForegroundColor(255, 255, 255);
	close("B&C");
	selectWindow("Nuclei");
	setBatchMode("hide");
	for(i=0;i<nrNuclei;i++) {
		roiManager("select",i);
		getSelectionBounds(x, y, roiWidth, roiHeight);
		setSelectionLocation((roiWidthMax-roiWidth)/2, (roiHeightMax-roiHeight)/2);
		Stack.setFrame(i+1);	//select the correct nucleus in the stack
		for(z=1;z<=slices;z++) {
			Stack.setSlice(z);
			for(c=1;c<channels;c++) {	//Skip blue channel -> drawing becomes yellow
				Stack.setChannel(c);
				run("Draw", "slice");
			}
		}
	}
	run("Select None");
	removeNaNs("Nuclei");
	if(display_images == true) setBatchMode("show");
	saveAs("Tiff", outputSubfolder + File.separator + name + "_nuclei_3D");
	run("RGB Color", "slices frames keep");
	run("Z Project...", "projection=[Max Intensity] all");
	if(display_images == true) setBatchMode("show");
	run("In [+]");
	run("In [+]");
	run("In [+]");
	saveAs("Tiff", outputSubfolder + File.separator + name + "_nuclei_3D_Zprojection");

	fraction = newArray(nrNuclei);
	run("Clear Results");
	for(i=0;i<nrNuclei;i++) {
		fraction[i] = CentromerSignalInNucleoli[i] / CentromerSignalTotal[i];
		setResult("Nucleus",i,i+1);
		setResult("Nucleus_area ("+unit+"^2)",i,nucleusArea[i]);
		setResult("Nucleoli_volume (micron^3)",i,nucleoliVolume[i]);
		setResult("Centr_threshold",i,d2s(CentromerThreshold,1));
		setResult("Centr_int_in_nucleoli",i,d2s(CentromerSignalInNucleoli[i]/1000000,2));
		setResult("Centr_int_total",i,d2s(CentromerSignalTotal[i]/1000000,2));
		setResult("Fraction",i,fraction[i]);
	}
	saveAs("Results", outputSubfolder + File.separator + name + ".tsv");

	endtime = getTime();
	processtime = processtime+(endtime-starttime)/60000;


}


function segmentNuclei(image) {
	selectWindow(image);
	run("Duplicate...", "duplicate title=DAPI channels="+chNuc);
	run("Subtract Background...", "rolling="+nucSize / pixelWidth / 4+" sliding stack");
	run("Median 3D...", "x=2 y=2 z=2");
	run("Z Project...", "projection=[Max Intensity]");
	if(thresholdMethod=="local") {
		run("8-bit");
		setAutoThreshold("Percentile");
		run("Create Selection");
		List.setMeasurements();
		median = List.getValue("Median");
		std = List.getValue("StdDev");
		resetThreshold();
		run("Select None");
		parameter_1=-(4*std);					//empirical choice
		parameter_1=-minOf(-parameter_1,10);	//minimum value
		parameter_1=-maxOf(-parameter_1,2);		//maximum value
		run("Auto Local Threshold", "method=Mean radius="+(nucSize/pixelWidth/1.5)+" parameter_1="+parameter_1+" parameter_2=0 white");
		print("\\Update:");	//remove statement from auto local threshold plugin
		getThreshold(min,max);
		resetThreshold();
		setThreshold(min,max);
	}
	else if (autoThreshold==true) setAutoThreshold("Otsu dark no-reset");
	else {
		print("threshold method: "+thresholdMethod);
		setBatchMode("show");
		run("Threshold...");
		selectWindow("Threshold");
		waitForUser("Set threshold for segmentation of nuclei and press OK");
	}

	run("Convert to Mask", "  black");
	run("Fill Holes");

	if(watershed==true) run("Watershed");
	setThreshold(127, 255);
	if(excludeEdges==true) run("Analyze Particles...", "size="+PI/4*nucSize/3*nucSize/3+"-"+PI/4*nucSize*3*nucSize*3+" circularity=0.20-1.00 show=Nothing display exclude add");
	else run("Analyze Particles...", "size="+PI/4*nucSize/3*nucSize/3+"-"+PI/4*nucSize*3*nucSize*3+" circularity=0.20-1.00 show=Nothing display add");

	//Manually edit detected nuclei
	if (editNuclei==true) edit_ROIs(image);
	
	//Measure nuclei again
	run("Clear Results");
	roiManager("Measure");
	nucleusArea = Table.getColumn("Area", "Results");

	for(j=0;j<roiManager("count");j++) {
		roiManager("select", j);
		roiManager("rename", "nucleus_"+j+1);
	}
	return nucleusArea;
}


function thresholdSignal(image) {	//Determine the threshold on the slice with the highest intensity. Alternative: just take the MAX projection.
	selectWindow(image);			//The signal image with the nucleus selection still active
	setAutoThreshold("Otsu dark stack");	//use the stack histogram
	getThreshold(lower, upper);
	//print(lower);

	run("Z Project...", "projection=[Max Intensity]");
	setThreshold(lower, upper);
	run("Create Selection");
	run("Clear Outside", "stack");	//Threshold will be done in a z-cylinder
	roiManager("Add");
	run("Restore Selection");
}


//Returns the sum of all elements of an arrays, neglecting NaNs
function sumArray(array) {
	sum=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
	}
	return sum;
}


//Replace NaNs by zeros in the active image
function removeNaNs(image) {
	selectWindow(image);
	for(c=1;c<=channels;c++) {
		Stack.setChannel(c);
		showStatus("Removing NaNs...");
		showProgress(c/channels);
		for(z=1;z<=slices;z++) {
			Stack.setSlice(z);
			for(f=1;f<=frames;f++) {
				Stack.setFrame(f);
				changeValues(NaN,NaN,0);
			}
		}
	}
}


//Divides the elements of two arrays and returns the new array
function divideArrays(array1, array2) {
	divArray=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		divArray[a]=array1[a]/array2[a];
	}
	return divArray;
}


//Divides all elements of an array by a scalar
function divideArraybyScalar(array, scalar) {
	divided_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		divided_array[a]=array[a]/scalar;
	}
	return divided_array;
}


//Manual editing of detected nuclei
function edit_ROIs(image1) {
	shift=1;
	ctrl=2; 
	rightButton=4;
	alt=8;
	leftButton=16;
	insideROI = 32;
	
	flags=-1;
	//x2=-1; y2=-1; z2=-1; flags2=-1;
	
	selectWindow(image1);
	roiManager("Show All without labels");
	setOption("DisablePopupMenu", true);
	setBatchMode(true);
	resetMinAndMax();
	setBatchMode("show");
	color_ROIs();
	print("\\Clear");
	print("Delete, combine and draw new ROIs. \n- Left clicking while holding CTRL deletes a ROI.\n- Select multiple ROIs with shift-left mouse button and right-click to merge them into one ROI. \n- Draw new ROIs using the magic wand or the freehand tool and press 't' to add. \n- Press space bar when finished editing.\n");
	showMessage("Delete, combine and draw new ROIs. \n- Left clicking while holding CTRL deletes a ROI.\n- Select multiple ROIs with shift-left mouse button and right-click to merge them into one ROI. \n- Draw new ROIs using the magic wand or the freehand tool and press 't' to add. \n- Press space bar when finished editing.\n\nThis information is also printed in the log window.");
	print("\nStarting editing "+roiManager("count")+" ROIs...");
	
	setTool("freehand");
	roiManager("Show All Without Labels");
	setOption("DisablePopupMenu", true);
	
	nROIs = roiManager("Count");
	
	while(!isKeyDown("space")) {		//exit by pressing space bar
		getCursorLoc(x, y, z, flags);
		if(flags==17 || flags==18)	{	//(de)select multiple ROIs with shift-leftclick; delete ROI with rightclick
			for(i=0;i<roiManager("Count");i++) {
				roiManager("Select",i);
				if(Roi.contains(x, y)==true) {
				selected = Roi.getProperty("selected");
					//click to select a single ROI
					if(flags==17 && selected==false) {		//select ROI
						//print("selecting ROI "+i);
						Roi.setStrokeColor("red");
						Roi.setProperty("selected",true);
					}
					else if(flags==17 && selected==true) {	//deselect ROI
						//print("deselecting ROI "+i);
						Roi.setStrokeColor("cyan");
						//Roi.setFillColor("1900ffff");
						Roi.setProperty("selected",false);
					}
					else if(flags==18) {	//delete ROI
						roiManager("Delete");
						for(j=0;j<roiManager("Count");j++) {	//deselect all ROIs and rename
							roiManager("Select",j);
							roiManager("Rename", "ROI "+j);
						}
					}
				}
			}
			roiManager("Deselect");
			run("Select None");
			updateDisplay();
		}
	
		if(flags==4) {	//right button: combine selected ROIs
			selected_ROI_array = newArray(roiManager("Count"));	//create array with indices of selected ROIs
			j=0;
			for(i=0;i<roiManager("Count");i++) {
				roiManager("select",i);
				selected = Roi.getProperty("selected");
				if(selected==true) {
					selected_ROI_array[j] = i;
					j++;
					//print(j);
				}
			}
			//check if more than one ROI is selected. If yes, combine the selected ROIs and update the list
			selected_ROI_array = Array.trim(selected_ROI_array,j);
			//print(selected_ROI_array.length + " ROIs selected");
			if(selected_ROI_array.length > 1) {
				roiManager("Select",selected_ROI_array);
				roiManager("Combine");
				roiManager("Update");
				to_delete_array = Array.copy(selected_ROI_array);								//selecting and deleting redundant ROIs
				to_delete_array = Array.slice(selected_ROI_array,1,selected_ROI_array.length);	//create array without the first element
				roiManager("Deselect");
				roiManager("select", to_delete_array);
				roiManager("Delete");
				roiManager("Select",selected_ROI_array[0]);
				run("Enlarge...", "enlarge=1 pixel");			//remove wall between ROIs by enlarging and shrinking with 1 pixel
				run("Enlarge...", "enlarge=-1 pixel");
				roiManager("Update");
				
				setKeyDown("none");
				
				color_ROIs();
			}
		}
	
	
		if(nROIs!=roiManager("Count")) {	//change in the number of ROIs 
			run("Select None");
			color_ROIs();
			nROIs = roiManager("Count");
		}
	
		else wait(50);
	}	//end of while loop
	
	//Deselect and rename all ROIs once more
	color_ROIs();
}


function color_ROIs() {
	run("Remove Overlay");

	for(j=0;j<roiManager("Count");j++) {	//fill all ROIs
		roiManager("Select",j);
		roiManager("Rename", "ROI "+j+1);
		Roi.setProperty("selected",false);
		//Roi.setFillColor("1900ffff");	//10% cyan fill
	}
	roiManager("Deselect");
	if(roiManager("count")>0) run("From ROI Manager");	//Add overlay containing the ROI fill
	roiManager("Select All");
//	roiManager("Set Color", "cyan");
	roiManager("Deselect");
	roiManager("Show All");
	updateDisplay();
}
