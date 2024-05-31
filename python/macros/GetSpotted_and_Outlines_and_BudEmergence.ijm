//// runs 3D maxima finder on every timeframe in Channel2 (mtDNA) and saves the results in the SpotCoordinates folder,
///which was created in Channel 1
/// then you have to select the outline of each cell (+daughter cell) by clicking with the point tool at the edges of the cell (hold SHIFT)
/// After that it creates hyperstack of graychannel timeframes, from which you can now select the timeframe where the bud starts emerging.
/// saves both csv files in new folders in channel1


_RootFolder = getDirectory("Choose a Directory");
allFiles = getFileList(_RootFolder);
directories = getDirectories(allFiles);

// 3D maxima finder 
/*
setBatchMode(true);
for (directory = 0; directory < directories.length; directory++){
	thisfile = directories[directory];
	//thisfile = replace(thisfile, "/", "");
	print(thisfile+ File.separator + "Channel1");
	File.makeDirectory(thisfile+ File.separator + "Channel1" + File.separator + "Cell_outline");
	File.makeDirectory(thisfile+ File.separator + "Channel1" + File.separator + "SpotCoordinates");
	File.makeDirectory(thisfile+ File.separator + "Channel1" + File.separator + "Bud_Emergence");
	Channel2 = directories[directory] + File.separator + "Channel2";
	//Channel2 = replace(Channel2, "/","");
	print("Nur channel: " + Channel2);
	print("beides: " + _RootFolder + Channel2);
	list = getTifFileList(Channel2);
	for(i=0; i < list.length; i++){
		filename= Channel2 + File.separator + list[i];
		Imagename = list[i] + ".csv";
		Imagename = replace(Imagename, ".tif","");
		print("filename: " + filename);
		print("Imagename: " + Imagename);
	   	open(filename);
	   	run("3D Maxima Finder", "minimmum=0 radiusxy=1.50 radiusz=1.50 noise=500");
		saveAs("Results", thisfile+ File.separator + "Channel1" + File.separator + "SpotCoordinates"+ File.separator  + Imagename);
		close();
	}
}
run("Close All");
setBatchMode(false);
print("hallooo");

//creates hyperstack of graychannel timeframes, from which you can now select the timeframe where the bud starts emerging
for (directory = 0; directory < directories.length; directory++){
	thisfile = directories[directory];
	//thisfile = replace(thisfile, "/", "");
	Channel0 = directories[directory]  + "Channel0/";
	//Channel0 = replace(Channel0, "/","");
	run("Clear Results");
	print("\\Clear");
	run("Image Sequence...", "select="+Channel0 +"dir="+Channel0+" sort");
	run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=47 frames=31 display=Grayscale");
	run("Z Project...", "projection=[Max Intensity] all");
	zProject = getImageID;
	getDimensions(width, height, channels, slices, frames);
	//print(frames);
	for (frame = 0; frame < frames; frame++) {
		setSlice(frame+1);
		run("Set Measurements...", "mean redirect=None decimal=3");
		run("Measure");
		//print("slice no" + frame);
	}
	means = newArray(nResults);
	for(i=0; i<nResults; i++) {
		means[i] = getResult("Mean", i);
	}
	//Array.print(means);
	Array.getStatistics(means, min, max, mean, stdDev);
	//print("MAX: "+mean);
	for (frame = 0; frame < frames; frame++) {
		value =  mean/means[frame];
		//print(value);
		setSlice(frame+1);
		run("Multiply...", "value=" + value + " slice");
		run("Measure");
		run("Clear Results");
		//print("slice no" + frame);
	}
	//run("Enhance Contrast", "saturated=0.35");
	run("Clear Results");
	setSlice(1);
	waitForUser("Find the point of first bud emergence!");
	result = getNumber("which timeframe?", 1 );
	print(result);
	waitForUser("Find the point of second bud emergence!");
	result = getNumber("which timeframe?", 1 );
	print(result);
	selectWindow("Log");
	saveAs("Text", thisfile+ File.separator + "Channel1"+ File.separator + "Bud_Emergence"+File.separator + "BudEmergence.csv");
	run("Close All");

}

*/

// Get outlines

item = 0;
for (directory = 0; directory < directories.length; directory++){
	thisfile = directories[directory];
	//thisfile = replace(thisfile, "/", "");
	print(thisfile+ File.separator + "Channel1");
	Channel0 = directories[directory] + "Channel0";
	//Channel0 = replace(Channel0, "/","");
	print(Channel0);
	_List = getTifFileList(Channel0);
	run("Clear Results");
	for(i=0; i < _List.length; i++){
			open(Channel0 + File.separator + _List[i]);
			run("Set Scale...", "distance=0 known=0 unit=pixel");
			original = getImageID;
			run("Z Project...", "projection=[Max Intensity]");
			zProject = getImageID;
			run("Enhance Contrast", "saturated=0.35");
			run("In [+]");
			run("In [+]");
			//run("In [+]");
			setTool("point");
			waitForUser("Please select points to describe outline of mother cell");
			run("Set Measurements...", "  redirect=None decimal=3");
			run("Measure");
			print(Channel0 + File.separator + replace(_List[i], ".tif", ".csv"));
			saveAs("Results", thisfile+ File.separator + "Channel1"+ File.separator +"Cell_outline"+ File.separator + replace(_List[i], ".tif", ".csv"));
			close();
			run("Clear Results");
			selectImage(original);
			close();
			
	}
}


//goes through BF folder, duplicates each image containing only a single plane 
// (somewhere in the middle of the stack, to get the cell in the focus

item = 0;
for (directory = 0; directory < directories.length; directory++){
	thisfile = directories[directory];
	thisfile = replace(thisfile, "/", "");
	print(thisfile+ File.separator + "Channel1");
	Channel0 = directories[directory] + File.separator + "Channel0";
	Channel0 = replace(Channel0, "/","");
	print(Channel0);
	_List = getTifFileList(Channel0);
	run("Clear Results");
	for(i=0; i < _List.length; i++){
			open(Channel0 + File.separator + _List[i]);
			original = getImageID;
			run("Duplicate...", "duplicate range=22-22 use");
			single= getImageID();
			print("SIngle:" +single);
			print(Channel0 + File.separator + _List[i]);
			selectImage(single);
			saveAs("Tiff", thisfile+ File.separator + "Channel1"+ File.separator +"BF_Single_plane"+ File.separator + _List[i]);
			close();
			selectImage(original);
			close();
	}
}



// goes through folder, selects only TIF-files and stors them in a new array called tifFileList
function getTifFileList(directory){
	item = 0;
	fileList = getFileList(directory);
	tifFileList = newArray();
	while (item < fileList.length)  {
		if (endsWith(fileList[item],".tif") ) {		
			tifFileList = Array.concat(tifFileList, fileList[item]);
			}
	item += 1;
	}
	return tifFileList;
}


function getDirectories(files){
	directoryList = newArray();
	for (file = 0; file < files.length; file++){
		if (File.isDirectory(_RootFolder + files[file])){
			directoryList = Array.concat(directoryList, _RootFolder + files[file]);
		}
	}
	return directoryList;
}
