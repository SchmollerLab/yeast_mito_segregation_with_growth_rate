/*
 * Finds ROIs of budding cells  with YeastMate and measures intensity in fluorescent channel of MaxProjection (Channel_2). 
 * Saves intensity values in csv file.
   IMPORTANT: YeastMate + FIJI plugin must be installed and backends activated in order for the script to work 
 * 
 */

_RootFolder = getDirectory("Choose a Directory");
run("Set Measurements...", "area mean integrated redirect=None decimal=3");

allFiles = getFileList(_RootFolder);
directories = getDirectories(allFiles);
print(directories.length);
setBatchMode(true);
for (directory = 0; directory < directories.length; directory++){
	print(directories[directory]);
	thisfile = directories[directory];	
	print(thisfile);
	tifFiles = getTifFileList(directories[directory]);
	run("Clear Results");
	for (i = 0; i < tifFiles.length; i++) {
		open(directories[directory] + tifFiles[i]);
		original = getTitle();
		cell = substring(original, (lengthOf(original)-29), lengthOf(original)-6);
		print(cell);
		cell = replace(cell, "_", "");
		run("Z Project...", "projection=[Max Intensity]");
		zProjMito=getTitle();
		selectWindow(original);
		currentimage= getImageID();
		Stack.setPosition(1, 10,0);
		close("ROI Manager");
		run("YeastMate", "scorethresholdsingle=0.9400000000000001 scorethresholdmating=0.75 scorethresholdbudding=0.77 minnormalizationqualtile=0.015 maxnormalizationqualtile=0.985 addsinglerois=false addmatingrois=false addbuddingrois=true showsegmentation=false onlyselectedclassesinmask=false processeveryframe=false mintrackingoverlap=0.25 ipadress=127.0.0.1:11005");
		wait(5);
		selectWindow(zProjMito);
		for (roi = 0; roi < roiManager("count"); roi++) {
			rName = RoiManager.getName(roi);
		//print(rName);
			if (endsWith(rName, "budding")) {
				print(rName);
			print("ROI: " + roi);
			roiManager("Select",roi);
			Stack.setPosition(2, 10,0);
			run("Measure");
			//run("Convex Hull");
			}}
		selectWindow(original);
		close();
		selectWindow(zProjMito);
		close();
		roiManager("reset");
		close("ROI Manager");
	
		run("Close All");
		

	}
	saveAs("Results", _RootFolder + cell + "_Results.csv");
	close("ROI Manager");

}
setBatchMode(false);
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