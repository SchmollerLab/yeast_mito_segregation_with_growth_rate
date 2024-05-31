var _RootFolder = getDirectory("Choose a Directory");
File.makeDirectory(_RootFolder + "csv_files");
var _csv = _RootFolder + "/csv_files/";
tifFiles = getTifFileList();

for (file=0; file<tifFiles.length; file++){
	run("Clear Results");
	open(_RootFolder + "/" + tifFiles[file]);
	current = getImageID();
	run("Smooth");
	setThreshold(20116, 65535, "raw");
	run("Convert to Mask");
	makeOval(456, 237, 1647, 1671);
	waitForUser("Select Plate Area");
	run("Dilate");
	run("Watershed");
	run("Analyze Particles...", "size=0.0001-10000 circularity=0.7-1.00 display exclude clear add");
	saveAs("Results", _csv + tifFiles[file] + ".csv");
	selectImage(current);
	close();
	
	}
	


function getTifFileList(){
	item = 0;
	fileList = getFileList(_RootFolder);
	tifFileList = newArray();
	while (item < fileList.length)  {
		if (endsWith(fileList[item],".Tif") ) {
			print(fileList[item]);
		
			tifFileList = Array.concat(tifFileList, fileList[item]);
			}
	item += 1;
	}
	return tifFileList;
}