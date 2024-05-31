

// Matheus Viana - vianamp@gmail.com - 7.29.2013
// ==========================================================================

// This macro must be used to generate a stack of max projection from a set
// of microscope frames. In general, the max projection stack is then used
// to drawn ROIs around the cells that are going to be further analysed. 

// Selecting the folder that contains the TIFF frame files

_RootFolder = getDirectory("Choose a directory");
allFiles = getFileList(_RootFolder);
directories = getDirectories(allFiles);
setBatchMode(true);
for (directory = 0; directory < directories.length; directory++){
	thisfile = directories[directory];
	//thisfile = replace(thisfile, "/", "");
	print(thisfile+  "Hyperstack");
	item = 0;
	ntiff = 0;
	_List = getFileList(thisfile + "Hyperstack/");
	while (item < _List.length)  {
		if ( endsWith(_List[item],".tif") ) {
			if (ntiff==0) {
				open(thisfile + "Hyperstack/" + _List[item]);
				w = getWidth();
				h = getHeight();
				close();
			}
			ntiff++;
		}
		item++;
	}
	if (ntiff== 0) {
		showMessage("No TIFF files were found.");
	} else {
		print("Number of TIFF files: " + ntiff);
	}
	
	// Generating the max projection stack
	
	newImage("MaxProjs", "16-bit black", w, h, ntiff);
	
	item = 0; im = 0;
	while (item < _List.length)  {
		if ( endsWith(_List[item],".tif") ) {
			im++;
			open(thisfile + "Hyperstack/" + _List[item]);
			_FileName = split(_List[item],"."); 
			_FileName = _FileName[0];
			print(_FileName);
			// the graychannel needs to be removed for the maximum projection
			run("Split Channels");
			selectWindow("C1-" + _FileName +  ".tif");
			close();
			run("Merge Channels...", "c1=C2-" + _FileName  + ".tif c2=C3-" + _FileName + ".tif create");
			run("Z Project...", "start=1 stop=500 projection=[Max Intensity]");
			run("Copy");
			close();

			
			selectWindow("MaxProjs");
			setSlice(im);
			run("Paste");
			setMetadata("Label",_FileName);
		}
		item++;
	}
	
	// Saving max projection stack
	
	run("Save", "save=" +  thisfile + "Hyperstack/" + "MaxProjs.tif");
	run("Close All");

}
setBatchMode(false);
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