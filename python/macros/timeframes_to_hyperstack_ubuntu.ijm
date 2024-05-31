/// merges single timeframes into a hyperstack. rootfolder: the folder where every series is in its own folder.
/// the macro then runs through every folder, merges the timeframes and saves them in a new folder 
/// called "Hyperstack"

_RootFolder = getDirectory("Choose a Directory");
allFiles = getFileList(_RootFolder);
directories = getDirectories(allFiles);

Dialog.create("Settings");
  Dialog.addNumber("Channels:", 2);
  Dialog.addNumber("Slices:", 21);
  Dialog.addNumber("Timepoints", 57);
  Dialog.show();
_channels =Dialog.getNumber();
_slices = Dialog.getNumber();
_frames = Dialog.getNumber();



getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
begin = (60*hour) + minute; 


setBatchMode(true);
for (directory = 0; directory < directories.length; directory++){

	thisfile = directories[directory];
	print(thisfile);
	//thisfile = replace(thisfile, "/", "");
	File.makeDirectory(thisfile +  "Hyperstack");
	_List = getTifFileList(thisfile);
	StackName = _List[0];
	StackName= replace(StackName, "_T0_decon", "");
	print(StackName);
	run("Image Sequence...", "select="+thisfile +"dir="+thisfile+" sort");
	
	run("Stack to Hyperstack...", "order=xyczt(default) channels=_channels slices=_slices frames=_frames display=Grayscale");
	saveAs("Tiff", thisfile + File.separator + "Hyperstack" + File.separator +StackName);
	close();
	}



getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
end =(60*hour +  minute) - begin;
print("It took me " + end + " minutes to finish the job");

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