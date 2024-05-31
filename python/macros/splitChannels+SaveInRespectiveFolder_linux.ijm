// This macro goes through folders in directory, splits each picture into the
// 3 different channels (gray, green, red) and saves them into new folders, respectively
//The output of the PipelineForMovies macro can be used as input for this macro



_RootFolder = getDirectory("Choose a directory"); 
print(_RootFolder);
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
begin = (60*hour) + minute; 

channels = 3;
allFiles = getFileList(_RootFolder);
directories = getDirectories(allFiles);

setBatchMode(true);
//goes to any folder in the _RootFolder directory and creates 3 new folders for each channel
for (directory = 0; directory < directories.length; directory++){
	print(directories[directory]);
	File.makeDirectory(directories[directory] + "Channel0" );
	File.makeDirectory(directories[directory] + "Channel1" );
	File.makeDirectory(directories[directory] + "Channel2" );
	list = getTifFileList(directories[directory]);
	print(list.length);
	//now it opens all the pictures in the folder and splits the channels
	for(i=0; i < list.length; i++){
		filename= directories[directory] + list[i];
		print("Filename: " + filename);
		run("Bio-Formats", "open=" + filename + " autoscale color_mode=Default rois_import=[ROI manager] split_channels view=[Standard ImageJ] stack_order=Default");
		// every channel will be saved in its respective folder (channel 0,1,2)
		for(j=0; j < channels; j++){
			pic = getTitle();
			print("Pic: " + pic);
			selectWindow(pic);
			pic = list[i] + "_"+ substring(pic, (lengthOf(pic)-3),lengthOf(pic));
			pic = replace(pic, ".tif", "");
			print("Pic final: " + pic);
			Channelnumber = substring(pic, (lengthOf(pic)-1), lengthOf(pic));
			final_destination = directories[directory];
			//final_destination = replace(final_destination, "/", "");
			saveAs("TIFF",final_destination + "Channel" + Channelnumber  + File.separator + pic);
			close();	
		}
		
	}

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