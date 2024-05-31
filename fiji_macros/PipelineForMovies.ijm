_RootFolder = getDirectory("Choose a Directory"); 


Dialog.create("Settings");
  Dialog.addNumber("Channels:", 3);
  Dialog.addNumber("Slices:", 41);
  Dialog.addNumber("Timepoints", 18);
  Dialog.show();
CHANNELS =Dialog.getNumber();
_slices = Dialog.getNumber();
_frames = Dialog.getNumber();

lookUpTables = newArray("Grays", "Red", "Green");
XY_DIMENSIONS = 500;


fileList = getTifFileList(_RootFolder);

for (i=0; i<fileList.length; i++){
	print(fileList[i]);
	filename = substring(fileList[i], 0, lengthOf(fileList[i])-4);
	print(filename);
	File.makeDirectory(_RootFolder + filename + File.separator);
	open(_RootFolder + fileList[i]);
	setBatchMode(true);
	if (_frames > 1) {
		run("Properties...", "channels=CHANNELS slices=_slices frames=_frames unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
		run("Arrange Channels...", "new=123");
		}
	
	else {
		//creates a multichanel image
		// WARNING: This is only optimized for pictures. As soon as I have a video, this will be adapted for those. Check Script 1 for if else
		run("Properties...", "channels=CHANNELS slices=_slices frames=_frames unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
		run("Arrange Channels...", "new=123");}
	
	run("Bio-Formats Exporter", "save=" + _RootFolder + filename + File.separator + filename +".tif write_each_timepoint compression=Uncompressed");
	close();
	setBatchMode(false);
   }
   

allFiles = getFileList(_RootFolder);
directories = getDirectories(allFiles);


setBatchMode(true);

getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
begin = (60*hour) + minute; 
// the next three lines generate tifs which contain three additional slices
// at the beginning and three slices at the end, where the signal fades out.
// this is neccessary for Mitograph. The images are saved under the 'Added_Slices'
// folder. Also, a ProjectionSummary.tif is generated that contains Z-Projections
// of all images.

for (directory = 0; directory < directories.length; directory++){ 
	tifFiles = getTifFileList(directories[directory]);
	addThreeSlices(CHANNELS, directories[directory], tifFiles);
}

setBatchMode(false);

end =(60*hour +  minute) - begin;
print("It took me " + end + " minutes to finish the job");

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

function addThreeSlices(CHANNELS, folder, tifFiles){
	for (file = 0; file < tifFiles.length; file++){
		
		fileName = tifFiles[file];
		open(folder + fileName);
		image = getImageID();
		for (channel=1; channel<CHANNELS+1; channel++){
			setSlice(channel);
			run(lookUpTables[channel-1]);
		}
		selectImage(image);
		run("Split Channels");
		channelNames = newArray("C1-" + fileName, "C2-" + fileName, "C3-" + fileName, "C4-" + fileName, "C5-" + fileName, "C6-" + fileName);
		for (channel=1; channel<CHANNELS+1; channel++){ 
			selectWindow(channelNames[channel-1]);
			for (i = 0; i < 3; i++){
				setSlice(nSlices);
				run("Copy");
				run("Add Slice", "add=slice");
				run("Paste");
				setSlice(1);
				run("Copy");
				run("Add Slice", "add=slice prepend");
				run("Paste");}
			for (slice=1; slice<4; slice++){
				if (slice==1){
					setSlice(1);
					run("Multiply...", "value=0.05 slice");
					setSlice(nSlices);
					run("Multiply...", "value=0.05 slice");}
				if (slice==2){
					setSlice(2);
					run("Multiply...", "value=0.4 slice");
					setSlice(nSlices-1);
					run("Multiply...", "value=0.4 slice");}
				if (slice==3){
					setSlice(3);
					run("Multiply...", "value=0.8 slice");
					setSlice(nSlices-2);
					run("Multiply...", "value=0.8 slice");}
			}
			
		}
		run("Merge Channels...", "c1=" + channelNames[0] + " c2=" + channelNames[1] + " c3=" + channelNames[2] + " create");
		
		for (channel=1; channel<CHANNELS+1; channel++){ 
			Stack.setChannel(channel);
			run("Clear Results");
			run("Measure Stack...", "slices order=czt(default)");
			run("Summarize");
			max = getResult("Max");
			setMinAndMax(0, max);
		}
		save(folder + fileName);
		close();
	}
}
