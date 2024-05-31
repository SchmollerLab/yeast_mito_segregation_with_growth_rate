_RootFolder = getDirectory("Choose a Directory");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
begin = (60*hour) + minute; 
setBatchMode(true);
channels = getNumber("How many wavelengths?", 4);
// channel names are given to the images in the macro, which is needed to merge them later
// if you have more than 6 channels you will need extend the macro
channel_names = newArray("a", "b", "c", "d", "e", "f");
File.makeDirectory(_RootFolder + "Merged/");
pre = "";
post = ".ome_cmle_";


//***************************************************
//************  I M P O R T A N T  ******************
//***************************************************
// the script will only work if all files have the same number of channels
// if you have multiple files with different amounts of channesl separate them first into different folders and 
// run the script on each folder

// the following two lines get an Array of unique filenames from the chosen folder
// these commands only get the file stem and remove the string added by huygens, e.g. ".ome_cmle_ch001.tif"
files = getTifFileList();
files = ArrayUnique(files);

for (i=0; i<files.length; i++) {
	for (channel=0; channel<channels; channel++){
		open(_RootFolder + files[i] + post + "ch0" + channel + ".tif");
		rename(channel_names[channel]);}
	//print("c1=" + imageb + " c2=" + imagec + " c3=" + imaged +" c4=" + imagea + " create");
	if(channels==2){
		run("Merge Channels...", "c1=a c2=b create");}
	if(channels==3){
		run("Merge Channels...", "c1=a c2=b c3=c create");}
	if(channels==4){
		run("Merge Channels...", "c1=a c2=b c3=c c4=d create");}
	if(channels==5){
		run("Merge Channels...", "c1=a c2=b c3=c c4=d c5=e create");}
	if(channels==6){
		run("Merge Channels...", "c1=a c2=b c3=c c4=d c5=e c6=f create");}
		
	getDimensions(w, h, channels, slices, frames);
	Stack.setSlice(slices/2);
	for (j=1; j<=channels; j++) {
		Stack.setChannel(j);
		run("Clear Results");
		run("Select All");
		run("Measure");
		min = getResult("Min", 0);
		max = getResult("Max", 0);
		print(min);
		print(max);
		setMinAndMax(min, max);	
		}
	saveAs("Tiff", _RootFolder + "Merged/" + files[i] + "_decon.tif");
	close();


	
}

setBatchMode(false);
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
end =(60*hour +  minute) - begin;
print("It took me " + end + " minutes to finish the job");

function getTifFileList(){
	item = 0;
	fileList = getFileList(_RootFolder);
	tifFileList = newArray();
	while (item < fileList.length)  {
		if (endsWith(fileList[item],".tif") ) {
			print(fileList[item]);
			new_file = substring(fileList[item], 0, lengthOf(fileList[item])-18);
			print(new_file);
			tifFileList = Array.concat(tifFileList, new_file);
			}
		item += 1;
	}
	return tifFileList;
}

function ArrayUnique(array) {
	array 	= Array.sort(array);
	array 	= Array.concat(array, 999999);
	uniqueA = newArray();
	i = 0;	
   	while (i<(array.length)-1) {
		if (array[i] == array[(i)+1]) {
			//print("found: "+array[i]);			
		} else {
			uniqueA = Array.concat(uniqueA, array[i]);
		}
   		i++;
   	}
	return uniqueA;
}