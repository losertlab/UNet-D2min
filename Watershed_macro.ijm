dir1 = getDirectory();
format = "TIFF";
dir2 = getDirectory();
list = getFileList(dir1);
setBatchMode(true);
for  (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	open(dir1+list[i]);
	//run("Brightness/Contrast...");
	setMinAndMax(0, 1);
	run("Apply LUT");
	run("Convert to Mask");
	run("Watershed");
	saveAs("Tiff", dir2+list[i]);
	close();
}