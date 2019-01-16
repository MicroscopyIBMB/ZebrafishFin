// FinSegmentation. January 11th 2019.
/*Authors:  Hanh Nguyen*, Jaume Boix-Fabrés# , Nadine Peyriéras* and Elena Kardash*
* BioEmergences Laboratory (USR 3695), CNRS, University Paris-Saclay, 91190, Gif-sur-Yvette,
France
# Molecular Imaging Platform, Molecular Biology Institute of Barcelona (IBMB), Spanish National
Research Council (CSIC); 08028, Barcelona, Spain.

Fiji lifeline 22 Dec 2015

This macro obtains the AP and DV axis length and sagital projections from a fluorescence 
microscopy XYZt image of a Zebrafish pectoral fin. Segmentation is performed obtaining
the average Z-projection and automatic thresholding frame by frame. Fin is selected by size 
filtering in particle analysis. Then, it calculates the AP and DV axis by estimation of 
the Feret's diameters of the segmented object. Also it projects the sagital views of both 
AP and DV axis. Binary masks obtained after segmentatin are saved as a verification control
*/


// DIALOG FOR MACRO PARAMETERS
dir=getDirectory("Choose a directory to save the final images:");
Dialog.create("Segmentation and Scaling Parameters:");
Dialog.addCheckbox("Flip Vertically",false);
Dialog.show;
flip=Dialog.getCheckbox();


// PRE-PROCESSING STEPS  (substract background, Gaussian filter, remove scale)
id=getImageID();  
getVoxelSize(widthV,heightV,depth,unit);
run("Duplicate...", "title=fin duplicate");
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
getDimensions(width,height,channels,slices,frames);


// ARRAYS TO STORE FERET DIAMETERS
APaxis=newArray(frames);
DVaxis=newArray(frames);
XC=newArray(frames);
YC=newArray(frames);
Angle=newArray(frames);


// PREPROCESSING STEPTS

selectWindow("fin");
run("Z Project...", "projection=[Sum Slices] all");
rename("Zproj");
run("Subtract Background...", "rolling=100 stack");
run("Median...", "radius=2 stack");

// SECOND PART: SEGMENT THE PROJECTION, FIND THE FERET AXIS, DRAW IT AND RESLICE 

setBatchMode(true);

for(i=1; i<frames+1; i++) {

// SELECT A FRAME, THRESHOLD IT AND ADJUST THE OBJECT TO AN ELLIPSE
selectWindow("Zproj");
setSlice(i);
run("Duplicate...","title=mask");
setAutoThreshold("Default dark");
run("Convert to Mask");

run("Set Measurements...", "centroid feret's redirect=None decimal=3");
run("Analyze Particles...", "size=1000-Infinity pixel clear add");
roiManager("select",0);
run("Fit Ellipse");
roiManager("Add");
roiManager("Select", 0);
roiManager("Delete");
roiManager("Set Color", "red");
roiManager("Set Line Width", 2);
run("RGB Color");
//roiManager("Select",0);
run("Flatten");
rename("check"+i);

// MEASURE FERET AXIS, FERET ANGLE AND CENTROID POSITION. 

roiManager("Select", 0);
roiManager("Measure");
APaxis[i-1]=getResult("Feret");
DVaxis[i-1]=getResult("MinFeret");

XC[i-1]=getResult("X");
YC[i-1]=getResult("Y");
Angle[i-1]=getResult("FeretAngle")*PI/180;

selectWindow("mask");
run("Close");
roiManager("reset");
}

// CLOSE DUPLICATED 4D STACK 
selectWindow("fin");
run("Close");

// GET MAXIMUM FERET DIAMETERS AND SPECIFY LENGTH OF LINES FOR RESLICE (XZ, YZ)
Array.getStatistics(APaxis,minAP,maxAP);
Array.getStatistics(DVaxis,minDV,maxDV);
APline=round(maxAP+(width-maxAP)/2);
DVline=round(maxDV+(width-maxDV)/2);

run("Images to Stack", "name=Segmentation_Check title=check use");
saveAs("tif",dir+"Segmentation_Check.tif");


// DRAW ARROWS IN THE Z-PROJECTED TIME-SERIES
selectWindow("Zproj");
run("RGB Color");
setForegroundColor(255, 255, 0);
for (i=1;i<frames+1;i++){
	setSlice(i);
	Arrow(XC[i-1],YC[i-1],APline,Angle[i-1]);
	Arrow(XC[i-1],YC[i-1],DVline,(Angle[i-1]+PI/2));
}
saveAs("tif",dir+"Zproj_series.tif");
run("Close");

// RESLICE AP AXIS
for (i=1;i<frames+1;i++) {
	selectImage(id);
	run("Make Substack...", "slices=1-"+slices+" frames="+i+"");
	rename("temp");
	Reslice(XC[i-1],YC[i-1],APline,Angle[i-1],depth);
	rename("XZ"+i);
	selectWindow("temp");
	run("Close");
}

// COMBINE ALL INDIVIDUAL XZ FRAMES IN A TIME SERIES AND SAVE TO A FILE
selectWindow("XZ1");
rename("XZcomb");
for (i=1;i<frames;i++) {
	run("Concatenate...", "  title=XZcomb image1=XZcomb image2=XZ"+(i+1)+" image3=[-- None --]");
}

selectWindow("XZcomb");
if (flip==true) {run("Flip Vertically");}
saveAs("tif",dir+"AP.tif");

// RESLICE DV AXIS
for (i=1;i<frames+1;i++) {
	selectImage(id);
	run("Make Substack...", "slices=1-"+slices+" frames="+i+"");
	rename("temp");
	Reslice(XC[i-1],YC[i-1],DVline,(Angle[i-1]+PI/2),depth);
	rename("YZ"+i);
	selectWindow("temp");
	run("Close");
}

// COMBINE ALL INDIVIDUAL YZ FRAMES IN A TIME SERIES AND SAVE TO A FILE
selectWindow("YZ1");
rename("YZcomb");
for (i=1;i<frames;i++) {
	run("Concatenate...", "  title=YZcomb image1=YZcomb image2=YZ"+(i+1)+" image3=[-- None --]");
}

selectWindow("YZcomb");
if (flip==true) {run("Flip Vertically");}
saveAs("tif",dir+"DV.tif");
setBatchMode(false);

// CREATE TABLE FOR FERET DIAMETERS AND SAVE IT

run("Table...", "name=[Ferets] width=400 height=300 menu");
print("[Ferets]", "\\Headings:"+"Frame \t AP axis \t DV axis");
for(i=0; i<frames; i++){
	print("[Ferets]", ""+(i+1)+ "\t" + APaxis[i]*heightV + "\t" + DVaxis[i]*heightV);
}

//SAVE RESULTS TABLE AS .XLS IN RESULTS FOLDER
selectWindow("Ferets");
saveAs("Text", dir+"Ferets.xls");
run("Close");



// DEFINE FUNCTION TO DRAW ARROWS ON THE IMAGE
function Arrow(CenterX,CenterY,Length,Angle) {
	X1=CenterX-(Length/2)*cos(Angle);
	Y1=CenterY-(Length/2)*sin(-Angle);
	X2=CenterX+(Length/2)*cos(Angle);
	Y2=CenterY+(Length/2)*sin(-Angle);
	makeArrow(X1, Y1, X2, Y2, "filled");
	Roi.setStrokeWidth(2);
	Roi.setStrokeColor("yellow");
	run("Draw","slice");
}

function Reslice(CenterX,CenterY,Length,Angle,Zdepth) {
	X1=CenterX-(Length/2)*cos(Angle);
	Y1=CenterY-(Length/2)*sin(-Angle);
	X2=CenterX+(Length/2)*cos(Angle);
	Y2=CenterY+(Length/2)*sin(-Angle);
	makeLine(X2,Y2,X1,Y1);
	run("Reslice [/]...", "output="+d2s(Zdepth,3)+" slice_count=1");
}
	
