
/*
Advanced Optical Microscopy Unit
Scientific and Technological Centers. Cl�nic Campus
University of Barcelona
C/ Casanova 143
Barcelona 08036 
Tel: 34 934037159
Fax: 34 934024484
mail: confomed@ccit.ub.edu
------------------------------------------------
Gemma Martin, Maria Calvo.
------------------------------------------------

Name of Macro: 143_Cell_Intensity_Time_Lapse.ijm

Date: 12/09/2023

version: 11 
															
Changes from last version: measurements without membrane

*/



if(isOpen("Results")){
    IJ.deleteRows(0, nResults);
}
run("ROI Manager...");
roiManager("reset"); //to delete previous ROIs
IJ.deleteRows(0, nResults);



dir = getDirectory("Choose images folder");
list=getFileList(dir);
dirRes=dir+"Results"+File.separator;
File.makeDirectory(dirRes);
run("ROI Manager...");

run("Set Measurements...", "area mean min integrated display redirect=None decimal=5");

run("Options...", "iterations=1 count=1 do=Nothing");



for(i=0;i<list.length;i++){
	
	if(endsWith(list[i],".czi")){

		run("Bio-Formats Importer", "open=[" + dir + list[i] + "] autoscale color_mode=Default open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

		numseries=nImages;
		
		run("Close All");

		
		for (k=1; k<=numseries; k++) {

			
			run("Bio-Formats Importer", "open=[" + dir + list[i] + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+k);

			title=getTitle();
			t=title;

			path = dirRes+t+"_Results.txt";
			
			File.append( "frame \t cell \t Area cell (um2) \t Red Intensity Cell \t points cell \t Area cell - membrane (um2) \t Red Intensity Cell - membrane  \t points cell - membrane", path);



			ImageID=getImageID();
			run("Duplicate...", "duplicate");

			Stack.setChannel(3);
			run("Blue");

			run("RGB Color", "frames keep");

			run("StackReg", "transformation=[Rigid Body]");

			run("Split Channels");



			
	        for (j=1; j<=nImages; j++) {

	        	selectImage(j);

	        	title=getTitle();
			
	       		t=substring(title,0,lengthOf(title)-7);


		        if (matches(title,".*blue.*")==1){
		        	bf=getImageID();
		        	run("Grays");
		        	rename("BF");
	
				}	
				
		        if (matches(title,".*red.*")==1){
					red=getImageID();
					run("Red");
					rename("Red");

				}

		        if (matches(title,".*green.*")==1){
					green=getImageID();
					run("Green");
					rename("Green");

				}

	        }


			// ********************************* Z project red ***********************run("Random");**************
	
			selectImage(green);

			slices=nSlices;


			IJ.deleteRows(0, nResults);

			run("Duplicate...", "duplicate");

			run("Median...", "radius=5 stack");
			
			rename("videogreen");


			selectImage(red);

			IJ.deleteRows(0, nResults);

			run("Duplicate...", "duplicate");

			run("Median...", "radius=3 stack");

			rename("videored");

			
			
			imageCalculator("Subtract create stack", "videogreen","videored");

			run("Z Project...", "projection=Median");
			run("Gamma...", "value=0.70");

			run("Find Maxima...", "prominence=10 output=[Segmented Particles]");

		
			run("Create Selection");
			roiManager("Add");

	
			roiManager("save", dirRes+"ROI_"+t +".zip")

			selectWindow("videogreen");
			run("Duplicate...", "duplicate");


			run("Convert to Mask", "method=Huang background=Dark calculate");
			run("Close-", "stack");
			run("Dilate", "stack");

			run("Analyze Particles...", "size=500-Infinity pixel show=Masks stack");

			imageCalculator("Min create stack", "Mask of videogreen-1","MED_Result of videogreen Segmented");
			
			// Check correct ROIs from the original movie

			run("Analyze Particles...", "size=500-Infinity pixel show=Masks stack");
			rename("maskcells");

			run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=10 redirect_to=none");
			run("3D Objects Counter", "threshold=128 slice=60 min.=10 max.=31457280 objects");

			selectWindow("Objects map of maskcells");
			rename("objects");

			run("Duplicate...", " ");

			run("Random");

			saveAs("Tiff", dirRes+"cellsegmentation_"+t +".tif");

			selectWindow("objects");

			IJ.deleteRows(0, nResults);
			run("Z Project...", "projection=[Max Intensity]");
			run("Measure");
			ncells=getResult("Max", 0);


			selectWindow("objects");
			run("Duplicate...", "duplicate");
			setOption("BlackBackground", false);
			run("Convert to Mask", "method=Percentile background=Dark calculate");
			rename("maskobjects");
			

			selectImage("Red");
			run("Duplicate...", "duplicate");
			run("Subtract Background...", "rolling=50 stack");
			rename("Red2");
			run("Duplicate...", "duplicate");
			rename("Red3");

			
			for (j=0;j<slices;j++) {

				IJ.deleteRows(0, nResults);
				roiManager("reset");

				
				areacell=newArray(ncells);
				intensityredcell=newArray(ncells);
				numpointscell=newArray(ncells);
				
				areamem=newArray(ncells);
				intensityredmem=newArray(ncells);
				numpointsmem=newArray(ncells);

				selectImage("Red3");
				setSlice(j+1);

				run("Find Maxima...", "prominence=40 strict output=[Single Points]");
				rename("Red4"+j);
				run("Properties...", "channels=1 slices=1 frames=1 pixel_width=1.0000000 pixel_height=1.0000000 voxel_depth=1.0000000 frame=[11.84 sec]");
				
				selectImage("maskobjects");

				

				setSlice(j+1);

				run("Analyze Particles...", "size=500-Infinity pixel show=Nothing add slice");

				
				nROIS=roiManager("count");

				for (l = 0; l < nROIS; l++) {
		
					selectImage("objects");

					roiManager("select", l);

					roiManager("measure");

					cell=getResult("Max", 0);

					selectImage("Red2");
					
					roiManager("select", l);;

					roiManager("measure");
			
					intensityredcell[cell-1]=getResult("RawIntDen",1);
					areacell[cell-1]=getResult("Area",1);


					run("Set Measurements...", "area integrated limit display redirect=None decimal=5");
					
					IJ.deleteRows(0, nResults);
	
					selectImage("Red4"+j);

					roiManager("select", l);

					setThreshold(1, 255);

					run("Measure");

					//waitForUser("");
					
					
					numpointscell[cell-1]=getResult("Area",0);

		
					//File.append( l + "\t" + areacell + "\t" + Areacelltotal + "\t" + Intensitytotal + "\t" + numpointstotal + "\t" + Areaimagetotal + "\t" + Intensityimagetotal , path);

					IJ.deleteRows(0, nResults);

					run("Set Measurements...", "area mean min integrated display redirect=None decimal=5");


				}
				

				print(nROIS);
				print(roiManager("count"));
				for (l = 0; l < nROIS; l++) {
					roiManager("Select", l);
					run("Enlarge...", "enlarge=-4 pixel");
					//waitForUser("");
					roiManager("Update")
				}


				
				for (l = 0; l < nROIS; l++) {
		
					selectImage("objects");

					roiManager("select", l);

					roiManager("measure");

					cell=getResult("Max", 0);

					selectImage("Red2");
					
					roiManager("select", l);;

					roiManager("measure");
			
					intensityredmem[cell-1]=getResult("RawIntDen",1);
					areamem[cell-1]=getResult("Area",1);


					run("Set Measurements...", "area integrated limit display redirect=None decimal=5");
					
					IJ.deleteRows(0, nResults);
	
					selectImage("Red4"+j);

					roiManager("select", l);

					setThreshold(1, 255);

					run("Measure");

					//waitForUser("");
					
					
					numpointsmem[cell-1]=getResult("Area",0);

		
					//File.append( l + "\t" + areacell + "\t" + Areacelltotal + "\t" + Intensitytotal + "\t" + numpointstotal + "\t" + Areaimagetotal + "\t" + Intensityimagetotal , path);

					IJ.deleteRows(0, nResults);

					run("Set Measurements...", "area mean min integrated display redirect=None decimal=5");


				}
				


				for (l = 0; l < ncells; l++) {

		
					File.append( j+1 + "\t" + l+1 + "\t" + areacell[l] + "\t" + intensityredcell[l] + "\t" + numpointscell[l] + "\t" + areamem[l] + "\t" + intensityredmem[l] + "\t" + numpointsmem[l], path);

				}
			
			}


			//File.append( t + "\t" + cells + "\t" + j+1 + "\t" + Areatotal + "\t" + Areacelltotal + "\t" + Intensitytotal + "\t" + numpointstotal + "\t" + Areaimagetotal + "\t" + Intensityimagetotal , path);

			IJ.deleteRows(0, nResults);


		closeImagesWindows();

		}

		
	
	}


}




waitForUser("Macro has finished");

function closeImagesWindows(){
	run("Close All");
	if(isOpen("Results")){
		selectWindow("Results");
		run("Close");
	}
	if(isOpen("ROI Manager")){
		selectWindow("ROI Manager");
		run("Close");
	}
	if(isOpen("Threshold")){
		selectWindow("Threshold");
		run("Close");
	}
	if(isOpen("Summary")){
		selectWindow("Summary");
		run("Close");
	}
	if(isOpen("B&C")){
		selectWindow("B&C");
		run("Close");
	}

}


	
