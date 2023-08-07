import ij.IJ
import ij.ImagePlus
import ij.gui.PolygonRoi
import ij.gui.Roi
import ij.gui.ShapeRoi
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.frame.RoiManager
import inra.ijpb.binary.BinaryImages
import loci.plugins.BF
import loci.plugins.in.ImporterOptions
import org.apache.commons.compress.utils.FileNameUtils
import inra.ijpb.morphology.Strel


// INPUT UI
//
#@File(label = "Input File Directory", style = "directory") inputFilesDir
#@File(label = "Output directory", style = "directory") outputDir
#@Integer(label = "Nuclei Channel", value = 0) nucleiChannel
#@Integer(label = "Cyto Channel", value = 4) cytoChannel
//#@Boolean(label = "Apply DAPI?") applyDAPI


// IDE
//
//
//def headless = true;
//new ImageJ().setVisible(true);

IJ.log("-Parameters selected: ")
IJ.log("    -inputFileDir: " + inputFilesDir)
IJ.log("    -outputDir: " + outputDir)
IJ.log("                                                           ");
def rm = new RoiManager();
/** Get files (images) from input directory */
def listOfFiles = inputFilesDir.listFiles();
/** Define table to store results */
def tableTotal = new ResultsTable()
def analyzeLif = null
if (listOfFiles[0].getName().contains(".lif")) {
    analyzeLif = true
}
for (def i = 0; i < listOfFiles.length; i++) {

    if (!listOfFiles[i].getName().contains("DS")) {
        IJ.log("Analyzing file: " + listOfFiles[i].getName());
        /** Define output directory per file */

        IJ.log("    -Creating output dir in " + outputDir.getAbsolutePath());

        if (FileNameUtils.getExtension(listOfFiles[i].getName()).contains("lif")) {


            /** Importer options for .lif file */
            def options = new ImporterOptions();
            options.setId(inputFilesDir.getAbsolutePath() + File.separator + listOfFiles[i].getName());
            options.setSplitChannels(false);
            options.setSplitTimepoints(false);
            options.setSplitFocalPlanes(false);
            options.setAutoscale(true);
            options.setStackFormat(ImporterOptions.VIEW_HYPERSTACK);
            options.setStackOrder(ImporterOptions.ORDER_XYCZT);
            options.setColorMode(ImporterOptions.COLOR_MODE_COMPOSITE);
            options.setCrop(false);
            options.setOpenAllSeries(true);
            def imps = BF.openImagePlus(options);
            /** Define results table per image */
            def tableLif = new ResultsTable();

            for (int j = 0; j < imps.length; j++) {
                IJ.log("        -Analyzing serie: " + (j + 1).toString())

                /** Declare each image to process within input directory */
                def imp = imps[j];
                def impTitleSerie = null;
                if (imp.getTitle().contains("/")) {
                    impTitleSerie = imp.getTitle().replaceAll("/", "_");
                } else {
                    impTitleSerie = imp.getTitle();
                }
              
                /** Define results table per slice */
                def tableSlice = new ResultsTable()

                /** Get calibration from non-transformed image. */
                def cal = imp.getCalibration();

                /** Define all variables to store results per slice */
                def nucleiSum = new ArrayList<Double>()
                def cytoSum = new ArrayList<Double>();
                def nucleiMean = new ArrayList<Double>()
                def cytoMean = new ArrayList<Double>();
                def cytoInverseSum = new ArrayList<Double>();
                def redSum = new ArrayList<Double>();
                def redMean = new ArrayList<Double>();

                /** Split channels */
                def channels = ChannelSplitter.split(imp);
                /** Get channel Blue */
                def chNuclei = channels[nucleiChannel.intValue()];
                /** Get channel Red */
                def chCyto = channels[cytoChannel.intValue()];
                
                for (def s = 1; s <= chNuclei.stackSize; s++) {
                    /** Analysis of nuclei slices */
                    def ipNucleiSlice = chNuclei.getImageStack().getProcessor(s);
                    def chNucleiSlice = new ImagePlus("", ipNucleiSlice);
                    def chNucleiSliceMeasure = chNucleiSlice.duplicate();
                    /** Segment Nuclei Areas */
                    IJ.run(chNucleiSlice, "Auto Threshold", "method=Otsu ignore_black white");
                    /** Apply median filter to remove background areas below
                     5 pixels */
                    IJ.run(chNucleiSlice, "Median...", "radius=3")
                    /** Apply area opening in binary images */
                    chNucleiSlice = new ImagePlus("", BinaryImages.areaOpening(chNucleiSlice.getProcessor(), 500));
                    //chNucleiSlice.show()
                    /** Apply area closing in binary images */
                    Strel strel = Strel.Shape.DISK.fromRadius(4);
                    chNucleiSlice = new ImagePlus("", strel.closing(chNucleiSlice.getProcessor()));
                    /** Create contour selection on nuclei areas */
                    IJ.run(chNucleiSlice, "Create Selection", "");
                    //chNucleiSlice.show()
                    /** Get contour as ROI */
                    def roiNucleiSlice = chNucleiSlice.getRoi();
                    //rm.addRoi(roiNucleiSlice);

                    /** Analysis of cyto slices */
                    def ipCytoSlice = chCyto.getImageStack().getProcessor(s);
                    def chCytoSlice = new ImagePlus("", ipCytoSlice);
                    def chCytoSliceMeasure = chCytoSlice.duplicate();
                    /** Segment Cyto Areas */
                    IJ.run(chCytoSlice, "Auto Threshold", "method=Otsu ignore_black white");
                    /** Apply median filter to remove background areas below
                     5 pixels */
                    IJ.run(chCytoSlice, "Median...", "radius=3")
                    /** Apply area closing in binary images */
   
                    chCytoSlice = new ImagePlus("", strel.closing(chCytoSlice.getProcessor()));
                    /** Create contour selection on cyto areas */
                    IJ.run(chCytoSlice, "Create Selection", "");
                    /** Get contour as ROI */
                    def roiCytoSlice = chCytoSlice.getRoi();
                    //rm.addRoi(roiCytoSlice)
                    //chCytoSlice.show()
                    /** NUCLEI ANALYSIS (measure nuclei on red areas) */
                    def roiNuclei = null;
                    def meanNuclei = null;
                    def sumNuclei = null;
                    if (roiNucleiSlice == null) {
                        roiNuclei = null;
                        meanNuclei = 0.0;
                        sumNuclei = 0.0
                    } else {
                        roiNuclei = roiNucleiSlice;
                        /** Measure nuclei intensity on cyto channel */
                        chCytoSliceMeasure.setRoi(roiNuclei);
                        /** Store statistics relative to intensity (mean,max,min,sum...) */
                        meanNuclei = chCytoSliceMeasure.getStatistics().mean;
                        /** Store statistics sum intensity */
                        sumNuclei = meanNuclei * roiNuclei.getStatistics().area;
                    }
                    /** CYTO ANALYSIS (measure cyto (cyto without nuclei) on red areas) */
                    def roiCyto = null;
                    def meanCyto = null;
                    def sumCyto = null;
                    if (roiCytoSlice == null || roiNucleiSlice == null) {
                        roiCyto = null;
                        meanCyto = 0.0;
                        sumCyto = 0.0
                    } else {
                        def roiCytoAND = new ShapeRoi(roiNuclei).and(new ShapeRoi(roiCytoSlice));
                        roiCyto = roiCytoAND.xor(new ShapeRoi(roiCytoSlice)).shapeToRoi();
                        //rm.addRoi(roiCyto)
                        def areaCyto = null;
                        if(roiCyto == null){
                        	areaCyto = 0.0
                        }else{
                        	areaCyto = roiCyto.getStatistics().area;
                        	}
                        //rm.addRoi(roiCyto)
                        /** Measure nuclei intensity on cyto channel */
                        //chCytoSliceMeasure.show()
                        chCytoSliceMeasure.setRoi(roiCyto);
                        /** Store statistics relative to intensity (mean,max,min,sum...) */
                        meanCyto = chCytoSliceMeasure.getStatistics().mean;
                        /** Store statistics sum intensity */
                        sumCyto = meanCyto * areaCyto;
                    }
                    /** RED ANALYSIS (measure red on red areas) */
                    def roiRed = null;
                    def sumRed = null;
                    def meanRed = null;
                    if (roiCytoSlice == null) {
                        roiRed = null;
                        sumRed = 0.0
                    } else {
                        roiRed = roiCytoSlice;
                        /** Measure red intensity on red channel */
                        chCytoSliceMeasure.setRoi(roiRed);
                        /** Store statistics sum intensity */
                        sumRed = chCytoSliceMeasure.getStatistics().mean * roiRed.getStatistics().area;
                        meanRed = chCytoSliceMeasure.getStatistics().mean;
                    }

                    /** NUCLEI ANALYSIS CONTAINED ON RED AREAS (measure nuclei contained on red areas) */
                    def sumNucleiSplit = new ArrayList<Double>();
                    def meanNucleiSplit = new ArrayList<Double>();
                    if(roiNuclei != null || roiCytoSlice !=null) {
                          def roisNucleiSplit = new ShapeRoi(roiNuclei).getRois();
                        for (def n = 0.0.intValue(); n < roisNucleiSplit.length; n++) {
                                    if (roiCytoSlice.containsPoint(roisNucleiSplit[n].getContourCentroid()[0], roisNucleiSplit[n].getContourCentroid()[1])) {
                                        chCytoSliceMeasure.setRoi(roisNucleiSplit[n]);
                                        sumNucleiSplit.add(chCytoSliceMeasure.getStatistics().mean * roisNucleiSplit[n].getStatistics().area);
                                         meanNucleiSplit.add(chCytoSliceMeasure.getStatistics().mean);
                                
                            }
                        }
                    }else{
                        sumNucleiSplit.add(0.0.toDouble());
                    }
					sumRed = sumNucleiSplit.stream().mapToDouble(a -> a).sum().toDouble() + sumCyto;
                    /**Create results table per image (per slice results) */
                    tableSlice.incrementCounter();
                    tableSlice.setValue("Plane (Slice)", s, s.toString())
                    tableSlice.setValue("Cyto Sum Intensity", s, sumCyto.toString())
                    tableSlice.setValue("Cyto Mean Intensity", s, meanCyto.toString())
                    tableSlice.setValue("Positive Nuclei Sum Intensity", s, sumNucleiSplit.stream().mapToDouble(a -> a).sum().toDouble())
                    tableSlice.setValue("Positive Nuclei Mean Intensity", s, meanNucleiSplit.stream().mapToDouble(a -> a).sum().toDouble())
                    tableSlice.setValue("Total Sum Intensity", s, sumRed.toString())
                     tableSlice.setValue("Total Mean Intensity", s, meanRed.toString())
                    /**Store measurements to fit table per image */
                    //nucleiSum.add(sumNuclei.toDouble())
                    cytoSum.add(sumCyto.toDouble())
                    cytoMean.add(meanCyto.toDouble())
                    redSum.add(sumRed.toDouble())
                    redMean.(meanRed.toDouble())
                    cytoInverseSum.add(sumNucleiSplit.stream().mapToDouble(a -> a).sum().toDouble())
                    nucleiMean.add(meanNucleiSplit.stream().mapToDouble(a -> a).average().orElse(0.0))
                    
				  //rm.runCommand("Save",
                       //outputDir.getAbsolutePath() + File.separator +imp.getTitle().replaceAll(".lif", "")+"slice_"+(s+1) + "_RoiSet.zip")
					
                }
                //rm.reset()
                  tableSlice.setValue("Cyto Sum Intensity", chNuclei.stackSize+1, cytoSum.stream()
                        .mapToDouble(a -> a)
                        .sum().toDouble().toString())
                 tableSlice.setValue("Cyto Mean Intensity", chNuclei.stackSize+1, cytoMean.stream()
                        .mapToDouble(a -> a)
                        .average().orElse(0.0).toString())
                tableSlice.setValue("Positive Nuclei Sum Intensity", chNuclei.stackSize+1, cytoInverseSum.stream()
                        .mapToDouble(a -> a)
                        .sum().toDouble().toString())
                        tableSlice.setValue("Positive Nuclei Mean Intensity", chNuclei.stackSize+1, nucleiMean.stream()
                        .mapToDouble(a -> a)
                        .average().orElse(0.0).toString())
                tableSlice.setValue("Total Sum Intensity", chNuclei.stackSize+1, redSum.stream()
                        .mapToDouble(a -> a)
                        .sum().toDouble().toString())
               tableSlice.setValue("Total Mean Intensity", chNuclei.stackSize+1, redMean.stream()
                        .mapToDouble(a -> a)
                        .average().orElse(0.0).toString())
                /**Save table per slice (plane) */
                def tablePath = new File(outputDir, impTitleSerie.replaceAll(".tif", "") + "_" + "table_slice_results" + ".csv").toString();
                IJ.log("Saving table per serie: " + tablePath + " in " + outputDir.getAbsolutePath());
                tableSlice.save(tablePath);
				
                /** Iterate through table per image (.Lif) */
               tableLif.setValue("Image Serie", j, imp.getTitle().replaceAll(".lif", ""))
                tableLif.setValue("N Slices (Planes)", j, chNuclei.stackSize.toString())
                tableLif.setValue("Cyto Sum Intensity", j, cytoSum.stream()
                        .mapToDouble(a -> a)
                        .sum().toDouble().toString())
                 tableLif.setValue("Cyto Mean Intensity", j, cytoMean.stream()
                        .mapToDouble(a -> a)
                        .average().orElse(0.0).toString())
                tableLif.setValue("Positive Nuclei Sum Intensity", j, cytoInverseSum.stream()
                        .mapToDouble(a -> a)
                        .sum().toDouble().toString())
               tableLif.setValue("Positive Nuclei Mean Intensity", j, nucleiMean.stream()
                        .mapToDouble(a -> a)
                        .average().orElse(0.0).toString())
                tableLif.setValue("Total Sum Intensity", j, redSum.stream()
                        .mapToDouble(a -> a)
                        .sum().toDouble().toString())
               tableLif.setValue("Total Mean Intensity", j, redMean.stream()
                        .mapToDouble(a -> a)
                        .average().orElse(0.0).toString())


            
            def tablePathLif = new File(outputDir, listOfFiles[i].getName().replaceAll(".lif", "") + "_" + "table_perimage_results" + ".csv").toString();
            IJ.log("Saving table per file: " + tablePathLif + " in " + outputDir.getAbsolutePath());
            tableLif.save(tablePathLif);
           }

        } else {

            /** Create results table to store results per slice (section or plane) */
            def tableSlice = new ResultsTable()
            /** Create image for each file in the input directory */
            def imps = new ImagePlus(listOfFiles[i].getAbsolutePath())

            /** Define all variables to store results per slice */
            def nucleiSum = new ArrayList<Double>()
            def cytoSum = new ArrayList<Double>();
            def cytoInverseSum = new ArrayList<Double>();
            def redSum = new ArrayList<Double>();


            IJ.log("        -Analyzing image: " + imps.getTitle())
            /** Get original calibration from image. */
            def cal = imps.getCalibration();
            /** Split channels */
            def channels = ChannelSplitter.split(imps);
            /** Get channel Blue */
            def chNuclei = channels[0];
            /** Get channel Red */
            def chCyto = channels[4];


            for (def s = 1; s <= chNuclei.stackSize; s++) {
                /** Analysis of nuclei slices */
                def ipNucleiSlice = chNuclei.getImageStack().getProcessor(s);
                def chNucleiSlice = new ImagePlus("", ipNucleiSlice);
                def chNucleiSliceMeasure = chNucleiSlice.duplicate();
                /** Segment Nuclei Areas */
                IJ.run(chNucleiSlice, "Auto Threshold", "method=RenyiEntropy ignore_black white");
                /** Apply median filter to remove background areas below
                 5 pixels */
                IJ.run(chNucleiSlice, "Median...", "radius=3")
                /** Apply area opening in binary images */
                chNucleiSlice = new ImagePlus("", BinaryImages.areaOpening(chNucleiSlice.getProcessor(), 500));
                //chNucleiSlice.show()
                /** Apply area closing in binary images */
                Strel strel = Strel.Shape.DISK.fromRadius(4);
                chNucleiSlice = new ImagePlus("", strel.closing(chNucleiSlice.getProcessor()));
                /** Create contour selection on nuclei areas */
                IJ.run(chNucleiSlice, "Create Selection", "");
                //chNucleiSlice.show()
                /** Get contour as ROI */
                def roiNucleiSlice = chNucleiSlice.getRoi();

                /** Analysis of cyto slices */
                def ipCytoSlice = chCyto.getImageStack().getProcessor(s);
                def chCytoSlice = new ImagePlus("", ipCytoSlice);
                def chCytoSliceMeasure = chCytoSlice.duplicate();
                /** Segment Cyto Areas */
                IJ.run(chCytoSlice, "Auto Threshold", "method=Otsu ignore_black white");
                /** Apply median filter to remove background areas below
                 5 pixels */
                IJ.run(chCytoSlice, "Median...", "radius=5")
                /** Create contour selection on cyto areas */
                IJ.run(chCytoSlice, "Create Selection", "");
                /** Get contour as ROI */
                def roiCytoSlice = chCytoSlice.getRoi();
                //chCytoSlice.show()
                /** NUCLEI ANALYSIS (measure nuclei on red areas) */
                def roiNuclei = null;
                def meanNuclei = null;
                def sumNuclei = null;
                if (roiNucleiSlice == null) {
                    roiNuclei = null;
                    meanNuclei = 0.0;
                    sumNuclei = 0.0
                } else {
                    roiNuclei = roiNucleiSlice;
                    /** Measure nuclei intensity on cyto channel */
                    chCytoSliceMeasure.setRoi(roiNuclei);
                    /** Store statistics relative to intensity (mean,max,min,sum...) */
                    meanNuclei = chCytoSliceMeasure.getStatistics().mean;
                    /** Store statistics sum intensity */
                    sumNuclei = meanNuclei * roiNuclei.getStatistics().area;
                }
                /** CYTO ANALYSIS (measure cyto (cyto without nuclei) on red areas) */
                def roiCyto = null;
                def meanCyto = null;
                def sumCyto = null;
                if (roiCytoSlice == null || roiNucleiSlice == null) {
                    roiCyto = null;
                    meanCyto = 0.0;
                    sumCyto = 0.0
                } else {
                    def roiCytoAND = new ShapeRoi(roiNuclei).and(new ShapeRoi(roiCytoSlice));
                    roiCyto = roiCytoAND.xor(new ShapeRoi(roiCytoSlice)).shapeToRoi();
                    //rm.addRoi(roiCyto)
                    /** Measure nuclei intensity on cyto channel */
                    //chCytoSliceMeasure.show()
                    chCytoSliceMeasure.setRoi(roiCyto);
                    /** Store statistics relative to intensity (mean,max,min,sum...) */
                    meanCyto = chCytoSliceMeasure.getStatistics().mean;
                    /** Store statistics sum intensity */
                    sumCyto = meanCyto * roiCyto.getStatistics().area;
                }
                /** RED ANALYSIS (measure red on red areas) */
                def roiRed = null;
                def sumRed = null;
                if (roiCytoSlice == null) {
                    roiRed = null;
                    sumRed = 0.0
                } else {
                    roiRed = roiCytoSlice;
                    /** Measure red intensity on red channel */
                    chCytoSliceMeasure.setRoi(roiRed);
                    /** Store statistics sum intensity */
                    sumRed = chCytoSliceMeasure.getStatistics().mean * roiRed.getStatistics().area;
                }

                /** NUCLEI ANALYSIS CONTAINED ON RED AREAS (measure nuclei contained on red areas) */
                def sumNucleiSplit = new ArrayList<Double>();
                if(roiNuclei != null) {
                        def roisNucleiSplit = new ShapeRoi(roiNuclei).getRois();
                        for (def n = 0.0.intValue(); n < roisNucleiSplit.length; n++) {
                            for (def x = 0.0.intValue(); x <roisNucleiSplit[n].getConvexHull().xpoints.length; x++) {
                                for (def y = 0.0.intValue(); y <roisNucleiSplit[n].getConvexHull().ypoints.length; y++) {
                                    if (roiCytoSlice.contains(roisNucleiSplit[n].getConvexHull().xpoints[x], roisNucleiSplit[n].getConvexHull()yxpoints[y])) {
                                        chCytoSliceMeasure.setRoi(roisNucleiSplit[n]);
                                        sumNucleiSplit.add(chCytoSliceMeasure.getStatistics().mean * roisNucleiSplit[n].getStatistics().area);
                                    }
                                }
                            }
                        }
                }else{
                    sumNucleiSplit.add(0.0.toDouble());
                }
                /**Iterate through table per slice */
                tableSlice.incrementCounter();
                tableSlice.setValue("Plane (Slice)", s, s.toString())
                tableSlice.setValue("Nuclei Sum Intensity", s, sumNuclei.toString())
                tableSlice.setValue("Cyto Sum Intensity", s, sumCyto.toString())
                tableSlice.setValue("Inverse Cyto Sum Intensity", s, sumNucleiSplit.stream().mapToDouble(a -> a).sum().toDouble())
                tableSlice.setValue("Red Sum Intensity", s, sumRed.toString())
                /**Store measurements to fit table per image */
                nucleiSum.add(sumNuclei.toDouble())
                cytoSum.add(sumCyto.toDouble())
                redSum.add(sumRed.toDouble())
                cytoInverseSum.add(sumNucleiSplit.stream().mapToDouble(a -> a).sum().toDouble())


            }
            /**Store table per slice (plane) */
            def tablePath = new File(outputDir, imps.getTitle().replaceAll(".tif", "") + "_" + "table_image_results" + ".csv").toString();
            IJ.log("Saving table per image: " + tablePath + " in " + outputDir.getAbsolutePath());
            tableSlice.save(tablePath);


            tableTotal.incrementCounter()
            /** Iterate per image */
            tableTotal.setValue("Image Title", i, imps.getTitle().replaceAll(".tif", "").toString())
            tableTotal.setValue("N Slices (Planes)", i, chNuclei.stackSize.toString())
            tableTotal.setValue("Nuclei Sum Intensity Average", i, nucleiSum.stream()
                    .mapToDouble(a -> a)
                    .average().orElse(0.0).toString())
            tableTotal.setValue("Cyto Sum Intensity Average", i, cytoSum.stream()
                    .mapToDouble(a -> a)
                    .average().orElse(0.0).toString())
            tableTotal.setValue("Cyto Inverse Sum Intensity Average", i, cytoInverseSum.stream()
                    .mapToDouble(a -> a)
                    .average().orElse(0.0).toString())
            tableTotal.setValue("Red Sum Intensity Average", i, redSum.stream()
                    .mapToDouble(a -> a)
                    .average().orElse(0.0).toString())


        }
    }

}
if (!analyzeLif) {
    /**Save table per image */
    def tablePath = new File(outputDir, "table_image_results" + ".csv").toString();
    IJ.log("Saving table per image: " + tablePath + " in " + outputDir.getAbsolutePath());
    tableTotal.save(tablePath);
}


IJ.log("Done!!!")


