import ij.IJ
import ij.ImagePlus
import ij.gui.ShapeRoi
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.frame.RoiManager
import inra.ijpb.binary.BinaryImages
import loci.plugins.BF
import loci.plugins.in.ImporterOptions
import java.io.File;

// INPUT UI
//
#
@File(label = "Input File Directory", style = "directory") inputFilesDir
#
@File(label = "Output directory", style = "directory") outputDir
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
        def outputImageDir = new File(outputDir.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".lif", ""));

        if (!outputImageDir.exists()) {
            def results = false;

            try {
                outputImageDir.mkdir();
                results = true;
            } catch (SecurityException se) {
            }
        }
        IJ.log("    -Creating output dir for image " + listOfFiles[i].getName().replaceAll(".lif", "") + " in " + outputDir.getAbsolutePath());
        def imps = null;
        if (listOfFiles[i].getName().contains(".lif")) {
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
            imps = BF.openImagePlus(options);
            /** Define results table per image */
            def tableLif = new ResultsTable();

            for (int j = 0; j < imps.length; j++) {
                IJ.log("        -Analyzing serie: " + (j + 1).toString())

                /** Declare each image to process within input directory */
                def imp = imps[j];
                def impTitleSerie = null;
                if (imp.getTitle().contains("/")) {
                    impTitleSerie = imp.getTitle().replaceAll("/", "");
                } else {
                    impTitleSerie = imp.getTitle();
                }
                /** Define results table per slice */
                def tableSlice = new ResultsTable()
                imps = imp;
                /** Get calibration from non-transformed image. */
                def cal = imps.getCalibration();

                /** Split channels */
                def channels = ChannelSplitter.split(imps);
                /** Get channel Blue */
                def chBlue = channels[0];
                /** Get channel Red */
                def chRed = channels[4];
                /** Define all variables to store results per slice */
                def nucleiMean = new ArrayList<Double>();
                def nucleiMax = new ArrayList<Double>();
                def nucleiMin = new ArrayList<Double>();
                def nucleiSum = new ArrayList<Double>();
                def nucleiArea = new ArrayList<Double>();
                def sumRedTotal = new ArrayList<Double>();
                def cytoMean = new ArrayList<Double>();


                for (def s = 1; s <= chBlue.stackSize; s++) {
                    /** Analysis of blue slices */
                    def ipBlueSlice = chBlue.getImageStack().getProcessor(s);
                    def chBlueSlice = new ImagePlus("", ipBlueSlice);
                    def chBlueSliceMeasure = chBlueSlice.duplicate();
                    /** Segment Blue Areas */
                    IJ.run(chBlueSlice, "Auto Threshold", "method=Huang ignore_black white");
                    /** Apply median filter to remove background areas below
                     5 pixels */
                    IJ.run(chBlueSlice, "Median...", "radius=5")
                    /** Apply area opening in binary images */
                    chBlueSlice = new ImagePlus("", BinaryImages.areaOpening(chBlueSlice.getProcessor(), 100));
                    /** Create contour selection on blue areas */
                    IJ.run(chBlueSlice, "Create Selection", "");
                    /** Get contour as ROI */
                    def roiBlueSlice = chBlueSlice.getRoi();

                    /** Analysis of red slices */
                    def ipRedSlice = chRed.getImageStack().getProcessor(s);
                    def chRedSlice = new ImagePlus("", ipRedSlice);
                    def chRedSliceMeasure = chRedSlice.duplicate();
                    /** Segment Red Areas */
                    IJ.run(chRedSlice, "Auto Threshold", "method=Otsu ignore_black white");
                    /** Apply median filter to remove background areas below
                     5 pixels */
                    IJ.run(chRedSlice, "Median...", "radius=3")
                    /** Create contour selection on red areas */
                    IJ.run(chRedSlice, "Create Selection", "");
                    /** Get contour as ROI */
                    def roiRedSlice = chRedSlice.getRoi();
                    /** NUCLEI ANALYSIS (blue overlapping red) */
                    /** Apply AND operator to keep those blue areas overlapping with red areas */
                    def roiBlueRed = new ShapeRoi(roiBlueSlice).and(new ShapeRoi(roiRedSlice)).shapeToRoi();
                    /** Measure total overlapping blue-red intensity on red channel */
                    chRedSliceMeasure.setRoi(roiBlueRed);
                    /** Store statistics relative to coloc intensity (mean and sum...) */
                    def meanColoc = chRedSliceMeasure.getStatistics().mean
                    def sumColoc = chRedSliceMeasure.getStatistics().mean * roiBlueRed.getStatistics().area;
                    def minColoc = chRedSliceMeasure.getStatistics().min
                    def maxColoc = chRedSliceMeasure.getStatistics().max

                    /** Measure total red intensity on red channel */
                    chRedSliceMeasure.setRoi(roiRedSlice)
                    /** Store statistics relative to red intensity (mean and sum...) */
                    def meanRed = chRedSliceMeasure.getStatistics().mean
                    def sumRed = chRedSliceMeasure.getStatistics().mean * roiRedSlice.getStatistics().area;

                    /**Create results table per image (per slice results) */
                    tableSlice.incrementCounter();
                    tableSlice.setValue("Plane (Slice)", s, s.toString())
                    tableSlice.setValue("Nuclei Mean Fluorescence", s, meanColoc.toString())
                    tableSlice.setValue("Nuclei Max Fluorescence", s, maxColoc.toString())
                    tableSlice.setValue("Nuclei Min Fluorescence", s, minColoc.toString())
                    tableSlice.setValue("Nuclei Sum Fluorescence", s, (meanColoc * roiBlueRed.getStatistics().area).toString())
                    tableSlice.setValue("Nuclei Area", s, roiBlueRed.getStatistics().area)
                    tableSlice.setValue("Sum Red", s, sumRed)
                    tableSlice.setValue("Cyto Sum Fluorescence", s, (sumRed - sumColoc).toString())
                    /**Store table per image (per slice results) */
                    nucleiMean.add(meanColoc.toDouble())
                    nucleiMax.add(maxColoc.toDouble())
                    nucleiMin.add(minColoc.toDouble())
                    nucleiSum.add((meanColoc * roiBlueRed.getStatistics().area).toDouble())
                    nucleiArea.add(roiBlueRed.getStatistics().area)
                    sumRedTotal.add(sumRed)
                    cytoMean.add((sumRed - sumColoc).toDouble())


                }
                /**Save table per slice (plane) */
                def tablePath = new File(outputDir, impTitleSerie.replaceAll(".tif", "") + "_" + "table_slice_results" + ".csv").toString();
                IJ.log("Saving table per serie: " + tablePath + " in " + outputDir.getAbsolutePath());
                tableSlice.save(tablePath);

                /** Iterate through table per image (.Lif) */
                tableLif.setValue("Image Serie", j, imps.getTitle().replaceAll(".lif",""))
                tableLif.setValue("N Slices (Planes)", j, chBlue.stackSize.toString())
                tableLif.setValue("Nuclei Mean Fluorescence", j,nucleiMean.stream()
                        .mapToDouble(a -> a)
                        .average().getAsDouble().toString())
                tableLif.setValue("Nuclei Max Fluorescence", j, nucleiMax.stream()
                        .mapToDouble(a -> a)
                        .average().getAsDouble())
                tableLif.setValue("Nuclei Min Fluorescence", j, nucleiMin.stream()
                        .mapToDouble(a -> a)
                        .average().getAsDouble().toString())
                tableLif.setValue("Nuclei Sum Fluorescence", j, nucleiSum.stream()
                        .mapToDouble(a -> a)
                        .average().getAsDouble().toString())
                tableLif.setValue("Nuclei Area", j, nucleiArea.stream()
                        .mapToDouble(a -> a)
                        .average().getAsDouble().toString())
                tableLif.setValue("Sum Red",j, sumRed.stream()
                        .mapToDouble(a -> a)
                        .average().getAsDouble().toString())
                tableLif.setValue("Cyto Mean Fluorescence", j, cytoMean.stream()
                        .mapToDouble(a -> a)
                        .average().getAsDouble().toString())

                def tablePathLif = new File(outputDir, listOfFiles[i].getName().replaceAll(".lif", "") + "_" + "table_perimage_results" + ".csv").toString();
                IJ.log("Saving table per file: " + tablePathLif + " in " + outputDir.getAbsolutePath());
                tableLif.save(tablePathLif);

            }



        } else {

            /** Create results table to store results per slice (section or plane) */
            def tableSlice = new ResultsTable()
            /** Create image for each file in the input directory */
            imps = new ImagePlus(listOfFiles[i].getAbsolutePath())

            /** Define all variables to store results per slice */
            def nucleiMean = new ArrayList<Double>();
            def nucleiMax = new ArrayList<Double>();
            def nucleiMin = new ArrayList<Double>();
            def nucleiSum = new ArrayList<Double>();
            def nucleiArea = new ArrayList<Double>();
            def sumRedTotal = new ArrayList<Double>();
            def cytoMean = new ArrayList<Double>();

            IJ.log("        -Analyzing image: " + imps.getTitle())
            /** Get original calibration from image. */
            def cal = imps.getCalibration();
            /** Split channels */
            def channels = ChannelSplitter.split(imps);
            /** Get channel Blue */
            def chBlue = channels[0];
            /** Get channel Red */
            def chRed = channels[4];


            for (def s = 1; s <= chBlue.stackSize; s++) {
                /** Analysis of blue slices */
                def ipBlueSlice = chBlue.getImageStack().getProcessor(s);
                def chBlueSlice = new ImagePlus("", ipBlueSlice);
                def chBlueSliceMeasure = chBlueSlice.duplicate();
                /** Segment Blue Areas */
                IJ.run(chBlueSlice, "Auto Threshold", "method=Huang ignore_black white");
                /** Apply median filter to remove background areas below
                 5 pixels */
                IJ.run(chBlueSlice, "Median...", "radius=5")
                /** Apply area opening in binary images */
                chBlueSlice = new ImagePlus("", BinaryImages.areaOpening(chBlueSlice.getProcessor(), 100));
                /** Create contour selection on blue areas */
                IJ.run(chBlueSlice, "Create Selection", "");
                /** Get contour as ROI */
                def roiBlueSlice = chBlueSlice.getRoi();

                /** Analysis of red slices */
                def ipRedSlice = chRed.getImageStack().getProcessor(s);
                def chRedSlice = new ImagePlus("", ipRedSlice);
                def chRedSliceMeasure = chRedSlice.duplicate();
                /** Segment Red Areas */
                IJ.run(chRedSlice, "Auto Threshold", "method=Otsu ignore_black white");
                /** Apply median filter to remove background areas below
                 5 pixels */
                IJ.run(chRedSlice, "Median...", "radius=3")
                /** Create contour selection on red areas */
                IJ.run(chRedSlice, "Create Selection", "");
                /** Get contour as ROI */
                def roiRedSlice = chRedSlice.getRoi();
                /** NUCLEI ANALYSIS (blue overlapping red) */
                /** Apply AND operator to keep those blue areas overlapping with red areas */
                def roiBlueRed = new ShapeRoi(roiBlueSlice).and(new ShapeRoi(roiRedSlice)).shapeToRoi();
                /** Measure intensity on overlapping Blue-Red on red channel */
                chRedSliceMeasure.setRoi(roiBlueRed);
                /** Store statistics relative to intensity (mean,max,min,sum...) */
                def meanColoc = chRedSliceMeasure.getStatistics().mean
                def sumColoc = chRedSliceMeasure.getStatistics().mean * roiBlueRed.getStatistics().area;
                def minColoc = chRedSliceMeasure.getStatistics().min
                def maxColoc = chRedSliceMeasure.getStatistics().max
                /** Measure total red intensity on red channel */
                chRedSliceMeasure.setRoi(roiRedSlice)
                /** Store statistics relative to red intensity (mean and sum...) */
                def meanRed = chRedSliceMeasure.getStatistics().mean
                def sumRed = chRedSliceMeasure.getStatistics().mean * roiRedSlice.getStatistics().area;

                /**Iterate through table per slice */
                tableSlice.incrementCounter();
                tableSlice.setValue("Plane (Slice)", s, s.toString())
                tableSlice.setValue("Nuclei Mean Fluorescence", s, meanColoc.toString())
                tableSlice.setValue("Nuclei Max Fluorescence", s, maxColoc.toString())
                tableSlice.setValue("Nuclei Min Fluorescence", s, minColoc.toString())
                tableSlice.setValue("Nuclei Sum Fluorescence", s, (meanColoc * roiBlueRed.getStatistics().area).toString())
                tableSlice.setValue("Nuclei Area", s, roiBlueRed.getStatistics().area)
                tableSlice.setValue("Sum Red", s, sumRed);
                tableSlice.setValue("Cyto Sum Fluorescence", s, (sumRed - sumColoc).toString())
                /**Store measurements to fit table per image */
                nucleiMean.add(meanColoc.toDouble())
                nucleiMax.add(maxColoc.toDouble())
                nucleiMin.add(minColoc.toDouble())
                nucleiSum.add((meanColoc * roiBlueRed.getStatistics().area).toDouble())
                nucleiArea.add(roiBlueRed.getStatistics().area)
                sumRedTotal.add(sumRed)
                cytoMean.add((sumRed - sumColoc).toDouble())


            }
            /**Store table per slice (plane) */
            def tablePath = new File(outputDir, imps.getTitle().replaceAll(".tif", "") + "_" + "table_image_results" + ".csv").toString();
            IJ.log("Saving table per image: " + tablePath + " in " + outputDir.getAbsolutePath());
            tableSlice.save(tablePath);


            tableTotal.incrementCounter()
            /** Iterate per image */
            tableTotal.setValue("Image Title", i, imps.getTitle().replaceAll(".tif", "").toString())
            tableTotal.setValue("N Slices (Planes)", i, chBlue.stackSize.toString())
            tableTotal.setValue("Nuclei Mean Fluorescence", i, nucleiMean.stream()
                    .mapToDouble(a -> a)
                    .average().getAsDouble().toString())
            tableTotal.setValue("Nuclei Max Fluorescence", i, nucleiMax.stream()
                    .mapToDouble(a -> a)
                    .average().getAsDouble().toString())
            tableTotal.setValue("Nuclei Min Fluorescence", i, nucleiMin.stream()
                    .mapToDouble(a -> a)
                    .average().getAsDouble().toString())
            tableTotal.setValue("Nuclei Sum Fluorescence", i, nucleiSum.stream()
                    .mapToDouble(a -> a)
                    .average().getAsDouble().toString())
            tableTotal.setValue("Nuclei Area", i, nucleiArea.stream()
                    .mapToDouble(a -> a)
                    .average().getAsDouble().toString())
            tableTotal.setValue("Sum Red", i, sumRedTotal.stream()
                    .mapToDouble(a -> a)
                    .average().getAsDouble().toString())
            tableTotal.setValue("Cyto Mean Fluorescence", i, cytoMean.stream()
                    .mapToDouble(a -> a)
                    .average().getAsDouble().toString())


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


