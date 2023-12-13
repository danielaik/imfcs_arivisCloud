package com.apeer.modules;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Formatter;
import java.io.FileWriter;
//import java.awt.Rectangle;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.NonSquareMatrixException;
import org.apache.commons.math3.linear.NonSymmetricMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
//import org.json.simple.JSONObject;
//import org.json.simple.JSONArray;
import org.json.JSONArray;
import org.json.JSONObject;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.AbstractCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.special.Erf;

/**
 * autocorrelations of image stacks and data fitting for Fluorescence
 * Correlation Spectroscopy (FCS)
 * 
 * @author Daniel Aik
 * @version 1.0
 * @date 19/07/2020
 */

public class MultiTauCor {

    private boolean isGLS;
    private boolean isBayes;

    // settings
    private boolean isBcorrected; // is imagestack previously bleach corrected before performing ACF? 16 bit
                                  // implies no correction, 32 bit implies bleach corrected imagestack
    private int firstframe; // first frame to be used in the correlation
    private int lastframe; // last frame to be used in the correlation
    private double frametime; // acquisition time per frame
    private int frames; // number of frames in stack

    // parameters determined form the open image stack and defined by the user
    private int noSettings = 31;
    private String[] panelSettings = new String[noSettings]; // array to store the individual settings, the settings are
                                                             // stored in the same order as used to create the results
                                                             // table
    private boolean[] keyParam = { true, true, true, true, true, true, true, true, true, true, false, true, true, true,
            true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, false };
    private int correlatorp; // parameters for the correlator structure; correlatorp refers to the number of
                             // bins in the first channel group.
    private int correlatorq; // all higher groups have correlatorp/2 channels; correaltorq refers to the
                             // number of higher groups
    private double objmag; // microscope objective magnification; note the variable 'magnification' is used
                           // for the magnification of the stack image in ImageJ
    private double pixelsize; // pixel size in micrometer on camera chip (variable "a")
    private float background; // background is determined from the smallest count value from the stack; this
                              // can be maually changed in the control panel
    private int background2; // for FCCS background for the second region is determined as the minimum value
                             // for that region
    private int cfXshift; // x,y distance used in fitting of FCS correlations between two pixels; (0, 0)
                          // is for ACFs
    private int cfYshift;
    private boolean overlap = false; // can pixels overlap if binning is used?

    // parameters for cursor position in window for correlations
    private int pixbinX; // this is a multiplication factor to determine positions in the image window;
                         // it is 1 for overlap and equal to binning for non-overlap
    private int pixbinY;
    private int binningX;
    private int binningY;
    private int cfXDistance; // x,y distance between pixels to be correlated; in DC-FCCS this is the distance
                             // between corresponding pixels in red and green channels
    private int cfYDistance;
    private int pixelWidthX; // determines the number of pixels that can be correlated, depending on whether
                             // overlap is allowed for the width and height of the image
    private int pixelHeightY;
    private int roi1StartX = 0; // first pixel of the ROI used for correlation
    private int roi1StartY = 0;
    private int roi1WidthX = 0; // dimensions of the ROI used for correlations
    private int roi1HeightY = 0;

    // arrays and parameters used for computation of correlations and storage of
    // results
    private double[][] currentCovmats; // the regularized covariance matrix is at the moment not stored but only
                                       // updated for the current pixel if "GLS Fit" is selected
    private double[][][][] sdacf; // standerd deviation of the ACF; dimensions [ccf: ACF1, ACF2,
                                  // CCF][width][height][chanum]
    private double[][][][] acf; // ACF array to store all FCS functions for the pixels in the image; dimensions
                                // [ccf: ACF1, ACF2, CCF][width][height][chanum]
    private double[][][][] varacf; // standerd deviation of the ACF; dimensions [ccf: ACF1, ACF2,
                                   // CCF][width][height][chanum]
    private double[][][] blocked; // store yes/no whether blocking was succesful
    private float[][] filterArray; // filter array; intensity filters to select the pixels to be correlated
    private double[][] datac; // temporary array which stores the values along two pixel columns to be
                              // correlated through the stack
    private double[] intTrace1; // intensity traces
    private double[] intTime; // time for intensity traces
    private int[] lag; // lag for each correlation channel counted as multiples of the smallest basic
                       // time; independent of time; lagtime = lag * frametime
    private int base; // base = number of channels in first group
    private int hbase; // hbase = number of channels in all higher groups
    private double[] mtab; // number of samples for each correlation channel; dimensions [chanum]; defined
                           // in setParameters()
    private double[] samp; // sampletime (or bin width) of each channel; dimensions [chanum]; defined in
                           // setParameters()
    private int lagnum; // number of lag groups for the correlator
    private int chanum; // number of total channels of the correlator
    private double[] lagtime; // lagtime array which stores the correlation times; dimensions [chanum];
                              // defined in setParameters()
    private int nopit; // points in the shortened intensity traces for plotting
    private int blockIndS; // Index whether block Transformation is succesful (1) or maximum blocking is
                           // used (0)

    // Image window
    private int width; // width of loaded stack
    private int height; // height of loaded stack
    private float impmin; // minimum value in the stack
    private double maxsc; // scale for ACF plot
    private double minsc;
    private double imaxsc; // scale for intensity plot
    private double iminsc;// array with information whether a parameter is fit (true) or hold (false)

    // background image
    private boolean bgrloaded = false; // flag to indicate whether backgroudn image was loaded by user
    private double[][] bgrmean; // mean values (pixel, row, coumn) for background file

    // ImFCS fit panel
    private String fitModel;

    /*
     * //variable for fitting true (addition for fitting); previously program does
     * ACF only
     */
    // fitting mode
    private boolean doFit;
    private double[] paraminitval; // intitial values for the parameters used in the fit (can be the same as
                                   // initparam or can be continously set; both options are given but only the
                                   // first is used momentarily)
    private boolean[] paramfit; // array with information whether a parameter is fit (true) or hold (false)
    private double[] initparam;
    private int noparam = 11; // number of fit parameters; at the moment there are two fit models with a
                              // maximum of 11 parameters
    private double q2; // fit parameter, ratio of brightness of second component to first component (if
                       // it exists)
    private double q3; // fit parameter, ratio of brightness of third component to first component (if
                       // it exists)
    private double pixeldimx; // pixel size in object space in [m] used in fit
    private double pixeldimy; // pixel size in object space in [m] used in fit
    private double sigmaZ;
    private double fitobsvol; // normalization factor including the observation volume; for ACFs this is
                              // correct; note that N for spatial cross-correaltions has no physical meaning
    private double psfsize; // actual lateral PSF in [m] for laser 1, used in fit; this is the 1/e2 radius
    private double sigma;
    private double emlambda;
    private double NA; // numerical aperture
    private double lsthickness; // light sheet thickness for SPIM for laser 1 given as 1/e2 radius
    private int filterLL = 0; // intensity filter values with lower and upper limits LL and UL
    private int filterUL = 65536;

    private boolean[][][] pixfitted; // whether pixel has been successfully fitted or not; in the absence of
                                     // user-defined thresholds, this array determines pixvalid[][][]
    private double[][][][] fitres; // fit parameter results; [ccf: ACF1, ACF2, CCF][width][height][noparam]
    private double[][][] chi2; // chi2 values
    private String currentmodel;
    private final int fitMaxIterations = 2000; // maximum number of iterations allowed in the fits; the maximum number
                                               // is Integer.MAX_VALUE but that can lead to very long evaluations with
                                               // typically little improvement on the results
    private final int fitMaxEvaluations = 2000; // maximum number of evaluations allowed in the fits; the maximum number
                                                // is Integer.MAX_VALUE but that can lead to very long evaluations with
                                                // typically little improvement on the results
    private int fitstart = 1; // parameter for settign the range of points of the CF to be fitted
    private int fitend;
    private double[][][][] fitacf; // fit functions; [ccf: ACF1, ACF2, CCF][width][height][chanum]
    private double[][][][] res; // fit residuals; [ccf: ACF1, ACF2, CCF][width][height][chanum]
    private String $param[] = { "N", "D", "vx", "vy", "G", "F2", "D2", "F3", "D3", "Ftrip", "Ttrip", "reduced Chi2",
            "blocked", "valid pixels", "q map" };

    // GLS
    private RealMatrix lowerCholDecCovmats; // lower diagonal matrix from the Cholesky decomposition of the regularized
                                            // covariance matrix covamts
    private double[][] transTheoreticalGradientACF; // temporary array with transformed theoretical model function
                                                    // L\y(x)
    private double[] transTheoreticalACF; // temporary array with transformed theoretical model function L\y(x)
    private RealVector transACF; // transformed data L\y

    // Bayes
    private double modelprobP1; // model probability for 1 particle fits. will get updated once bayesian fit is
                                // done
    private double modelprobP2; // model probability for 2 particle fits

    // input tiff file
    public ImagePlus imp;
    public String $title;
    public String $filename;

    // output excel
    public File nf;

    // logger variable and methods
    public Formatter x;
    public FileWriter fw;
    public File testfile;

    public void openFile() {
        try {
            x = new Formatter("logger.txt");
        } catch (Exception e) {
            System.out.println("openFile error");
        }
    }

    public void addRecords(String inputImagepath, Double frametime, int binx, int biny, int p, int q, String abspath,
            String path) {
        x.format("%s\n", "this is a log file");
        x.format("\n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s", inputImagepath, frametime, binx, biny, p, q, abspath, path);
    }

    public void closeFile() {
        x.close();
    }

    public Map<String, Object> log(String inputImagePath, double Frametime, int BinX, int BinY, int p, int q) {
        var out = new HashMap<String, Object>();
        testfile = new File("test.csv");
        String b = testfile.getAbsolutePath();
        String c = testfile.getPath();
        openFile();
        addRecords(inputImagePath, Frametime, BinX, BinY, p, q, b, c);
        closeFile();

        out.put("output_file", "logger.txt");
        return out;

    }

    public Formatter x1;
    public JSONObject json_toHtml;

    public Map<String, Object> getHtml() throws IOException {

        var out = new HashMap<String, Object>();

        // to replace static content in template_v7.html
        String jstring = fillJson();
        String title = "FCS";

        // Appending html
        File htmlTemplateFile = new File("template_v9.html"); // alternatively "/usr/src/app/template/template_v7.html"
        boolean isExist = htmlTemplateFile.exists();
        String htmlString = FileUtils.readFileToString(htmlTemplateFile);
        htmlString = htmlString.replace("$title", title);
        htmlString = htmlString.replace("$json", jstring);

        // creating updated html to be exported
        String fileName = $title + ".html";
        File newHtmlFile = new File(fileName);
        FileUtils.writeStringToFile(newHtmlFile, htmlString);

        // logger start
        try {
            x1 = new Formatter("logger1.txt");
        } catch (Exception e) {
            System.out.println("openFile error");
        }

        x1.format("%s\n", "this is log file 2");
        x1.format("\n%s \n%s \n%s \n%s \n%s", isExist, htmlTemplateFile.getPath(), htmlTemplateFile.getAbsolutePath(),
                4, 4);
        x1.format("\n%s \n%s \n%s \n%s \n%s \n%s \n%s", title, jstring, htmlString, 4, 4, 4, 4);
        x1.close();
        // logger end

        // out.put("output_original", "logger1.txt");
        // out.put("output_original", htmlTemplateFile.getPath());
        out.put("output_new", newHtmlFile.getPath());

        return out;
    }

    /*
     * 
     * <!--
     * https://mvnrepository.com/artifact/com.googlecode.json-simple/json-simple -->
     * <dependency> <groupId>com.googlecode.json-simple</groupId>
     * <artifactId>json-simple</artifactId> <version>1.1.1</version> </dependency>
     * 
     * org.json.simple.JSONObject uses raw type collections internally. Alternative:
     * <dependency> <groupId>org.json</groupId> <artifactId>json</artifactId>
     * <version>20180813</version> </dependency>
     */

    private String fillJson() {

        json_toHtml = new JSONObject();
        json_toHtml.put("exposureTime", frametime);
        json_toHtml.put("correlatorP", correlatorp);
        json_toHtml.put("correlatorQ", correlatorq);
        json_toHtml.put("binningX", binningX);
        json_toHtml.put("binningY", binningY);
        json_toHtml.put("width", width);
        json_toHtml.put("height", height);
        json_toHtml.put("totalFrame", frames);
        json_toHtml.put("maxY", maxsc);
        json_toHtml.put("minY", minsc);
        json_toHtml.put("doFit", doFit);
        json_toHtml.put("objmag", objmag);
        json_toHtml.put("pixelsize", pixelsize);
        json_toHtml.put("NA", NA);
        json_toHtml.put("sigma", sigma);
        json_toHtml.put("emlambda", emlambda);
        json_toHtml.put("fitModel", fitModel);
        json_toHtml.put("isGLS", isGLS);
        json_toHtml.put("isBayes", isBayes);
        json_toHtml.put("title", $filename);
        json_toHtml.put("isBleachCorrected", isBcorrected);
        json_toHtml.put("imaxsc", imaxsc);
        json_toHtml.put("iminsc", iminsc);
        json_toHtml.put("filterLL", filterLL);
        json_toHtml.put("filterUL", filterUL);
        json_toHtml.put("background", background);

        // filling timelag
        JSONArray timeLag = new JSONArray(); // include correlation zero chanum = 73 not 72 for 16,8 correlator
        for (int i = 1; i < chanum; i++) {
            // timeLag.add(lagtime[i]);
            timeLag.put(lagtime[i]);

        }
        json_toHtml.put("timeLag", timeLag);

        // filling intTrace and intTime
        JSONArray jintTrace = new JSONArray();
        JSONArray jintTime = new JSONArray();
        for (int i = 0; i < nopit; i++) {
            jintTrace.put(intTrace1[i]);
            jintTime.put(intTime[i]);
            IJ.log("inside fill json: " + intTrace1[i] + " " + intTime[i]);
        }
        json_toHtml.put("intTrace", jintTrace);
        json_toHtml.put("intTime", jintTime);

        // filling G(tau) and sd
        JSONArray listOfAcf = new JSONArray();
        JSONArray listOfSd = new JSONArray();

        for (int y = 0; y < height / binningY; y++) { // convert 2D to 1D byrow
            for (int x = 0; x < width / binningX; x++) {
                JSONArray jacf = new JSONArray();
                JSONArray jsd = new JSONArray();
                for (int i = 1; i < chanum; i++) {
                    jacf.put(acf[0][x][y][i]);
                    jsd.put(sdacf[0][x][y][i]);
                }
                listOfAcf.put(jacf);
                listOfSd.put(jsd);
            }
        }
        json_toHtml.put("G(tau)", listOfAcf);
        json_toHtml.put("sd", listOfSd);

        // filling fitact if fitting is selected
        if (doFit) {
            // filling fitacf
            JSONArray listOfFitAcf = new JSONArray();
            for (int y = 0; y < height / binningY; y++) { // convert 2D to 1D byrow
                for (int x = 0; x < width / binningX; x++) {
                    JSONArray jfitacf = new JSONArray();
                    for (int i = 1; i < chanum; i++) {
                        jfitacf.put(fitacf[0][x][y][i]);
                    }
                    listOfFitAcf.put(jfitacf);
                }
            }
            json_toHtml.put("fitacf", listOfFitAcf);

            // filling D and isfitted
            JSONArray Dcoeff = new JSONArray(); // include correlation zero chanum = 73 not 72 for 16,8 correlator
            JSONArray jisFitted = new JSONArray();
            for (int h = 0; h < height / binningY; h++) {
                for (int w = 0; w < width / binningX; w++) {
                    if (!pixfitted[0][w][h]) { // set D to 0 if not fitted correctly
                        fitres[0][w][h][1] = 0.0;
                    }
                    Dcoeff.put(fitres[0][w][h][1]);
                    jisFitted.put(pixfitted[0][w][h]);
                }
            }
            json_toHtml.put("Dcoeff", Dcoeff);
            json_toHtml.put("isFitted", jisFitted);

        }

        // filling model probabiliy
        if (doFit && isBayes) {
            json_toHtml.put("modelprobP1", modelprobP1);
            json_toHtml.put("modelprobP2", modelprobP2);
        }

        return (json_toHtml.toString());
    }

    public Map<String, Object> run(String inputImagePath, double Frametime, int BinX, int BinY, int p, int q, int fit,
            double pixsize, double mag, double na, double inputEmLambda, int isGLSfit, int isBayesHypothesis,
            int filterLL, int filterUL) {

        var outputs = new HashMap<String, Object>();

        // reading tiff file
        imp = new ImagePlus(inputImagePath);
        $filename = imp.getTitle();
        $title = $filename.substring(0, $filename.lastIndexOf('.'));

        SetUI(Frametime, BinX, BinY, p, q, fit, pixsize, mag, na, inputEmLambda, isGLSfit, isBayesHypothesis, filterLL,
                filterUL);

        if (setParameters()) {

            int roi2StartX;
            int roi2WidthX;
            int roi2StartY;
            int roi2HeightY;
            if (cfXDistance > 0) {
                roi1StartX = 0;
                roi2StartX = cfXDistance;
            } else {
                roi1StartX = -cfXDistance;
                roi2StartX = 0;
            }
            if (cfYDistance > 0) {
                roi1StartY = 0;
                roi2StartY = cfYDistance;
            } else {
                roi1StartY = -cfYDistance;
                roi2StartY = 0;
            }
            if (overlap) {
                roi1WidthX = width - Math.abs(cfXDistance);
                roi1HeightY = height - Math.abs(cfYDistance);
            } else {
                roi1WidthX = (int) Math.floor((width - Math.abs(cfXDistance)) / binningX) * binningX;
                roi1HeightY = (int) Math.floor((height - Math.abs(cfYDistance)) / binningY) * binningY;
            }
            roi2WidthX = roi1WidthX;
            roi2HeightY = roi1HeightY;

            Roi impRoi1 = new Roi(roi1StartX, roi1StartY, roi1WidthX, roi1HeightY);
            correlateROI(impRoi1);
            plotCF(impRoi1, 2, false);
        }

        // creation of file
        String xlsxFN = $title + ".xlsx";
        nf = new File(xlsxFN);

        // writing file and closing
        writeExperimentACF(nf, "Failed to write data files", true);
        outputs.put("output_excel", xlsxFN);

        return outputs;

    }

    public void plotCF(Roi plotroi, int cormode, boolean map) {

        int ct = 0;
        int cpx1;
        int cpy1;
        int cpxf;
        int cpyf;
        int cpx2;
        int cpy2;

        /*
         * //unable to import Rectangle Rectangle plotrect = plotroi.getBounds(); cpx1 =
         * (int) Math.ceil(plotrect.getX() / pixbinX); cpy1 = (int)
         * Math.ceil(plotrect.getY() / pixbinY); cpxf = (int)
         * Math.floor((plotrect.getX() + plotrect.getWidth() - binningX) / pixbinX);
         * cpyf = (int) Math.floor((plotrect.getY() + plotrect.getHeight() - binningY) /
         * pixbinY);
         */

        cpx1 = (int) Math.ceil(0.0 / pixbinX);
        cpy1 = (int) Math.ceil(0.0 / pixbinY);
        cpxf = (int) Math.floor((0.0 + width - binningX) / pixbinX);
        cpyf = (int) Math.floor((0.0 + height - binningY) / pixbinY);
        cpx2 = (cpx1 * pixbinX + cfXDistance);
        cpy2 = (cpy1 * pixbinY + cfYDistance);

        maxsc = acf[ct][cpx1][cpy1][1]; // minimum and maximum setting for plot window
        minsc = 1000; // set minsc to a high value to make sure it is not 0

        if (cormode == 2) { // if multiple ACF or CCF
            for (int x = cpx1; x <= cpxf; x++) { // find minimum and maximum values in correlation functions that will
                                                 // be plotted
                for (int y = cpy1; y <= cpyf; y++) {
                    for (int z = 1; z <= (chanum - 1); z++) {
                        if ((!map && plotroi.contains(x * pixbinX, y * pixbinY)
                                && plotroi.contains(x * pixbinX, y * pixbinY + binningY - 1)
                                && plotroi.contains(x * pixbinX + binningX - 1, y * pixbinY)
                                && plotroi.contains(x * pixbinX + binningX - 1, y * pixbinY + binningY - 1))
                                || (map == true && plotroi.contains(x, y))) {
                            if (maxsc < acf[0][x][y][z]) {
                                maxsc = acf[0][x][y][z];
                            }
                            if (acf[0][x][y][z] != 0 && minsc > acf[0][x][y][z]) { // make sure minsc is not set to 0
                                                                                   // because of a missing CF
                                minsc = acf[0][x][y][z];
                            }
                        }
                    }
                }
            }

            // maximum scales are to be 10% larger than maximum value and 10% smaller than
            // minimum value
            minsc -= minsc * 0.1;
            maxsc += maxsc * 0.1;
        }

    }

    // Calculate correlations for a ROI
    public void correlateROI(Roi improi) {

        IJ.showStatus("Correlating ROI");

        // CPU calculation

        /*
         * //TODO: unable to import java.awt.* Rectangle imprect = improi.getBounds();
         * 
         * int startXmap = (int) Math.ceil(imprect.getX() / pixbinX); int startYmap =
         * (int) Math.ceil(imprect.getY() / pixbinY); int endXmap = (int)
         * Math.floor((imprect.getX() + imprect.getWidth() - binningX) / pixbinX); int
         * endYmap = (int) Math.floor((imprect.getY() + imprect.getHeight() - binningY)
         * / pixbinY);
         */

        int startXmap = (int) Math.ceil(0.0 / pixbinX);
        int startYmap = (int) Math.ceil(0.0 / pixbinY);
        int endXmap = (int) Math.floor((0.0 + width - binningX) / pixbinX);
        int endYmap = (int) Math.floor((0.0 + height - binningY) / pixbinY);

        int startX = startXmap * pixbinX;
        int startY = startYmap * pixbinY;
        int endX = endXmap * pixbinX;
        int endY = endYmap * pixbinY;
        int pixcount = 0;
        double q1;
        double q2;

        filterArray = new float[width][height]; // calculate the mean image of the stack

        for (int x1 = startX; x1 <= endX; x1 = x1 + pixbinX) {
            for (int x2 = startY; x2 <= endY; x2 = x2 + pixbinY) {
                for (int x3 = 0; x3 < binningX; x3++) {
                    for (int x4 = 0; x4 < binningY; x4++) {
                        if (improi.contains(x1, x2) && improi.contains(x1, x2 + binningY - 1)
                                && improi.contains(x1 + binningX - 1, x2)
                                && improi.contains(x1 + binningX - 1, x2 + binningY - 1)) {
                            filterArray[x1][x2] += imp.getStack().getProcessor(firstframe).getf(x1 + x3, x2 + x4);
                            pixcount++;
                        } else {
                            filterArray[x1][x2] = Float.NaN;
                        }
                    }
                }
            }
        }

        // do the FCS or DC-FCCS evaluation
        if (fitModel == "FCS") {
            if (doFit) {
                prepareFit();
                for (int x = startXmap; x <= endXmap; x++) {
                    for (int y = startYmap; y <= endYmap; y++) {
                        if (filterArray[x * pixbinX][y * pixbinY] >= filterLL * binningX * binningY
                                && filterArray[x * pixbinX][y * pixbinY] <= filterUL * binningX * binningY) {
                            correlate(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance,
                                    y * pixbinY + cfYDistance, 0, firstframe, lastframe);

                            if (!fit(x, y, 0, "FCS")) {
                                IJ.log("!fit = true");
                                return;
                            }

                        }
                    }
                    IJ.showProgress(x - startXmap, startXmap - endXmap);
                }
            } else {
                for (int x = startXmap; x <= endXmap; x++) {
                    for (int y = startYmap; y <= endYmap; y++) {
                        if (!Float.isNaN(filterArray[x * pixbinX][y * pixbinY])) {
                            correlate(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance,
                                    y * pixbinY + cfYDistance, 0, firstframe, lastframe);
                        }
                        IJ.showProgress(x - startXmap, startXmap - endXmap);
                    }
                }
            }
        }

        calcAverageIntensityTrace(filterArray, startX, startY, endX, endY, firstframe, lastframe);

    }

    // DAN: writing ACF1 into xls
    public void writeExperimentACF(File sfile, String $exception, boolean showlog) {
        int nov = width * height;
        // There is a limit of 1,048,576 rows and 16,384 (corresponds to 128 x 128)
        // columns in xlsx file as of 2019.
        if (nov > 16384 - 1) { // minus 1 because there is an additional (1st) column in Fit Parameters tab
                               // that list respective parameters
            IJ.log("Unable to save the data. Output data is larger than permissible 16,384 columns in .xlsx file.");
            return;
        }

        int t;
        File file;
        int nofpv = 35;
        String[] paramname = new String[nofpv];
        String[] paramsave = new String[nofpv];

        String $sfile = sfile.toString();
        int dotind = $sfile.lastIndexOf('.');
        if (dotind != -1) {
            $sfile = $sfile.substring(0, dotind);
        }
        file = new File($sfile + ".xlsx");
        XSSFWorkbook wb = new XSSFWorkbook();
        Sheet paraSheet = wb.createSheet("Panel Parameters");
        Sheet tauSheet = wb.createSheet("lagtime");
        Sheet acf0Sheet = wb.createSheet("ACF1");
        Sheet sd0Sheet = wb.createSheet("SD (ACF1)");
        Sheet fitacf0Sheet = wb.createSheet("Fit functions (ACF1)");
        Sheet fitres0Sheet = wb.createSheet("Fit Parameters (ACF1)");

        Row row;
        Row row0;
        // write Imaging FCS Panel parameters
        t = 0;
        paramname[t++] = "Image width";
        paramname[t++] = "Image height";
        paramname[t++] = "First frame";
        paramname[t++] = "Last frame";
        paramname[t++] = "Frame time";
        paramname[t++] = "Binning X";
        paramname[t++] = "Binning Y";
        paramname[t++] = "CF X distance";
        paramname[t++] = "CF Y distance";
        paramname[t++] = "Correlator P";
        paramname[t++] = "Correlator Q";
        paramname[t++] = "Fit model";
        paramname[t++] = "FCCS Display";
        paramname[t++] = "PixelSize";
        paramname[t++] = "Overlap";
        paramname[t++] = "Magnification";
        paramname[t++] = "NA";
        paramname[t++] = "Em Lambda";
        paramname[t++] = "Em Lambda2";
        paramname[t++] = "Sigma";
        paramname[t++] = "SigmaZ";
        paramname[t++] = "Sigma2";
        paramname[t++] = "SigmaZ2";
        paramname[t++] = "Background";
        paramname[t++] = "Background2";
        paramname[t++] = "Background File";
        paramname[t++] = "Bleach Correction";
        paramname[t++] = "Sliding Window Length";
        paramname[t++] = "Polynomial Order";
        paramname[t++] = "Filter";
        paramname[t++] = "filterUL";
        paramname[t++] = "filterLL";
        paramname[t++] = "Threshold";
        paramname[t++] = "Fit";

        t = 0;
        paramsave[t++] = Integer.toString(width);
        paramsave[t++] = Integer.toString(height);
        paramsave[t++] = Integer.toString(firstframe);
        paramsave[t++] = Integer.toString(lastframe);
        paramsave[t++] = Double.toString(frametime);
        paramsave[t++] = Integer.toString(binningX);
        paramsave[t++] = Integer.toString(binningY);
        paramsave[t++] = Integer.toString(cfXDistance);
        paramsave[t++] = Integer.toString(cfYDistance);
        paramsave[t++] = Integer.toString(correlatorp);
        paramsave[t++] = Integer.toString(correlatorq);
        paramsave[t++] = fitModel;
        paramsave[t++] = "Not Applicable";
        paramsave[t++] = Double.toString(pixelsize);
        paramsave[t++] = Boolean.toString(overlap);
        paramsave[t++] = Double.toString(objmag);
        paramsave[t++] = Double.toString(NA);
        paramsave[t++] = Double.toString(emlambda);
        paramsave[t++] = "Not Applicable";
        paramsave[t++] = Double.toString(sigma);
        paramsave[t++] = Double.toString(sigmaZ);
        paramsave[t++] = "Not Applicable";
        paramsave[t++] = "Not Applicable";
        paramsave[t++] = Float.toString(background);
        paramsave[t++] = "Not Applicable";
        paramsave[t++] = "Not Implemented";
        paramsave[t++] = Boolean.toString(isBcorrected);
        paramsave[t++] = "Not Applicable";
        String polyOrdeInfo = null;
        if (isBcorrected) {
            polyOrdeInfo = "Refer to parent Apeer Module";
        } else {
            polyOrdeInfo = "Not Applicable";
        }
        paramsave[t++] = polyOrdeInfo;
        paramsave[t++] = "Not Applicable";
        paramsave[t++] = Integer.toString(filterUL);
        paramsave[t++] = Integer.toString(filterLL);
        paramsave[t++] = "Not Applicable";
        String fitInfo = null;
        if (doFit) {
            if (isGLS) {
                if (isBayes) {
                    fitInfo = "Fit on, GLS true, Bayes true";
                } else {
                    fitInfo = "Fit on, GLS true, Bayes false";
                }
            } else {
                if (isBayes) {
                    fitInfo = "Fit on, GLS false, Bayes true";
                } else {
                    fitInfo = "Fit on, Standard least-square fitting";
                }
            }
        } else {
            fitInfo = "false";
        }
        paramsave[t++] = fitInfo;

        for (int i = 0; i < nofpv; i++) {
            row = paraSheet.createRow(i);
            row.createCell(0).setCellValue(paramname[i]);
            row.createCell(1).setCellValue(paramsave[i]);
        }

        // write tau
        row = tauSheet.createRow(0);
        row.createCell(0).setCellValue("S/N");
        row.createCell(1).setCellValue("Lagtime");
        row.createCell(2).setCellValue("Bin width");
        for (int i = 0; i < chanum; i++) {
            row = tauSheet.createRow(i + 1);
            row.createCell(0).setCellValue(i);
            row.createCell(1).setCellValue(lagtime[i]);
            row.createCell(2).setCellValue(samp[i]);
        }

        // write ACF, SD, fitacf, res
        Row rowacf0 = acf0Sheet.createRow(0);
        Row rowsd0 = sd0Sheet.createRow(0);
        Row rowfitacf0 = fitacf0Sheet.createRow(0);

        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                String $tmp = "(" + Integer.toString(w) + ", " + Integer.toString(h) + ")";
                rowacf0.createCell(w + h * width).setCellValue($tmp);
                rowsd0.createCell(w + h * width).setCellValue($tmp);
                rowfitacf0.createCell(w + h * width).setCellValue($tmp);
            }
        }

        for (int i = 0; i < chanum; i++) {
            acf0Sheet.createRow(i + 1);
            sd0Sheet.createRow(i + 1);
            fitacf0Sheet.createRow(i + 1);
        }

        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                for (int i = 0; i < chanum; i++) {
                    acf0Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(acf[0][w][h][i]);
                    sd0Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(sdacf[0][w][h][i]);
                    fitacf0Sheet.getRow(i + 1).createCell(w + h * width).setCellValue(fitacf[0][w][h][i]);
                }
            }
        }

        // add fit results
        row0 = fitres0Sheet.createRow(0);
        row0.createCell(0).setCellValue("Parameter");
        fitres0Sheet.createRow(1).createCell(0).setCellValue("fitted");
        for (int q = 0; q < noparam + 3; q++) {
            fitres0Sheet.createRow(q + 2).createCell(0).setCellValue($param[q]);
        }
        for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
                fitres0Sheet.getRow(0).createCell(w + h * width + 1)
                        .setCellValue("(" + Integer.toString(w) + ", " + Integer.toString(h) + ")");
                fitres0Sheet.getRow(1).createCell(w + h * width + 1).setCellValue(Boolean.toString(pixfitted[0][w][h]));
                fitres0Sheet.getRow(noparam + 2).createCell(w + h * width + 1).setCellValue(chi2[0][w][h]);
                fitres0Sheet.getRow(noparam + 3).createCell(w + h * width + 1).setCellValue(blocked[0][w][h]);
                fitres0Sheet.getRow(noparam + 4).createCell(w + h * width + 1).setCellValue(pixfitted[0][w][h]);
            }
        }
        for (int i = 0; i < noparam; i++) {
            for (int h = 0; h < height; h++) {
                for (int w = 0; w < width; w++) {
                    fitres0Sheet.getRow(i + 2).createCell(w + h * width + 1).setCellValue(fitres[0][w][h][i]);
                }
            }
        }

        // fitres0Sheet.autoSizeColumn(0); //java.lang.UnsatisfiedLinkError:
        // /usr/local/openjdk-11/lib/libfontmanager.so: libfreetype.so.6: cannot open
        // shared object file: No such file or directory

        // Write the output to a file
        try {
            FileOutputStream fileOut = new FileOutputStream(file);
            wb.write(fileOut);
            fileOut.close();
        } catch (IOException e) {
            throw new RuntimeException($exception, e);
        }

        if (showlog) {
            IJ.log("File(s) saved");
        }

    }

    // Set Param
    public boolean setParameters() {
        int index = 0;
        boolean resetResults = false; // whether Result arrays need to be reset
        boolean proceed = true; // whether program should proceed resetting the Results
        boolean onlySigmaOrBinChanged = true; // whether sigma0 is the only parameter changed in the panel - PSF
                                              // calibration is not reset in that case
        boolean onlyBinChanged = true; // whether binning is the only parameter changed in the panel - diff law is not
                                       // reset in that case
        String[] newPanelSettings = new String[noSettings]; // an array to temporarily hold the settings from the Panel

        // read settings from the panel and store them temporarily in newPanelSettings
        int tmpct = 0;
        newPanelSettings[tmpct++] = Integer.toString(firstframe); // index 0
        newPanelSettings[tmpct++] = Integer.toString(lastframe);
        newPanelSettings[tmpct++] = Double.toString(frametime);
        newPanelSettings[tmpct++] = Integer.toString(binningX);
        newPanelSettings[tmpct++] = Integer.toString(binningY);
        newPanelSettings[tmpct++] = "0";
        newPanelSettings[tmpct++] = "0";
        newPanelSettings[tmpct++] = Integer.toString(correlatorp);
        newPanelSettings[tmpct++] = Integer.toString(correlatorq);
        newPanelSettings[tmpct++] = "FCS";
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = Double.toString(pixelsize);
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = Double.toString(objmag);
        newPanelSettings[tmpct++] = Double.toString(NA);
        newPanelSettings[tmpct++] = Double.toString(emlambda);
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = Double.toString(sigma);
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = Float.toString(impmin);
        newPanelSettings[tmpct++] = "1";
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = "no bleach";
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = "4";
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = "";
        newPanelSettings[tmpct++] = "";

        if (proceed) {
            // assign the values to the variables used in the calculations

            try {
                firstframe = Integer.parseInt(newPanelSettings[0]);
                lastframe = Integer.parseInt(newPanelSettings[1]);
                frametime = Double.parseDouble(newPanelSettings[2]);
                binningX = Integer.parseInt(newPanelSettings[3]);
                binningY = Integer.parseInt(newPanelSettings[4]);
                cfXDistance = Integer.parseInt(newPanelSettings[5]);
                cfYDistance = Integer.parseInt(newPanelSettings[6]);
                correlatorp = Integer.parseInt(newPanelSettings[7]);
                correlatorq = Integer.parseInt(newPanelSettings[8]);
                pixelsize = Double.parseDouble(newPanelSettings[11]);
                objmag = Double.parseDouble(newPanelSettings[13]);
                NA = Double.parseDouble(newPanelSettings[14]);
                emlambda = Double.parseDouble(newPanelSettings[15]);
                sigma = Double.parseDouble(newPanelSettings[17]);
                background = Float.parseFloat(newPanelSettings[21]);
                background2 = Integer.parseInt(newPanelSettings[22]);
            } catch (NumberFormatException nfe) {
                IJ.showMessage("A parameter in the panel has an invalid format");
                throw new NumberFormatException("Number format error.");
            }

            if (fitModel == "FCS") {
                cfXshift = cfXDistance;
                cfYshift = cfYDistance;
            } else {
                cfXshift = 0;
                cfYshift = 0;
            }

            // set maximum, and minimum cursor positions possible in the image, depending on
            // binning and whether overlap is allowed
            // parameters need to be checked according to these settings as allowed
            // parameter ranges differ in overlap and non-overlap mode
            if (overlap) {
                pixelWidthX = width - binningX; // these values are the correct maximum if counting from 0
                pixelHeightY = height - binningY;
                pixbinX = 1;
                pixbinY = 1;
            } else {
                pixelWidthX = (int) Math.floor(width / binningX) - 1;
                pixelHeightY = (int) Math.floor(height / binningY) - 1;
                pixbinX = binningX;
                pixbinY = binningY;
            }

            // set common arrays and parameters according to user settings in the panel
            if ((lastframe - firstframe + 1) < 1000) { // use 1000 points for the intensity, except when less than 1000
                                                       // frames are present
                nopit = (lastframe - firstframe + 1);
            } else {
                nopit = 1000;
            }

            // if sliding window correction is needed then the correlator structure has to
            // be adapted to the smaller number of lagtimes
            int num; // total number of frames to be correlated; sliding window length or
                     // lastframe-firstframe+1

            chanum = (int) (correlatorp + (correlatorq - 1) * correlatorp / 2 + 1);
            lagnum = (int) correlatorq;
            num = (lastframe - firstframe + 1);

            // initialize arrays required for calculations; they change with each new
            // paramter setting
            base = (int) correlatorp; // base = number of channels in first group
            hbase = (int) correlatorp / 2; // hbase = number of channels in all higher groups
            mtab = new double[chanum]; // number of samples for each correlation channel
            lag = new int[chanum]; // lag for each correlation channel; indepenednet of time
            samp = new double[chanum]; // sampletime (or bin width) of each channel
            lagtime = new double[chanum]; // lagtime = lag*frametime; this is the actual lagtime in seconds for each
                                          // channel

            for (int x = 0; x <= hbase; x++) { // calculate lag and lagtimes for the 0 lagtime channel and the first 8
                                               // channels
                lag[x] = x;
                lagtime[x] = lag[x] * frametime;
            }

            for (int x = 1; x <= lagnum; x++) { // calculate lag and lagtimes for all higher channels
                for (int y = 1; y <= hbase; y++) {
                    lag[x * hbase + y] = (int) (Math.pow(2, x - 1) * y + (base / 4) * Math.pow(2, x));
                    lagtime[x * hbase + y] = lag[x * hbase + y] * frametime;
                }
            }

            for (int x = 0; x <= base; x++) { // calculate sampletimes (bin width) for the 0 lagtime channel and the
                                              // first 8 channels
                samp[x] = 1;
            }

            for (int x = 2; x <= lagnum; x++) { // calculate sampletimes (bin width) for all higher channels
                for (int y = 1; y <= hbase; y++) {
                    samp[x * hbase + y] = Math.pow(2, x - 1);
                }
            }

            // calculate the number of samples for each channel including 0 lagtime; this
            // differs for sliding window correction
            // the variable num takes care of this
            for (int x = 0; x <= (chanum - 1); x++) {
                mtab[x] = (int) Math.floor((num - lag[x]) / samp[x]);
            }

            // set parameter in fit windows...
            pixeldimx = (pixelsize * 1000 / objmag * binningX) / Math.pow(10, 9);
            pixeldimy = (pixelsize * 1000 / objmag * binningY) / Math.pow(10, 9);
            psfsize = (sigma * emlambda / NA) / Math.pow(10, 9);
            lsthickness = (sigmaZ * emlambda / NA) / Math.pow(10, 9);

            return true;
        } else {
            return false;
        }
    }

    // setting input UI
    public void SetUI(double ft, int bx, int by, int cp, int cq, int ifit, double ipixsize, double imag, double ina,
            double iemlambda, int iGLS, int iBayes, int filter_LL, int filter_UL) {

        // get width, height, number of frames in stack, and magnification of the image
        // window
        width = imp.getWidth();
        height = imp.getHeight();
        frames = imp.getStackSize();
        firstframe = 1;
        lastframe = imp.getImageStackSize();
        frametime = ft; // default=0.001
        binningX = bx; // default=1
        binningY = by; // default=1
        correlatorp = cp; // 16 or 32 ; default=16
        correlatorq = cq; // default=8

        if (ifit == 0) {
            doFit = false;
        } else {
            doFit = true;
        }

        // doFit = ifit;

        if (iGLS == 0) {
            isGLS = false;
        } else {
            isGLS = true;
        }

        if (iBayes == 0) {
            isBayes = false;
        } else {
            isBayes = true;
        }

        pixelsize = ipixsize;// in um
        objmag = imag;
        sigmaZ = 1000000;
        sigma = 0.8;
        emlambda = iemlambda;
        fitModel = "FCS";
        NA = ina;
        q2 = 1;// TODO set default to 1. briughtness ratio
        q3 = 1; // TODO set default to 1
        filterLL = filter_LL;
        filterUL = filter_UL;

        if (imp.getBitDepth() == 16) {
            isBcorrected = false;
        } else {
            isBcorrected = true; // input image came from Apeer bleach correction module
        }

        // setting background: bg=0 if 32bit(input file came from bleach corrected
        // module, alread subtracte by impmin), bg=impmin if 16bit(input file does not
        // came from bleach correction module)
        if (isBcorrected) {
            impmin = 0;// minDetermination(imp);
        } else {
            impmin = 0;// minDetermination(imp); // calculate the minimum of the image stack; this will
                       // be used as the
                       // default background value
        }

        chanum = (int) (correlatorp + (correlatorq - 1) * correlatorp / 2 + 1);
        lagnum = (int) correlatorq;

        initializeArrays();

        // initialize arrays for fitting
        paraminitval = new double[noparam];
        paramfit = new boolean[noparam];

        // initial paramter settings
        initparam = new double[noparam];
        initparam[0] = 1.0; // remember starting values
        initparam[1] = 1.0 / Math.pow(10, 12);
        initparam[2] = 0 / Math.pow(10, 6);
        initparam[3] = 0 / Math.pow(10, 6);
        initparam[4] = 0.0;
        initparam[5] = 0;
        initparam[6] = 0;
        initparam[7] = 0;
        initparam[8] = 0;
        initparam[9] = 0;
        initparam[10] = 0;

    }

    // determine minimum value in stack
    public float minDetermination(ImagePlus image) {
        float min;
        min = image.getStack().getProcessor(1).getf(1, 1);
        for (int z = 1; z <= frames; z++) {
            for (int x = 1; x < width; x++) {
                for (int y = 1; y < height; y++) {
                    if (image.getStack().getProcessor(z).getf(x, y) < min) {
                        min = image.getStack().getProcessor(z).getf(x, y);
                    }
                }
            }
        }

        return (float) Math.floor(min);
    }

    public void initializeArrays() {
        acf = new double[3][width][height][chanum];
        varacf = new double[3][width][height][chanum];
        sdacf = new double[3][width][height][chanum];
        fitacf = new double[3][width][height][chanum];
        res = new double[3][width][height][chanum];
        blocked = new double[3][width][height];

        initializeFitres(3, width, height, noparam); // initialize fitres and pixfiltered
        storeTempFit.init(noparam); // for bayesian temp storage
    }

    // calcualte the average intensity for the image
    public void calcAverageIntensityTrace(float[][] filterArray, int startX, int startY, int endX, int endY,
            int initialframe, int finalframe) {
        // introi: roi over which average is to be determined
        // initialframe and finalframe provide the range of frames to be used

        int ave;
        int pixcount = 0;
        ave = (int) Math.floor((finalframe - initialframe + 1) / nopit);
        intTrace1 = new double[nopit]; // reset before updating.
        intTime = new double[nopit];

        for (int i = 0; i < nopit; i++) {
            for (int j = firstframe + i * ave; j <= firstframe + (i + 1) * ave - 1; j++) {
                pixcount = 0;
                for (int x1 = startX; x1 <= endX; x1 = x1 + pixbinX) {
                    for (int x2 = startY; x2 <= endY; x2 = x2 + pixbinY) {
                        for (int x4 = 0; x4 < binningX; x4++) {
                            for (int x5 = 0; x5 < binningY; x5++) {
                                if (!Float.isNaN(filterArray[x1][x2])) {
                                    intTrace1[i] += imp.getStack().getProcessor(j).getf(x1 + x4, x2 + x5) - background;
                                    pixcount++;
                                    intTime[i] = frametime * (i + 0.5) * ave;
                                }
                            }
                        }
                    }
                }
            }
            intTrace1[i] /= (ave * pixcount); // calculate average intensity for the 'ave' points
        }

        iminsc = intTrace1[1]; // minimum and maximum setting for plot window
        imaxsc = intTrace1[1];

        for (int x = 0; x < nopit; x++) { // find minimum and maximum values in intensity trace 1 that will be plotted
            if (imaxsc < intTrace1[x]) {
                imaxsc = intTrace1[x];
            }
            if (iminsc > intTrace1[x]) {
                iminsc = intTrace1[x];
            }
        }

        iminsc -= iminsc * 0.1; // maximum scales are to be 10% larger than maximum value and 10% smaller than
                                // minimum value
        imaxsc += imaxsc * 0.1;
    }

    // correlate one pixel with itself or two pixels with each other
    public void correlate(ImagePlus image, int px1, int py1, int px2, int py2, int kcf, int initialframe,
            int finalframe) {
        // image: the imp to be used
        // px1, py1, px2, py2: pixel cooredinates for pixel 1 and pixel 2 which are to
        // be correalted
        // if px1 = px2 AND py1 = py2: then a autocorrelation is calculated
        // kcf (0, 1, or 2) determines whether ACF1, ACF2, or CCF is calculated
        // initialframe and finalframe provide the range of frames to be used for the
        // correlation
        int num = (finalframe - initialframe + 1); // total number of frames to be correlated
        int numofsw; // number of sliding windows
        int swinitialframe; // if sliding window (bleach) correction is selected these are the initial and
                            // final frames of the sub-windows
        int swfinalframe;
        int pxm1; // pixel coordinates on the binned grid used to store the output and map it to
                  // the parameter map
        int pym1;

        if (kcf == 1 && px1 == px2 && py1 == py2) { // if red channel in the DC-FCCS mode, map the output to the
                                                    // corresponding green-channel pixels on the binned grid
            pxm1 = (int) (px1 - cfXDistance) / pixbinX;
            pym1 = (int) (py1 - cfYDistance) / pixbinY;
        } else { // otherwise map to the pixel on the pixel on the binned grid
            pxm1 = (int) px1 / pixbinX;
            pym1 = (int) py1 / pixbinY;
        }

        // Sliding window is not selected, correlate the full intensity trace
        datac = new double[2][num + 1]; // get the intensity data for the correlation
        datac[0] = getIntensity(image, px1, py1, 1, initialframe, finalframe); // getIntensity for first pixel; performs
                                                                               // a bleach correction if indicated in
                                                                               // the panel
        datac[1] = datac[0]; // perform an autocorrelation

        Map result;
        result = correlator(datac, initialframe, finalframe); // correlate the data
        acf[kcf][pxm1][pym1] = (double[]) result.get("corav"); // acf
        varacf[kcf][pxm1][pym1] = (double[]) result.get("blockvar"); // variance of the ACF; blocked
        sdacf[kcf][pxm1][pym1] = (double[]) result.get("blocksd"); // standard deviation of the ACF; blocked
        currentCovmats = (double[][]) result.get("covmats");

        blocked[kcf][pxm1][pym1] = blockIndS; // store whether blocking worked successfully for the pixel

    }

    // correlator calculates correlation functions
    public Map correlator(double[][] intcor, int initialframe, int finalframe) {
        // intcor contains the array of intensity values to be correlated for pixels 1
        // and 2
        // initialframe and finalframe provide the range of frames to be used for the
        // correlation
        int num = (finalframe - initialframe + 1); // total number of frames to be correlated
        int blockIndex; // index at which optimal blocking is reached; if it fails maximal blocking is
                        // used

        blockIndex = blockTransform(intcor, num, 1); // perform blocking on the first channel to determine when
                                                     // intensity bins are independent

        Map result;
        result = calculateCF(intcor, num, blockIndex, 1); // perform optimal blocking and return the CF, SD and
                                                          // covariance matrix

        return result;
    }

    // calculate the standard deviation by blocking
    public Map calculateCF(double[][] intcor, int num, int ind, int blocklag) {
        // intcor is the array of intensity values for the two traces which are
        // correlated
        // num is the number of frames which are correlated
        // ind is the blockindex, at which the SD has converged, previously found in
        // blockSD()
        // blocklag defines for which lag the blocking will be done; typically we use
        // the smalles, i.e. 1
        int numbin = num; // number of data points when they are binned
        int del; // delay or correlation time expressed in lags
        int currentIncrement;
        int ctbin = 0; // count how often the data was binned
        int binct;
        int pnum = (int) Math
                .floor(mtab[chanum - 1] / Math.pow(2, Math.max(ind - Math.log(samp[chanum - 1]) / Math.log(2), 0))); // minimum
                                                                                                                     // number
                                                                                                                     // of
                                                                                                                     // prodcuts
                                                                                                                     // given
                                                                                                                     // the
                                                                                                                     // used
                                                                                                                     // correlator
                                                                                                                     // structure
                                                                                                                     // and
                                                                                                                     // the
                                                                                                                     // blockIndex
                                                                                                                     // ind
        int[] prodnum = new int[chanum];
        double sumprod; // sum of all intensity products; divide by num to get the average
                        // <i(n)i(n+del)>
        double sumprod2; // sum of all intensity products squared; divide by num to get the average
                         // <(i(n)i(n+del))^2>
        double[] directm = new double[chanum]; // direct monitor required for ACF normalization
        double[] delayedm = new double[chanum]; // delayed monitor required for ACF normalization
        double[] blockvar;
        double[] blocksd;
        double[][] intblock;
        double[][] prod = new double[chanum][num];
        double[] corav = new double[chanum];
        double[] mcov = new double[chanum];
        double[] diagcovmat = new double[chanum];
        double[][] covmat = new double[chanum][chanum];
        double[][] covmats = new double[chanum - 1][chanum - 1]; // the final results does not contain information about
                                                                 // the zero lagtime channel
        double[][] cormat = new double[chanum][chanum];
        double[][] denomshrink = new double[chanum][chanum];
        double numerator = 0;
        double denominator = 0;
        double median = 0;
        double lamvar;
        double lamcov;

        intblock = new double[2][num];
        blockvar = new double[chanum];
        blocksd = new double[chanum];

        for (int x = 0; x < numbin; x++) { // re-initialize the intensities
            intblock[0][x] = intcor[0][x + 1];
            intblock[1][x] = intcor[1][x + 1];
        }

        currentIncrement = 1; // at the moment we always do blocking for smallest lag
        blocksd[0] = 0; // we do not calcualte the SD for the 0 lagtime as it is not used for fitting
                        // (shot noise)

        for (int x = 0; x < chanum; x++) { // run over all channels except the 0 lag time

            if (currentIncrement != samp[x]) { // check whether the channel width has changed
                numbin = (int) Math.floor(numbin / 2); // if yes, correct the number of actual data points
                currentIncrement = (int) samp[x]; // set the currentIncrement accordingly
                ctbin++; // count how often the data was binned
                for (int y = 0; y < numbin; y++) { // and bin the data according to the width of the current channel
                    intblock[0][y] = (intblock[0][2 * y] + intblock[0][2 * y + 1]);
                    intblock[1][y] = (intblock[1][2 * y] + intblock[1][2 * y + 1]);
                }
            }

            del = lag[x] / currentIncrement; // calculate the delay, i.e. the correlation time ...
            prodnum[x] = numbin - del; // and the number of products for that delay; //(int)
                                       // (mtab[chanum-1]*(samp[chanum-1]/samp[x]));//IJ.log(Double.toString(prodnum[x]));
            for (int y = 0; y < prodnum[x]; y++) { // calculate the ...
                directm[x] += intblock[0][y]; // direct and ...
                delayedm[x] += intblock[1][y + del]; // delayed monitor
            }
            directm[x] /= prodnum[x]; // calculate average of direct and delayed monitor, i.e. the average intensity
                                      // <n(0)> and <n(tau)>
            delayedm[x] /= prodnum[x];

            sumprod = 0;
            sumprod2 = 0;

            for (int y = 0; y < prodnum[x]; y++) { // calculate the correlation
                prod[x][y] = intblock[0][y] * intblock[1][y + del] - delayedm[x] * intblock[0][y]
                        - directm[x] * intblock[1][y + del] + delayedm[x] * directm[x];
                sumprod += prod[x][y]; // calculate the sum of prod, i.e. the raw correlation value ...
                sumprod2 += Math.pow(prod[x][y], 2); // ... and the sum of the squares
            }

            corav[x] = sumprod / (prodnum[x] * directm[x] * delayedm[x]); // calculate the ACF, i.e. the mean for the
                                                                          // later calculations of the
                                                                          // variance-covariance matrix

            binct = ind - ctbin; // determine whether data needs to be further binned or is already exceeding the
                                 // blocking number
            sumprod = 0;
            sumprod2 = 0;

            for (int y = 1; y <= binct; y++) { // bin the data until block time is reached
                prodnum[x] = (int) Math.floor(prodnum[x] / 2); // for each binning the number of data points is halfed
                for (int z = 0; z < prodnum[x]; z++) { // do the binning and divide by 2 so that average value does not
                                                       // change
                    prod[x][z] = (prod[x][2 * z] + prod[x][2 * z + 1]) / 2;
                }
            }

            prodnum[x] = pnum; // use only the minimal number of products to achieve a symmetric variance
                               // matrix
            for (int z = 0; z < prodnum[x]; z++) {
                sumprod += prod[x][z]; // calculate the sum of prod, i.e. the raw correlation value ...
                sumprod2 += Math.pow(prod[x][z], 2); // ... and the sum of the squares
            }

            blockvar[x] = (sumprod2 / prodnum[x] - Math.pow(sumprod / prodnum[x], 2))
                    / ((prodnum[x] - 1) * Math.pow(directm[x] * delayedm[x], 2)); // variance after blocking; extra
                                                                                  // division by prodnum to obtain SEM
            blocksd[x] = Math.sqrt(blockvar[x]); // standard deviation after blocking
        }

        // if GLS is selected then calulate the regularized covariance matrix
        if (isGLS) {
            // Calculate the mean of the products used for the variance-covariance matrix
            for (int x = 1; x < chanum; x++) {
                for (int z = 0; z < pnum; z++) {
                    mcov[x] += prod[x][z] / (directm[x] * delayedm[x]);
                }
                mcov[x] /= pnum; // normalize by the number of products
            }

            // Calculate the variance-covariance matrix
            for (int x = 1; x < chanum; x++) {
                for (int y = 1; y <= x; y++) { // calculate only the upper triangular part as the matrix is symmetric
                    for (int z = 0; z < pnum; z++) {
                        covmat[x][y] += (prod[x][z] / (directm[x] * delayedm[x]) - mcov[x])
                                * (prod[y][z] / (directm[y] * delayedm[y]) - mcov[y]);
                    }
                    covmat[x][y] /= (pnum - 1); // normalize by the number of products
                    covmat[y][x] = covmat[x][y]; // lower triangular part is equal to upper triangular part
                }
            }

            // Regularize variance-covariance matrix
            // first determine the shrinkage weight for the variance
            for (int x = 1; x < chanum; x++) { // get the variance (diagonal of covariance matrix) ...
                diagcovmat[x] = covmat[x][x];
            }

            Arrays.sort(diagcovmat); // ... and determine the median
            double pos1 = Math.floor((diagcovmat.length - 1.0) / 2.0);
            double pos2 = Math.ceil((diagcovmat.length - 1.0) / 2.0);
            if (pos1 == pos2) {
                median = diagcovmat[(int) pos1];
            } else {
                median = (diagcovmat[(int) pos1] + diagcovmat[(int) pos2]) / 2.0;
            }

            double tmpnum; // determine the variance of the variance
            for (int x = 1; x < chanum; x++) {
                tmpnum = 0;
                for (int z = 0; z < pnum; z++) {
                    tmpnum += Math.pow((Math.pow(prod[x][z] / (directm[x] * delayedm[x]) - mcov[x], 2) - covmat[x][x]),
                            2);
                }
                tmpnum *= (pnum) / Math.pow(pnum - 1, 3);
                numerator += tmpnum;
                denominator += Math.pow(covmat[x][x] - median, 2);
            }
            lamvar = Math.min(1, numerator / denominator); // shrinkage weight for the variance
            lamvar = Math.max(lamvar, 0);

            // determine the shrinkage weight for the covariance
            for (int x = 1; x < chanum; x++) { // calculate the sample correlation matrix
                for (int y = 1; y < chanum; y++) {
                    cormat[x][y] = covmat[x][y] / Math.sqrt(covmat[x][x] * covmat[y][y]);
                }
            }

            numerator = 0;
            denominator = 0;
            double cmx; // tmp variables to simplify ...
            double cmy; // ... in the loop
            for (int x = 1; x < chanum; x++) { // determine the variance of the covariance
                tmpnum = 0;
                for (int y = 1; y < x; y++) { // sum only over the upper triangle as the matrix is symmetric
                    for (int z = 0; z < pnum; z++) {
                        cmx = (prod[x][z] / (directm[x] * delayedm[x]) - mcov[x]) / Math.sqrt(covmat[x][x]);
                        cmy = (prod[y][z] / (directm[y] * delayedm[y]) - mcov[y]) / Math.sqrt(covmat[y][y]);
                        tmpnum += Math.pow(cmx * cmy - cormat[x][y], 2);
                    }
                    tmpnum *= (pnum) / Math.pow(pnum - 1, 3);
                    numerator += tmpnum;
                    denominator += Math.pow(cormat[x][y], 2); // sum of squares of off-diagonal elements of correlation
                                                              // matrix
                }
            }
            lamcov = Math.min(1, numerator / denominator); // shrinkage weight for the covariance
            lamcov = Math.max(lamcov, 0);

            // calculate the off-diagonal elements of the regularized variance-covariance
            // matrix
            for (int x = 1; x < chanum; x++) { // do not include zero lagtime channel as we don't use it for fitting
                for (int y = 1; y < x; y++) {
                    cmx = lamvar * median + (1 - lamvar) * covmat[x][x];
                    cmy = lamvar * median + (1 - lamvar) * covmat[y][y];
                    covmats[x - 1][y - 1] = (1 - lamcov) * cormat[x][y] * Math.sqrt(cmx * cmy) / pnum;
                    covmats[y - 1][x - 1] = covmats[x - 1][y - 1];
                }
            }
            for (int x = 1; x < chanum; x++) { // diagonal elements of the regularized variance-covariance matrix
                covmats[x - 1][x - 1] = (lamvar * median + (1 - lamvar) * covmat[x][x]) / pnum;
            }
        } // end of if GLS selected statement

        Map<String, Object> map = new HashMap<>();
        if (isGLS) { // hand over either the correlation function corav or the actual function used
                     // to calcualte the covariance matrix; they differ only slightly
            map.put("corav", mcov);
        } else {
            map.put("corav", corav);
        }
        map.put("blockvar", blockvar);
        map.put("blocksd", blocksd);
        map.put("covmats", covmats);
        return map;

    }

    public int blockTransform(double[][] intcor, int num, int blocklag) {
        // intcor is the array of intensity values for the two traces which are
        // correlated
        // num is the number of frames which are correlated
        // blocklag defines for which lag the blocking will be done; typically we use
        // the smalles, i.e. 1
        int blocknum = (int) Math.floor(Math.log(mtab[blocklag]) / Math.log(2)) - 2; // number of blocking operations
                                                                                     // that can be performed given
                                                                                     // blocklag
        int numbin = num; // number of data points when they are binned
        int del; // delay or correlation time expressed in lags
        int currentIncrement;
        int crwin = 2; // 3 points that fit the error bar overlap criterion
        double sumprod = 0.0; // sum of all intensity products; divide by num to get the average
                              // <i(n)i(n+del)>
        double sumprod2 = 0.0; // sum of all intensity products squared; divide by num to get the average
                               // <(i(n)i(n+del))^2>
        double directm = 0.0; // direct monitor required for ACF normalization
        double delayedm = 0.0; // delayed monitor required for ACF normalization
        double[][] intblock;
        double[] prod = new double[num];
        double[][] varblock;
        double[] upper;
        double[] lower;
        double[][] blockpoints;
        double[] blocksd;
        int[] crt;
        int[] cr12; // do neighbouring points have overlapping error bars; together with crwin=2
                    // this tests for three points that have overlapping erro bars
        int[] cr3;
        int[] diffpos;
        int last0 = 0;
        int ind = 0;
        double[] prodnum;
        double minblock;
        double maxblock;

        varblock = new double[3][blocknum];
        prodnum = new double[blocknum];
        intblock = new double[2][num];
        blocksd = new double[chanum];
        upper = new double[blocknum];
        lower = new double[blocknum];
        crt = new int[blocknum - 1];
        cr12 = new int[blocknum - 2];
        cr3 = new int[blocknum - 2];
        diffpos = new int[blocknum - 1];
        blockpoints = new double[3][3];

        for (int x = 0; x < numbin; x++) {
            intblock[0][x] = intcor[0][x];
            intblock[1][x] = intcor[1][x];
        }

        currentIncrement = blocklag; // at the moment we always do blocking for smallest lag which is 1 but in
                                     // general it can be used freely

        for (int x = 1; x < chanum; x++) { // run over all channels
            if (currentIncrement != samp[x]) { // check whether the channel width has changed
                currentIncrement = (int) samp[x]; // set the currentIncrement accordingly
                numbin = (int) Math.floor(numbin / 2); // and correct the number of actual data points accordingly
                for (int y = 0; y < numbin; y++) { // if yes, bin the data according to the width of the current channel
                    intblock[0][y] = (intblock[0][2 * y] + intblock[0][2 * y + 1]);
                    intblock[1][y] = (intblock[1][2 * y] + intblock[1][2 * y + 1]);
                }

            }

            if (x == blocklag) { // if the channel number is equal to the blocklag ...
                del = lag[x] / currentIncrement; // calculate the delay, i.e. the correlation time
                for (int y = 0; y < numbin - del; y++) { // calculate the ...
                    directm += intblock[0][y]; // direct and ...
                    delayedm += intblock[1][y + del]; // delayed monitor
                }
                prodnum[0] = numbin - del; // number of correlation products
                directm /= prodnum[0]; // calculate average of direct and delayed monitor,
                delayedm /= prodnum[0]; // i.e. the average intesity <n(0)> and <n(tau)>

                for (int y = 0; y < prodnum[0]; y++) { // calculate the correlation
                    prod[y] = intblock[0][y] * intblock[1][y + del] - delayedm * intblock[0][y]
                            - directm * intblock[1][y + del] + delayedm * directm;
                    sumprod += prod[y]; // calculate the sum of prod, i.e. the raw correlation value ...
                    sumprod2 += Math.pow(prod[y], 2); // ... and the sum of the squares
                }

                varblock[0][0] = currentIncrement * frametime; // the time of the block curve
                varblock[1][0] = (sumprod2 / prodnum[0] - Math.pow(sumprod / prodnum[0], 2))
                        / (prodnum[0] * Math.pow(directm * delayedm, 2)); // value of the block curve

                for (int y = 1; y < blocknum; y++) { // perform blocking operations
                    prodnum[y] = (int) Math.floor(prodnum[y - 1] / 2); // the number of samples for the blocking curve
                                                                       // decreases by a factor 2 with every step
                    sumprod = 0;
                    sumprod2 = 0;
                    for (int z = 0; z < prodnum[y]; z++) { // bin the correlation data and calculate the blocked values
                                                           // for the SD
                        prod[z] = (prod[2 * z] + prod[2 * z + 1]) / 2;
                        sumprod += prod[z];
                        sumprod2 += Math.pow(prod[z], 2);
                    }
                    varblock[0][y] = (currentIncrement * Math.pow(2, y)) * frametime; // the time of the block curve
                    varblock[1][y] = (sumprod2 / prodnum[y] - Math.pow(sumprod / prodnum[y], 2))
                            / (prodnum[y] * Math.pow(directm * delayedm, 2)); // value of the block curve
                }
            }
        }

        for (int x = 0; x < blocknum; x++) {
            varblock[1][x] = Math.sqrt(varblock[1][x]); // calculate the standard deviation
            varblock[2][x] = varblock[1][x] / Math.sqrt(2 * (prodnum[x] - 1)); // calculate the error
            upper[x] = varblock[1][x] + varblock[2][x]; // upper and lower quartile
            lower[x] = varblock[1][x] - varblock[2][x];
        }

        // determine index where blocking criteria are fulfilled
        for (int x = 0; x < blocknum - 1; x++) { // do neighboring points have overlapping error bars?
            if (upper[x] > lower[x + 1] && upper[x + 1] > lower[x]) {
                crt[x] = 1;
            }
        }

        for (int x = 0; x < blocknum - 2; x++) { // do three adjacent points have overlapping error bars?
            if (crt[x] * crt[x + 1] == 1) {
                cr12[x] = 1;
            }
        }

        for (int x = 0; x < blocknum - 1; x++) { // do neighboring points have a positive difference (increasing SD)?
            if (varblock[1][x + 1] - varblock[1][x] > 0) {
                diffpos[x] = 1;
            }
        }

        for (int x = 0; x < blocknum - 2; x++) { // do three neighboring points monotonically increase?
            if (diffpos[x] * diffpos[x + 1] == 1) {
                cr3[x] = 1;
            }
        }

        for (int x = 0; x < blocknum - 2; x++) { // find the last triple of points with monotonically increasing
                                                 // differences and non-overlapping error bars
            if ((cr3[x] == 1 && cr12[x] == 0)) {
                last0 = x;
            }
        }

        for (int x = 0; x <= last0; x++) { // indices of two pairs that pass criterion 1 an 2
            cr12[x] = 0;
        }

        cr12[blocknum - 3] = 0; // criterion 3, the last two points can't be part of the blocking triple
        cr12[blocknum - 4] = 0;

        for (int x = blocknum - 5; x > 0; x--) { // index of triplet with overlapping error bars and after which no
                                                 // other triplet has a significant monotonic increase
            if (cr12[x] == 1) { // or 4 increasing points
                ind = x + 1; // take the middle of the three points as the blocking limit
            }
        }

        if (ind == 0) { // if optimal blocking is not possible, use maximal blocking
            blockIndS = 0;
            if (blocknum - 3 > 0) {
                ind = blocknum - 3; // maximal blocking is performed for the 3rd last point in the blocking curve if
                                    // that exists
            } else {
                ind = blocknum - 1;
            }
        } else {
            blockIndS = 1;
        }

        ind = (int) Math.max(ind, correlatorq - 1); // block at least until maximum sample time

        return ind; // return the blockIndex
    }

    // get intensity data for correlate() and correct for bleaching if required;
    // note that you should call calcIntensityTrace() before to obtain intTrace1 and
    // 2
    public double[] getIntensity(ImagePlus image, int px, int py, int mode, int initialframe, int finalframe) {
        // image: imp form which intensity will be taken
        // px, py: coordinates of pixel within image
        // mode: determines whether intensity for pixel 1 or pixel 2 is read, and in
        // case of DC-FCCS,
        // whether background1 or background 2 is to be subtracted from the intensity
        // trace
        // initialframe and finalframe provide the range of frames to be used
        int num = (finalframe - initialframe + 1);
        double[] intdat = new double[num + 1];
        double[] res = new double[5];
        float bckg;

        bckg = background;

        for (int x = 1; x <= num; x++) { // read data from all relevant pixels, depending on the selected frames and
                                         // binning
            for (int i = 0; i < binningX; i++) {
                for (int k = 0; k < binningY; k++) {
                    if (bgrloaded) {
                        bckg = (int) Math.round(bgrmean[px + i][py + k]);
                    }
                    intdat[x] += image.getStack().getProcessor(initialframe + x - 1).getf(px + i, py + k) - bckg;
                }
            }
        }

        return (intdat);
    }

    // data fitting
    public void prepareFit() {

        // always start with the same initial values; the initial values are first
        // defined in createImFCSFit();
        // subsequently they are read in setParameters() from the fit window. One thus
        // should set the prameters to
        // realistic values, otherwise the fit might not converge and hang the plugin
        // for a long time
        paraminitval[0] = initparam[0];
        paramfit[0] = true;

        paraminitval[1] = initparam[1];
        paramfit[1] = true;

        paraminitval[2] = initparam[2];
        paramfit[2] = false;

        paraminitval[3] = initparam[3];
        paramfit[3] = false;

        paraminitval[4] = initparam[4];
        paramfit[4] = true;

        paraminitval[5] = initparam[5];
        paramfit[5] = false;

        paraminitval[6] = initparam[6];
        paramfit[6] = false;

        paraminitval[7] = initparam[7];
        paramfit[7] = false;

        paraminitval[8] = initparam[8];
        paramfit[8] = false;

        paraminitval[9] = initparam[9];
        paramfit[9] = false;

        paraminitval[10] = initparam[10];
        paramfit[10] = false;

        // set q2 and q2 in setUI in the beginning

        pixeldimx = (pixelsize * 1000 / objmag * binningX) / Math.pow(10, 9);
        pixeldimy = (pixelsize * 1000 / objmag * binningY) / Math.pow(10, 9);
        if (sigmaZ <= 0.0 || sigmaZ > 100) {
            fitobsvol = obsvolFCS_ST2D1p(2); // TODO: currently is the default. sigmaZ is fixed
        } else {
            fitobsvol = obsvolFCS_ST2D1p(3);
        }
    }

    // calculation of the observation area; this is used in the Diffusion Law Plot
    // as the y-axis
    // the calculation of the observation area/volume is provided on our website in
    // CDF files (http://www.dbs.nus.edu.sg/lab/BFL/index.html)
    public double obsvolFCS_ST2D1p(int dim) {
        // general parameters
        double pi = 3.14159265359;
        double sqrpi = Math.sqrt(pi);
        double ax = pixeldimx;
        double ay = pixeldimy;
        double s = psfsize;
        double sz = lsthickness;
        double psfz = 2 * emlambda / Math.pow(10, 9.0) * 1.33 / Math.pow(NA, 2.0); // size of PSF in axial direction
        double szeff = Math.sqrt(1 / (Math.pow(sz, -2.0) + Math.pow(psfz, -2.0))); // convolution of two Gaussians
                                                                                   // depending on illumination profile
                                                                                   // and detection PSF
        double rx = ax * cfXshift / binningX;
        double ry = ay * cfYshift / binningY;

        // help variables, for t = 0, to write the full fit function
        double p00 = s;
        double p1x0 = ax;
        double p2x0 = ax;
        double p1y0 = ay;
        double p2y0 = ay;
        double pexpx0 = 2 * Math.exp(-Math.pow(p1x0 / p00, 2)) - 2;
        double perfx0 = 2 * p1x0 * Erf.erf(p1x0 / p00);
        double pexpy0 = 2 * Math.exp(-Math.pow(p1y0 / p00, 2)) - 2;
        double perfy0 = 2 * p1y0 * Erf.erf(p1y0 / p00);

        // return (p00/sqrpi * pexpx0 + perfx0) * (p00/sqrpi * pexpy0 + perfy0) *
        // Math.pow(sz, 2);
        if (dim == 2) {
            return 4 * Math.pow(ax * ay, 2) / ((p00 / sqrpi * pexpx0 + perfx0) * (p00 / sqrpi * pexpy0 + perfy0));
        } else {
            // return sqrpi * szeff * 4 * Math.pow(ax*ay, 2)/( (p00/sqrpi * pexpx0 + perfx0)
            // * (p00/sqrpi * pexpy0 + perfy0) );
            return 4 * Math.pow(ax * ay, 2) / ((p00 / sqrpi * pexpx0 + perfx0) * (p00 / sqrpi * pexpy0 + perfy0));
        }

    }

    public boolean fit(int cpx, int cpy, int kcf, String model) {
        pixfitted[kcf][cpx][cpy] = false;
        if (isGLS) { // select whether GLS or OLS fits are used and whether Bayes is used or not
            if (isBayes) { // Bayes fitting will call fitting routines multiple times for model comparisons
                double testq2 = q2;
                double testq3 = q3;
                if (testq2 * testq3 <= 0) {
                    return false;
                }
                if (initparam[1] < 0 || initparam[6] < 0 || initparam[8] < 0 || initparam[5] < 0 || initparam[7] < 0) {
                    return false;
                }
                bayesFit(cpx, cpy, kcf, model);
            } else {
                if (initparam[1] < 0 || initparam[6] < 0 || initparam[8] < 0 || initparam[5] < 0 || initparam[7] < 0) {
                    // JOptionPane.showMessageDialog(null, "Fit parameters are out of range. Either
                    // D < 0 or F < 0.");
                    return false;
                }
                GLSFit glsfit = new GLSFit();
                glsfit.doFit(cpx, cpy, kcf, model);

            }
        } else {
            if (isBayes) {
                double testq2 = q2;
                double testq3 = q3;
                if (testq2 * testq3 <= 0) {
                    return false;
                }
                if (initparam[1] < 0 || initparam[6] < 0 || initparam[8] < 0 || initparam[5] < 0 || initparam[7] < 0) {
                    return false;
                }
                bayesFit(cpx, cpy, kcf, model);
            } else { // TODO: standard fitting is the default: standard fit refers to no GLS and no
                     // bayes
                if (initparam[1] < 0 || initparam[6] < 0 || initparam[8] < 0 || initparam[5] < 0 || initparam[7] < 0) {
                    // JOptionPane.showMessageDialog(null, "Fit parameters are out of range. Either
                    // D < 0 or F < 0.");
                    return false;
                }
                standardFit sfit = new standardFit();
                sfit.doFit(cpx, cpy, kcf, model);
            }
        }
        return true;
    }

    public void initializeFitres(int a, int b, int c, int d) { // initializes arrays for fit parameters, chi2 and
                                                               // cross-correlation amount with NaN
        fitres = new double[a][b][c][d];
        chi2 = new double[a][b][c];
        // CCFq = new double[b][c];
        pixfitted = new boolean[a][b][c];
        for (int q = 0; q < a; q++) {
            for (int r = 0; r < b; r++) {
                for (int s = 0; s < c; s++) {
                    chi2[q][r][s] = Double.NaN;
                    pixfitted[q][r][s] = false;
                    for (int t = 0; t < d; t++) {
                        fitres[q][r][s][t] = Double.NaN;
                    }
                }
            }
        }
        for (int r = 0; r < b; r++) {
            for (int s = 0; s < c; s++) {
                // CCFq[r][s] = Double.NaN;
            }
        }
    }

    // fit procedure for Bayes Hypothesis testing
    public void bayesFit(int cpx, int cpy, int kcf, String model) { // at the moment we will allow onle one and
                                                                    // two-component fits
        // cpx, cpy: pixel coordinates
        // kcf (0, 1, 2): ACF1, ACF2, CCF
        // model: which fit model to be used; the models are defined below
        ParametricUnivariateFunction fitfunction;
        ParametricUnivariateFunction plotfunction;
        double[] paramem = new double[noparam];
        boolean[] paramfitmem = new boolean[noparam];
        double[][] covpar;
        double[] residuals;
        double[] sigma;
        double[] modprob = new double[2];
        double normprob;
        double residsum;
        double logresid;
        double prodsigma;
        double det;
        double pi = 3.14159265359;
        int boxsize = 200;
        int numel;

        // take the existing fitparameters and do a one or two component fit; if the D2
        // etc. are not specified, specify them and take D2 as 10 time D1 etc.
        // remember the fit parameter settings so they can be reset after the Bayes

        if (cpx == 0 && cpy == 0) {// TODO: set initparam as firs prior
            paramem[0] = 1.0; // N
            paramem[1] = 1.0 / Math.pow(10, 12); // D
            paramem[2] = 0.0 / Math.pow(10, 6); // Vx
            paramem[3] = 0.0 / Math.pow(10, 6); // Vy
            paramem[4] = 0.0; // G
            paramem[5] = 0.0; // F2
            paramem[6] = 0.0 / Math.pow(10, 12); // D2

        } else {// use subsequnt fit result from 2P not 1P (need to clarify) why?
            paramem[0] = storeTempFit.getTempFitres(0); // N
            paramem[1] = storeTempFit.getTempFitres(1); // D
            paramem[2] = storeTempFit.getTempFitres(2); // Vx
            paramem[3] = storeTempFit.getTempFitres(3); // Vy
            paramem[4] = storeTempFit.getTempFitres(4); // G
            paramem[5] = storeTempFit.getTempFitres(5); // F2
            paramem[6] = storeTempFit.getTempFitres(6); // D2

        }

        System.arraycopy(paramfit, 0, paramfitmem, 0, noparam);

        // set initial values for one-component fit
        paramfit[0] = true;
        paramfit[1] = true;
        paramfit[2] = false;
        paramfit[3] = false;
        paramfit[4] = true;
        paramfit[5] = false;
        paramfit[6] = false;
        paramfit[7] = false;
        paramfit[8] = false;
        initparam[0] = paramem[0];
        initparam[1] = paramem[1];
        initparam[2] = 0;
        initparam[3] = 0;
        initparam[4] = paramem[4];
        initparam[5] = 0;
        initparam[6] = 0;
        initparam[7] = 0;
        initparam[8] = 0;
        initparam[9] = 0;
        initparam[10] = 0;

        residsum = 0;
        if (!isGLS) {
            standardFit sfit = new standardFit();
            Map map;
            map = sfit.doFit(cpx, cpy, kcf, model);

            covpar = (double[][]) map.get("covariance");
            residuals = (double[]) map.get("residuals");
            sigma = (double[]) map.get("sigma");

            for (int i = 1; i < residuals.length; i++) {
                residsum += Math.pow(residuals[i], 2);
            }

            logresid = -0.5 * residsum;
        } else {
            GLSFit glsfit = new GLSFit();
            Map map;
            map = glsfit.doFit(cpx, cpy, kcf, model);

            covpar = (double[][]) map.get("covariance");
            residuals = (double[]) map.get("residuals");
            sigma = (double[]) map.get("sigma");

            RealMatrix mat = MatrixUtils.createRealMatrix(currentCovmats);
            RealMatrix matT = mat.transpose();
            DecompositionSolver solver = new LUDecomposition(matT).getSolver();
            RealVector resvec = new ArrayRealVector(residuals);
            RealVector solution = solver.solve(resvec);
            logresid = -0.5 * solution.dotProduct(resvec);
        }

        numel = 3;
        det = determinant(covpar, numel);

        prodsigma = 1;
        for (int i = 0; i < sigma.length; i++) {
            prodsigma *= sigma[i];
        }

        modprob[0] = Math.exp(0.5 * numel * Math.log(2 * pi) + 0.5 * Math.log(det) + logresid
                - Math.log(prodsigma * Math.pow(2 * boxsize, numel)));

        // set initial values for two-component fit
        paramfit[0] = true;
        paramfit[1] = true;
        paramfit[2] = false;
        paramfit[3] = false;
        paramfit[4] = true;
        paramfit[5] = true;
        paramfit[6] = true;
        paramfit[7] = false;
        paramfit[8] = false;
        initparam[0] = paramem[0]; // N
        initparam[1] = paramem[1]; // D
        initparam[2] = 0;
        initparam[3] = 0;
        initparam[4] = paramem[4]; // G
        initparam[5] = 0.5; // F2 50:50 prior
        initparam[6] = 0.1 * paramem[1]; // D2 10x slower prior
        initparam[7] = 0;
        initparam[8] = 0;
        initparam[9] = 0;
        initparam[10] = 0;

        residsum = 0;
        if (!isGLS) {
            standardFit sfit = new standardFit();
            Map map;
            map = sfit.doFit(cpx, cpy, kcf, model);

            covpar = (double[][]) map.get("covariance");
            residuals = (double[]) map.get("residuals");
            sigma = (double[]) map.get("sigma");

            for (int i = 1; i < residuals.length; i++) {
                residsum += Math.pow(residuals[i], 2);
            }

            logresid = -0.5 * residsum;
        } else {
            GLSFit glsfit = new GLSFit();
            Map map;
            map = glsfit.doFit(cpx, cpy, kcf, model);

            covpar = (double[][]) map.get("covariance");
            residuals = (double[]) map.get("residuals");
            sigma = (double[]) map.get("sigma");

            RealMatrix mat = MatrixUtils.createRealMatrix(currentCovmats);
            RealMatrix matT = mat.transpose();
            DecompositionSolver solver = new LUDecomposition(matT).getSolver();
            RealVector resvec = new ArrayRealVector(residuals);
            RealVector solution = solver.solve(resvec);
            logresid = -0.5 * solution.dotProduct(resvec);
        }

        numel = 5;
        det = determinant(covpar, numel);

        prodsigma = 1;
        for (int i = 0; i < sigma.length; i++) {
            prodsigma *= sigma[i];
        }

        modprob[1] = Math.exp(0.5 * numel * Math.log(2 * pi) + 0.5 * Math.log(det) + logresid
                - Math.log(prodsigma * Math.pow(2 * boxsize, numel)));

        normprob = 0;
        for (int i = 0; i < modprob.length; i++) // calculate the normalization for the modle probabilities
        {
            normprob += modprob[i];
        }

        // store model probability for bayesian analysis
        modelprobP1 = modprob[0] / normprob;
        modelprobP2 = modprob[1] / normprob;

    }

    // standard Nonlinear Least Squares fit
    public class standardFit extends AbstractCurveFitter {

        private int numfreefitpar;

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess;

            int i = 0; // add data
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                weights[i] = point.getWeight();
                i++;
            }

            int numparfit = 0; // count how many paramters will be fit
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    numparfit++;
                }
            }

            numfreefitpar = numparfit;

            initialGuess = new double[numfreefitpar]; // setup the initial guesses for the parameters
            int num = 0;
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    initialGuess[num] = paraminitval[i];
                    num++;
                }
            }

            ParametricUnivariateFunction fitfunction;
            if (currentmodel.equals("FCS")) { // select the fit model to be used; extra fit models can be added here
                fitfunction = new FCS_3p(); // TODO: currently FCS_3p is set to default
            } else {
                fitfunction = new FCS_3p();
            }
            // else {
            // //fitfunction = new FCCS_2p();
            // }

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(
                    fitfunction, points);

            return new LeastSquaresBuilder().maxEvaluations(fitMaxIterations).maxIterations(fitMaxEvaluations)
                    .start(initialGuess).target(target).weight(new DiagonalMatrix(weights))
                    .model(model.getModelFunction(), model.getModelFunctionJacobian()).build();
        }

        public Map doFit(int cpx, int cpy, int kcf, String model) {
            // standardFit fitter = new standardFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();

            currentmodel = model; // set the presently used model

            ParametricUnivariateFunction fitfunction;
            if (currentmodel.equals("FCS")) { // select the fit model to be used; extra fit models can be added here
                fitfunction = new FCS_3p();
            } else {
                fitfunction = new FCS_3p();
            }

            // else {
            // fitfunction = new FCCS_2p();
            // }
            fitstart = 1;
            fitend = chanum - 1;// TODO: fix

            if (fitstart < 1) {
                IJ.showMessage("Illegal Value for Fit start");
                fitstart = 1;
            }
            if (fitend >= chanum) {
                IJ.showMessage("Illegal Value for Fit end");
                fitend = chanum - 1;
            }

            // Add points here
            for (int i = fitstart; i <= fitend; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1 / varacf[kcf][cpx][cpy][i], lagtime[i],
                        acf[kcf][cpx][cpy][i]);
                points.add(point);
            }

            Map<String, Object> map = new HashMap<>();
            try {
                final Optimum topt = getOptimizer().optimize(getProblem(points));
                double result[] = topt.getPoint().toArray();
                double tmpres[] = topt.getResiduals().toArray();
                double[] tres = new double[chanum - 1];

                // store results and create fit and residual function for the plot
                for (int i = 1; i < chanum; i++) {
                    if (i >= fitstart && i <= fitend) {
                        fitacf[kcf][cpx][cpy][i] = fitfunction.value(lagtime[i], result); // calculate the fit function
                        if (i == 0) {
                            res[kcf][cpx][cpy][i] = acf[kcf][cpx][cpy][i] - fitacf[kcf][cpx][cpy][i]; // calculate the
                                                                                                      // residuals
                        } else {
                            res[kcf][cpx][cpy][i] = tmpres[i - fitstart]; // use the weighted residuals
                            tres[i - 1] = acf[kcf][cpx][cpy][i] - fitacf[kcf][cpx][cpy][i]; // unweighted residuals for
                                                                                            // model probability
                                                                                            // calculations
                        }
                    } else {
                        fitacf[kcf][cpx][cpy][i] = 0;
                        res[kcf][cpx][cpy][i] = 0;
                        if (i > 0) {
                            tres[i - 1] = 0;
                        }
                    }
                }

                chi2[kcf][cpx][cpy] = 0; // initialize chi2 value for this pixel
                for (int i = fitstart; i <= fitend; i++) {
                    chi2[kcf][cpx][cpy] += Math.pow(res[kcf][cpx][cpy][i], 2)
                            / ((fitend - fitstart) - numfreefitpar - 1); // calculate chi2 value; do not include the 0
                                                                         // lagtime channel which contains shot noise
                }

                int num = 0; // store the fit results
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        fitres[kcf][cpx][cpy][i] = result[num]; // for free parameters use the fit results
                        storeTempFit.storeTempFitres(i, result[num]);
                        num++;
                    } else {
                        fitres[kcf][cpx][cpy][i] = paraminitval[i]; // for fixed parameters use the initial values
                        storeTempFit.storeTempFitres(i, paraminitval[i]);
                    }

                }
                pixfitted[kcf][cpx][cpy] = true;

                // update parameters in fit window (not applicable)
                map.put("results", result);
                map.put("covariance", topt.getCovariances(1).getData());
                map.put("residuals", tres);
                map.put("sigma", topt.getSigma(1).toArray());

            } catch (Exception e) {
                IJ.log(e.getClass().getName() + " at pixel " + Integer.toString(cpy) + " - " + Integer.toString(cpx));
                ArrayList<Double> result = new ArrayList<>();
                double[] tres = new double[chanum - 1];
                Arrays.fill(tres, 1);
                double[] tsig = new double[chanum - 1];
                Arrays.fill(tsig, 1);
                int num = 0;
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        fitres[kcf][cpx][cpy][i] = Double.NaN;
                        result.add(initparam[i]); // return the initial values
                        num++;
                    }
                }
                for (int i = 1; i < chanum; i++) {
                    fitacf[kcf][cpx][cpy][i] = 0.0;
                    res[kcf][cpx][cpy][i] = 0.0; // calculate the residuals
                }
                double[][] tcov = new double[num][num];
                for (double[] tcov1 : tcov) {
                    Arrays.fill(tcov1, 1);
                }
                map.put("results", result.toArray());
                map.put("covariance", tcov);
                map.put("residuals", tres);
                map.put("sigma", tsig);
                pixfitted[kcf][cpx][cpy] = false;
            }

            return map;
        }
    }

    // Generalized Least Squares fit; requires covariance matrix as calculated in
    // calculateCF()
    public class GLSFit extends AbstractCurveFitter {

        public int numfreefitpar;

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess;

            int i = 0; // add data
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                weights[i] = point.getWeight();
                i += 1;
            }

            int numparfit = 0; // count how many paramters will be fit
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    numparfit++;
                }
            }

            numfreefitpar = numparfit;

            initialGuess = new double[numparfit]; // setup the initial guesses for the parameters
            int num = 0;
            for (i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    initialGuess[num] = paraminitval[i];
                    num++;
                }
            }

            ParametricUnivariateFunction glsfitfunction;
            glsfitfunction = new GLS_fitFunction();

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(
                    glsfitfunction, points);

            return new LeastSquaresBuilder().maxEvaluations(fitMaxIterations).maxIterations(fitMaxEvaluations)
                    .start(initialGuess).target(target).weight(new DiagonalMatrix(weights))
                    .model(model.getModelFunction(), model.getModelFunctionJacobian()).build();
        }

        public Map doFit(int cpx, int cpy, int kcf, String model) {
            // GLSFit fitter = new GLSFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();

            dataTransform(acf[kcf][cpx][cpy], currentCovmats); // transform data
            transTheoreticalACF = new double[chanum - 1];
            transTheoreticalGradientACF = new double[noparam][chanum - 1];

            // Add points here
            for (int i = 1; i < chanum; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1, lagtime[i], transACF.getEntry(i - 1));
                points.add(point);
            }
            currentmodel = model;
            ParametricUnivariateFunction plotfunction;
            if (currentmodel.equals("FCS")) { // select the fit model to be used; extra fit models can be added here
                plotfunction = new FCS_3p();
            } else {
                plotfunction = new FCS_3p();
            }

            // else {
            // plotfunction = new FCCS_2p();
            // }
            Map<String, Object> map = new HashMap<>();
            try {
                final Optimum topt = getOptimizer().optimize(getProblem(points));
                double result[] = topt.getPoint().toArray();
                double tmpres[] = topt.getResiduals().toArray();
                double[] tres = new double[chanum - 1];

                // store results and create fit and residual function for the plot
                for (int i = 0; i < chanum; i++) {
                    fitacf[kcf][cpx][cpy][i] = plotfunction.value(lagtime[i], result); // calculate the fit function
                    if (i == 0) {
                        res[kcf][cpx][cpy][i] = acf[kcf][cpx][cpy][i] - fitacf[kcf][cpx][cpy][i]; // calculate the
                                                                                                  // residuals
                    } else {
                        res[kcf][cpx][cpy][i] = tmpres[i - 1];
                        tres[i - 1] = acf[kcf][cpx][cpy][i] - fitacf[kcf][cpx][cpy][i]; // unweighted residuals for
                                                                                        // model probability
                                                                                        // calculations
                    }
                }

                chi2[kcf][cpx][cpy] = 0; // initialize chi2 value for this pixel
                for (int i = 1; i < chanum; i++) {
                    chi2[kcf][cpx][cpy] += Math.pow(res[kcf][cpx][cpy][i], 2) / (chanum - numfreefitpar - 1); // calculate
                    // shot
                    // nosie
                }

                int num = 0; // store the fit results
                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        fitres[kcf][cpx][cpy][i] = result[num]; // for free parameters use the fit results
                        storeTempFit.storeTempFitres(i, result[num]);
                        num++;
                    } else {
                        fitres[kcf][cpx][cpy][i] = paraminitval[i]; // for fixed parameters use the initial values
                        storeTempFit.storeTempFitres(i, paraminitval[i]);
                    }
                }
                pixfitted[kcf][cpx][cpy] = true;

                // update parameters in fit window Not applicable
                map.put("results", result);
                map.put("covariance", topt.getCovariances(1).getData());
                map.put("residuals", tres);
                map.put("sigma", topt.getSigma(1).toArray());

            } catch (Exception e) {
                IJ.log(e.getClass().getName() + " at pixel " + Integer.toString(cpy) + " - " + Integer.toString(cpx));
                ArrayList<Double> result = new ArrayList<>();
                double[] tres = new double[chanum - 1];
                Arrays.fill(tres, 1);
                double[] tsig = new double[chanum - 1];
                Arrays.fill(tsig, 1);
                int num = 0;
                /*
                 * for (int i = 0; i < noparam; i++) { if ( paramfit[i] == true ) {
                 * result.add(initparam[i]); // return the initial values num++; } }
                 */

                for (int i = 0; i < noparam; i++) {
                    if (paramfit[i] == true) {
                        fitres[kcf][cpx][cpy][i] = Double.NaN;
                        result.add(initparam[i]); // return the initial values
                        num++;
                    }
                }
                for (int i = 0; i < chanum; i++) {
                    fitacf[kcf][cpx][cpy][i] = 0.0;
                    res[kcf][cpx][cpy][i] = 0.0; // calculate the residuals

                }
                double[][] tcov = new double[num][num];
                for (double[] tcov1 : tcov) {
                    Arrays.fill(tcov1, 1);
                }
                map.put("results", result.toArray());
                map.put("covariance", tcov);
                map.put("residuals", tres);
                map.put("sigma", tsig);
                pixfitted[kcf][cpx][cpy] = false;
            }

            return map;
        }
    }

    // Transform data using the regularized covariance matrix for a generalized
    // least squares fit
    public void dataTransform(double[] corav, double[][] covmats) { // this needs to be run before the Bayes fit or any
                                                                    // generalized least suqares
        // corav: correlation function
        // covmats: covariance matrix
        double[] cortmp = new double[chanum - 1];
        for (int x = 0; x < chanum - 1; x++) { // remove zero lagtime channel from the correlation as it is not fit
                                               // along
            cortmp[x] = corav[x + 1];
        }

        RealMatrix mat = MatrixUtils.createRealMatrix(covmats); // lower triangular matrix of the CholeskyDecomposition
                                                                // of covmats
        RealMatrix matL;

        try { // perfrom the Cholesky decomposition and determine the lower triangular matrix
            CholeskyDecomposition CD = new CholeskyDecomposition(mat);
            matL = CD.getL();
        } catch (NonSquareMatrixException | NonSymmetricMatrixException | NonPositiveDefiniteMatrixException ex) {
            IJ.log(ex.getMessage());
            throw ex;
        }

        DecompositionSolver solver = new LUDecomposition(matL).getSolver(); // solve for a new correlation vector with
                                                                            // independent elements
        RealVector constants = new ArrayRealVector(cortmp);
        RealVector solution = solver.solve(constants);
        lowerCholDecCovmats = matL;
        transACF = solution;
    }

    /*
     * Fit functions
     * 
     * class FCS_3p implements ParametricUnivariateFunction: FCS fit including
     * diffusion, flow and spatial cross-correlation for ITIR-FCS and SPIM-FCS;
     * accepts up to 3 components/particles with areas of rectangular shape class
     * FCCS_2p implements ParametricUnivariateFunction: DC-FCCS fit; assumes only
     * diffusion; accepts up to 2 components/particles class GLS_fitFunction
     * implements ParametricUnivariateFunction: generalized least squares fit;
     * transforms the model functions using the regularized covarince matrix public
     * double obsvolFCS_ST2D1p(int dim): Calculate the observation area/volume
     * public void dataTransform(double[] corav, double[][] covmats): transforms the
     * ACF using the regularized covarince matrix for the GLS fit
     * 
     */
    // FCS model: 3D fit assuming 3 components and flow in x and y direction; this
    // is the general fit formula; parameters can be set to 0 to obtain simpler
    // models
    // the models and their derivation are provided on our website in CDF files
    // (http://staff.science.nus.edu.sg/~chmwt/)
    class FCS_3p implements ParametricUnivariateFunction {
        // general parameters

        double pi = 3.14159265359;
        double sqrpi = Math.sqrt(pi);
        double ax = pixeldimx;
        double ay = pixeldimy;
        double s = psfsize;
        double sz = lsthickness;
        double psfz = 2 * emlambda / Math.pow(10, 9.0) * 1.33 / Math.pow(NA, 2.0); // size of PSF in axial direction
        // double szeff = Math.sqrt( 1 / ( Math.pow(sz, -2.0) + Math.pow(psfz, -2.0) )
        // ); // convolution of two Gaussians depending on illumination profile and
        // detection PSF
        double szeff = sz;
        double rx = ax * cfXshift / binningX;
        double ry = ay * cfYshift / binningY;

        @Override
        public double[] gradient(double x, double[] params) {
            double[] pareq = new double[noparam];
            int num = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i]) {
                    pareq[i] = params[num];
                    num++;
                } else {
                    pareq[i] = paraminitval[i];
                }
            }

            // note that x is used here as the time variable instead of t; that can be come
            // confusing as x and y are used in the names for the paramaters to indicate
            // spatial directions
            // pareq[0] = N
            // pareq[1] = D
            // pareq[2] = vx
            // pareq[3] = vy
            // pareq[4] = G
            // pareq[5] = F2
            // pareq[6] = D2
            // pareq[7] = F3
            // pareq[8] = D3
            // pareq[9] = Ftrip
            // pareq[10] = Ttrip
            // COMPONENT1
            // help variables, which are dependent on time, to write the full function
            double p0t = Math.sqrt(4 * pareq[1] * x + Math.pow(s, 2));
            double p1xt = ax + rx - pareq[2] * x;
            double p2xt = ax - rx + pareq[2] * x;
            double p3xt = rx - pareq[2] * x;
            double p4xt = 2 * Math.pow(ax, 2) + 3 * Math.pow(rx, 2) - 6 * x * rx * pareq[2]
                    + 3 * Math.pow(x * pareq[2], 2);
            double p5xt = Math.pow(p3xt, 2) + Math.pow(p1xt, 2);
            double p6xt = Math.pow(p3xt, 2) + Math.pow(p2xt, 2);
            double p7xt = 2 * (Math.pow(ax, 2) + Math.pow(rx, 2) - 2 * x * rx * pareq[2] + Math.pow(x * pareq[2], 2));
            double p1yt = ay + ry - pareq[3] * x;
            double p2yt = ay - ry + pareq[3] * x;
            double p3yt = ry - pareq[3] * x;
            double p4yt = 2 * Math.pow(ay, 2) + 3 * Math.pow(ry, 2) - 6 * x * ry * pareq[3]
                    + 3 * Math.pow(x * pareq[3], 2);
            double p5yt = Math.pow(p3yt, 2) + Math.pow(p1yt, 2);
            double p6yt = Math.pow(p3yt, 2) + Math.pow(p2yt, 2);
            double p7yt = 2 * (Math.pow(ay, 2) + Math.pow(ry, 2) - 2 * x * ry * pareq[3] + Math.pow(x * pareq[3], 2));
            double pexpxt = Math.exp(-Math.pow(p1xt / p0t, 2)) + Math.exp(-Math.pow(p2xt / p0t, 2))
                    - 2 * Math.exp(-Math.pow(p3xt / p0t, 2));
            double perfxt = p1xt * Erf.erf(p1xt / p0t) + p2xt * Erf.erf(p2xt / p0t) - 2 * p3xt * Erf.erf(p3xt / p0t);
            double dDpexpxt = 2 * Math.exp(-p4xt / Math.pow(p0t, 2)) * (Math.exp(p5xt / Math.pow(p0t, 2))
                    + Math.exp(p6xt / Math.pow(p0t, 2)) - 2 * Math.exp(p7xt / Math.pow(p0t, 2)));
            double dvxperfxt = (Erf.erf(p2xt / p0t) + 2 * Erf.erf(p3xt / p0t) - Erf.erf(p1xt / p0t)) * x;
            double pexpyt = Math.exp(-Math.pow(p1yt / p0t, 2)) + Math.exp(-Math.pow(p2yt / p0t, 2))
                    - 2 * Math.exp(-Math.pow(p3yt / p0t, 2));
            double dDpexpyt = 2 * Math.exp(-p4yt / Math.pow(p0t, 2)) * (Math.exp(p5yt / Math.pow(p0t, 2))
                    + Math.exp(p6yt / Math.pow(p0t, 2)) - 2 * Math.exp(p7yt / Math.pow(p0t, 2)));
            double dvyperfyt = (Erf.erf(p2yt / p0t) + 2 * Erf.erf(p3yt / p0t) - Erf.erf(p1yt / p0t)) * x;
            double perfyt = p1yt * Erf.erf(p1yt / p0t) + p2yt * Erf.erf(p2yt / p0t) - 2 * p3yt * Erf.erf(p3yt / p0t);

            // CF for the lateral dimension (x, y) and its derivative for D
            double plat = (p0t / sqrpi * pexpxt + perfxt) * (p0t / sqrpi * pexpyt + perfyt)
                    / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double dDplat = (1 / (sqrpi * p0t))
                    * (dDpexpyt * x * (p0t / sqrpi * pexpxt + perfxt) + dDpexpxt * x * (p0t / sqrpi * pexpyt + perfyt))
                    / (4 * Math.pow(ax * ay, 2) / fitobsvol);

            // CF for the axial dimension (z) and its derivative for D
            double pspim = 1 / Math.sqrt(1 + (4 * pareq[1] * x) / Math.pow(szeff, 2));
            double dDpspim = -4 * x
                    / (2 * Math.pow(szeff, 2) * Math.pow(Math.sqrt(1 + (4 * pareq[1] * x) / Math.pow(szeff, 2)), 3));

            double acf1 = plat * pspim;

            // COMPONENT2
            // help variables, which are dependent on time, to write the full function
            double p0t2 = Math.sqrt(4 * pareq[6] * x + Math.pow(s, 2));
            double p1xt2 = ax + rx - pareq[2] * x;
            double p2xt2 = ax - rx + pareq[2] * x;
            double p3xt2 = rx - pareq[2] * x;
            double p4xt2 = 2 * Math.pow(ax, 2) + 3 * Math.pow(rx, 2) - 6 * x * rx * pareq[2]
                    + 3 * Math.pow(x * pareq[2], 2);
            double p5xt2 = Math.pow(p3xt2, 2) + Math.pow(p1xt2, 2);
            double p6xt2 = Math.pow(p3xt2, 2) + Math.pow(p2xt2, 2);
            double p7xt2 = 2 * (Math.pow(ax, 2) + Math.pow(rx, 2) - 2 * x * rx * pareq[2] + Math.pow(x * pareq[2], 2));
            double p1yt2 = ay + ry - pareq[3] * x;
            double p2yt2 = ay - ry + pareq[3] * x;
            double p3yt2 = ry - pareq[3] * x;
            double p4yt2 = 2 * Math.pow(ay, 2) + 3 * Math.pow(ry, 2) - 6 * x * ry * pareq[3]
                    + 3 * Math.pow(x * pareq[3], 2);
            double p5yt2 = Math.pow(p3yt2, 2) + Math.pow(p1yt2, 2);
            double p6yt2 = Math.pow(p3yt2, 2) + Math.pow(p2yt2, 2);
            double p7yt2 = 2 * (Math.pow(ay, 2) + Math.pow(ry, 2) - 2 * x * ry * pareq[3] + Math.pow(x * pareq[3], 2));
            double pexpxt2 = Math.exp(-Math.pow(p1xt2 / p0t2, 2)) + Math.exp(-Math.pow(p2xt2 / p0t2, 2))
                    - 2 * Math.exp(-Math.pow(p3xt2 / p0t2, 2));
            double perfxt2 = p1xt2 * Erf.erf(p1xt2 / p0t2) + p2xt2 * Erf.erf(p2xt2 / p0t2)
                    - 2 * p3xt2 * Erf.erf(p3xt2 / p0t2);
            double dDpexpxt2 = 2 * Math.exp(-p4xt2 / Math.pow(p0t2, 2)) * (Math.exp(p5xt2 / Math.pow(p0t2, 2))
                    + Math.exp(p6xt2 / Math.pow(p0t2, 2)) - 2 * Math.exp(p7xt2 / Math.pow(p0t2, 2)));
            double dvxperfxt2 = (Erf.erf(p2xt2 / p0t2) + 2 * Erf.erf(p3xt2 / p0t2) - Erf.erf(p1xt2 / p0t2)) * x;
            double pexpyt2 = Math.exp(-Math.pow(p1yt2 / p0t2, 2)) + Math.exp(-Math.pow(p2yt2 / p0t2, 2))
                    - 2 * Math.exp(-Math.pow(p3yt2 / p0t2, 2));
            double dDpexpyt2 = 2 * Math.exp(-p4yt2 / Math.pow(p0t2, 2)) * (Math.exp(p5yt2 / Math.pow(p0t2, 2))
                    + Math.exp(p6yt2 / Math.pow(p0t2, 2)) - 2 * Math.exp(p7yt2 / Math.pow(p0t2, 2)));
            double dvyperfyt2 = (Erf.erf(p2yt2 / p0t2) + 2 * Erf.erf(p3yt2 / p0t2) - Erf.erf(p1yt2 / p0t2)) * x;
            double perfyt2 = p1yt2 * Erf.erf(p1yt2 / p0t2) + p2yt2 * Erf.erf(p2yt2 / p0t2)
                    - 2 * p3yt2 * Erf.erf(p3yt2 / p0t2);

            // CF for the lateral dimension (x, y) and its derivative for D
            double plat2 = (p0t2 / sqrpi * pexpxt2 + perfxt2) * (p0t2 / sqrpi * pexpyt2 + perfyt2)
                    / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double dDplat2 = (1 / (sqrpi * p0t2))
                    * (dDpexpyt2 * x * (p0t2 / sqrpi * pexpxt2 + perfxt2)
                            + dDpexpxt2 * x * (p0t2 / sqrpi * pexpyt2 + perfyt2))
                    / (4 * Math.pow(ax * ay, 2) / fitobsvol);

            // CF for the axial dimension (z) and its derivative for D
            double pspim2 = 1 / Math.sqrt(1 + (4 * pareq[6] * x) / Math.pow(szeff, 2));
            double dDpspim2 = -4 * x
                    / (2 * Math.pow(szeff, 2) * Math.pow(Math.sqrt(1 + (4 * pareq[6] * x) / Math.pow(szeff, 2)), 3));

            double acf2 = plat2 * pspim2;

            // COMPONENT3
            // help variables, which are dependent on time, to write the full function
            double p0t3 = Math.sqrt(4 * pareq[8] * x + Math.pow(s, 2));
            double p1xt3 = ax + rx - pareq[2] * x;
            double p2xt3 = ax - rx + pareq[2] * x;
            double p3xt3 = rx - pareq[2] * x;
            double p4xt3 = 2 * Math.pow(ax, 2) + 3 * Math.pow(rx, 2) - 6 * x * rx * pareq[2]
                    + 3 * Math.pow(x * pareq[2], 2);
            double p5xt3 = Math.pow(p3xt2, 2) + Math.pow(p1xt2, 2);
            double p6xt3 = Math.pow(p3xt2, 2) + Math.pow(p2xt2, 2);
            double p7xt3 = 2 * (Math.pow(ax, 2) + Math.pow(rx, 2) - 2 * x * rx * pareq[2] + Math.pow(x * pareq[2], 2));
            double p1yt3 = ay + ry - pareq[3] * x;
            double p2yt3 = ay - ry + pareq[3] * x;
            double p3yt3 = ry - pareq[3] * x;
            double p4yt3 = 2 * Math.pow(ay, 2) + 3 * Math.pow(ry, 2) - 6 * x * ry * pareq[3]
                    + 3 * Math.pow(x * pareq[3], 2);
            double p5yt3 = Math.pow(p3yt3, 2) + Math.pow(p1yt3, 2);
            double p6yt3 = Math.pow(p3yt3, 2) + Math.pow(p2yt3, 2);
            double p7yt3 = 2 * (Math.pow(ay, 2) + Math.pow(ry, 2) - 2 * x * ry * pareq[3] + Math.pow(x * pareq[3], 2));
            double pexpxt3 = Math.exp(-Math.pow(p1xt3 / p0t3, 2)) + Math.exp(-Math.pow(p2xt3 / p0t3, 2))
                    - 2 * Math.exp(-Math.pow(p3xt3 / p0t3, 2));
            double perfxt3 = p1xt3 * Erf.erf(p1xt3 / p0t3) + p2xt3 * Erf.erf(p2xt3 / p0t3)
                    - 2 * p3xt3 * Erf.erf(p3xt3 / p0t3);
            double dDpexpxt3 = 2 * Math.exp(-p4xt3 / Math.pow(p0t3, 2)) * (Math.exp(p5xt3 / Math.pow(p0t3, 2))
                    + Math.exp(p6xt3 / Math.pow(p0t3, 2)) - 2 * Math.exp(p7xt3 / Math.pow(p0t3, 2)));
            double dvxperfxt3 = (Erf.erf(p2xt3 / p0t3) + 2 * Erf.erf(p3xt3 / p0t3) - Erf.erf(p1xt3 / p0t3)) * x;
            double pexpyt3 = Math.exp(-Math.pow(p1yt3 / p0t3, 2)) + Math.exp(-Math.pow(p2yt3 / p0t3, 2))
                    - 2 * Math.exp(-Math.pow(p3yt3 / p0t3, 2));
            double dDpexpyt3 = 2 * Math.exp(-p4yt3 / Math.pow(p0t3, 2)) * (Math.exp(p5yt3 / Math.pow(p0t3, 2))
                    + Math.exp(p6yt3 / Math.pow(p0t3, 2)) - 2 * Math.exp(p7yt3 / Math.pow(p0t3, 2)));
            double dvyperfyt3 = (Erf.erf(p2yt3 / p0t3) + 2 * Erf.erf(p3yt3 / p0t3) - Erf.erf(p1yt3 / p0t3)) * x;
            double perfyt3 = p1yt3 * Erf.erf(p1yt3 / p0t3) + p2yt3 * Erf.erf(p2yt3 / p0t3)
                    - 2 * p3yt3 * Erf.erf(p3yt3 / p0t3);

            // TRIPLET
            double triplet = 1 + pareq[9] / (1 - pareq[9]) * Math.exp(-x / pareq[10]);
            double dtripletFtrip = Math.exp(-x / pareq[10])
                    * (1 / (1 - pareq[9]) + pareq[9] / Math.pow(1 - pareq[9], 2));
            double dtripletTtrip = Math.exp(-x / pareq[10]) * (pareq[9] * x)
                    / ((1 - pareq[9]) * Math.pow(pareq[10], 2));

            // CF for the lateral dimension (x, y) and its derivative for D
            double plat3 = (p0t3 / sqrpi * pexpxt3 + perfxt3) * (p0t3 / sqrpi * pexpyt3 + perfyt3)
                    / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double dDplat3 = (1 / (sqrpi * p0t3))
                    * (dDpexpyt3 * x * (p0t3 / sqrpi * pexpxt3 + perfxt3)
                            + dDpexpxt3 * x * (p0t3 / sqrpi * pexpyt3 + perfyt3))
                    / (4 * Math.pow(ax * ay, 2) / fitobsvol);

            // CF for the axial dimension (z) and its derivative for D
            double pspim3 = 1 / Math.sqrt(1 + (4 * pareq[8] * x) / Math.pow(szeff, 2));
            double dDpspim3 = -4 * x
                    / (2 * Math.pow(szeff, 2) * Math.pow(Math.sqrt(1 + (4 * pareq[8] * x) / Math.pow(szeff, 2)), 3));

            double acf3 = plat3 * pspim3;

            double pf1 = (1 - pareq[5] - pareq[7]) / (1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7]);
            double pf2 = (Math.pow(q2, 2) * pareq[5]) / (1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7]);
            double pf3 = (Math.pow(q3, 2) * pareq[7]) / (1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7]);
            double dfnom = Math.pow(1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7], 3);
            double df21 = 1 - pareq[5] - pareq[7] + q2 * pareq[5] - q3 * pareq[7] + 2 * q2 * pareq[7] - 2 * q2;
            double df22 = Math.pow(q2, 2) * (1 + pareq[5] - pareq[7] - q2 * pareq[5] + q3 * pareq[7]);
            double df23 = 2 * pareq[7] * Math.pow(q3, 2) * (1 - q2);
            double df31 = 1 - pareq[5] - pareq[7] - q2 * pareq[5] + 2 * q3 * pareq[5] - 2 * q3 + q3 * pareq[7];
            double df32 = 2 * pareq[5] * Math.pow(q2, 2) * (1 - q3);
            double df33 = Math.pow(q3, 2) * (1 - pareq[5] + pareq[7] + q2 * pareq[5] - q3 * pareq[7]);

            double pacf = (1 / pareq[0])
                    * ((1 - pareq[5] - pareq[7]) * acf1 + Math.pow(q2, 2) * pareq[5] * acf2
                            + Math.pow(q3, 2) * pareq[7] * acf3)
                    / Math.pow(1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7], 2) * triplet + pareq[4];

            double[] grad = new double[] {
                    (-1 / Math.pow(pareq[0], 2)) * (pf1 * acf1 + pf2 * acf2 + pf3 * acf3) * triplet,
                    (1 / pareq[0]) * pf1 * (plat * dDpspim + pspim * dDplat),
                    (1 / pareq[0]) * (pf1 * ((p0t / sqrpi * pexpyt + perfyt) * dvxperfxt) * pspim
                            / (4 * Math.pow(ax * ay, 2) / fitobsvol)
                            + pf2 * ((p0t2 / sqrpi * pexpyt2 + perfyt2) * dvxperfxt2) * pspim2
                                    / (4 * Math.pow(ax * ay, 2) / fitobsvol)
                            + pf3 * ((p0t3 / sqrpi * pexpyt3 + perfyt3) * dvxperfxt3) * pspim3
                                    / (4 * Math.pow(ax * ay, 2) / fitobsvol))
                            * triplet,
                    (1 / pareq[0]) * (pf1 * ((p0t / sqrpi * pexpxt + perfxt) * dvyperfyt) * pspim
                            / (4 * Math.pow(ax * ay, 2) / fitobsvol)
                            + pf2 * ((p0t2 / sqrpi * pexpxt2 + perfxt2) * dvyperfyt2) * pspim2
                                    / (4 * Math.pow(ax * ay, 2) / fitobsvol)
                            + pf3 * ((p0t3 / sqrpi * pexpxt3 + perfxt3) * dvyperfyt3) * pspim3
                                    / (4 * Math.pow(ax * ay, 2) / fitobsvol))
                            * triplet,
                    1, (1 / pareq[0]) * (1 / dfnom) * (df21 * acf1 + df22 * acf2 + df23 * acf3) * triplet,
                    (1 / pareq[0]) * pf2 * (plat2 * dDpspim2 + pspim2 * dDplat2) * triplet,
                    (1 / pareq[0]) * (1 / dfnom) * (df31 * acf1 + df32 * acf2 + df33 * acf3) * triplet,
                    (1 / pareq[0]) * pf3 * (plat3 * dDpspim3 + pspim3 * dDplat3) * triplet, dtripletFtrip * pacf,
                    dtripletTtrip * pacf };

            double[] gradret = new double[num]; // return the gradients of the fit model in respect to the fit
                                                // parameters
            num = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i] == true) {
                    gradret[num] = grad[i];
                    num++;
                }
            }

            return gradret;
        }

        @Override
        public double value(double x, double[] params) {
            double[] pareq = new double[noparam];
            int num = 0;
            for (int i = 0; i < noparam; i++) {
                if (paramfit[i]) {
                    pareq[i] = params[num];
                    num++;
                } else {
                    pareq[i] = paraminitval[i];
                }
            }

            // note that x is used here as the time variable instead of t; that can be come
            // confusing as x and y are used in the names for the paramaters to indicate
            // spatial directions
            // pareq[0] = N
            // pareq[1] = D
            // pareq[2] = vx
            // pareq[3] = vy
            // pareq[4] = G
            // pareq[5] = F2
            // pareq[6] = D2
            // pareq[7] = F3
            // pareq[8] = D3
            // pareq[9] = Ftrip
            // pareq[10] = Dtrip
            // q2 and q3, the brightness of the second and third components are fixed
            // parameters and have been globaly defined; see prepareFit()
            // COMPONENT 1
            // help variables, which are dependent on time, to write the full function
            double p0t = Math.sqrt(4 * pareq[1] * x + Math.pow(s, 2));
            double p1xt = ax + rx - pareq[2] * x;
            double p2xt = ax - rx + pareq[2] * x;
            double p3xt = rx - pareq[2] * x;
            double p1yt = ay + ry - pareq[3] * x;
            double p2yt = ay - ry + pareq[3] * x;
            double p3yt = ry - pareq[3] * x;
            double pexpxt = Math.exp(-Math.pow(p1xt / p0t, 2)) + Math.exp(-Math.pow(p2xt / p0t, 2))
                    - 2 * Math.exp(-Math.pow(p3xt / p0t, 2));
            double perfxt = p1xt * Erf.erf(p1xt / p0t) + p2xt * Erf.erf(p2xt / p0t) - 2 * p3xt * Erf.erf(p3xt / p0t);
            double pexpyt = Math.exp(-Math.pow(p1yt / p0t, 2)) + Math.exp(-Math.pow(p2yt / p0t, 2))
                    - 2 * Math.exp(-Math.pow(p3yt / p0t, 2));
            double perfyt = p1yt * Erf.erf(p1yt / p0t) + p2yt * Erf.erf(p2yt / p0t) - 2 * p3yt * Erf.erf(p3yt / p0t);

            double pplane1 = (p0t / sqrpi * pexpxt + perfxt) * (p0t / sqrpi * pexpyt + perfyt)
                    / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double pspim1 = 1 / Math.sqrt(1 + (4 * pareq[1] * x) / Math.pow(szeff, 2));
            double acf1 = pplane1 * pspim1;

            // COMPONENT 2
            // help variables, which are dependent on time, to write the full function
            double p0t2 = Math.sqrt(4 * pareq[6] * x + Math.pow(s, 2));
            double p1xt2 = ax + rx - pareq[2] * x;
            double p2xt2 = ax - rx + pareq[2] * x;
            double p3xt2 = rx - pareq[2] * x;
            double p1yt2 = ay + ry - pareq[3] * x;
            double p2yt2 = ay - ry + pareq[3] * x;
            double p3yt2 = ry - pareq[3] * x;
            double pexpxt2 = Math.exp(-Math.pow(p1xt2 / p0t2, 2)) + Math.exp(-Math.pow(p2xt2 / p0t2, 2))
                    - 2 * Math.exp(-Math.pow(p3xt2 / p0t2, 2));
            double perfxt2 = p1xt * Erf.erf(p1xt2 / p0t2) + p2xt2 * Erf.erf(p2xt2 / p0t2)
                    - 2 * p3xt2 * Erf.erf(p3xt2 / p0t2);
            double pexpyt2 = Math.exp(-Math.pow(p1yt2 / p0t2, 2)) + Math.exp(-Math.pow(p2yt2 / p0t2, 2))
                    - 2 * Math.exp(-Math.pow(p3yt2 / p0t2, 2));
            double perfyt2 = p1yt2 * Erf.erf(p1yt2 / p0t2) + p2yt2 * Erf.erf(p2yt2 / p0t2)
                    - 2 * p3yt2 * Erf.erf(p3yt2 / p0t2);

            double pplane2 = (p0t2 / sqrpi * pexpxt2 + perfxt2) * (p0t2 / sqrpi * pexpyt2 + perfyt2)
                    / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double pspim2 = 1 / Math.sqrt(1 + (4 * pareq[6] * x) / Math.pow(szeff, 2));
            double acf2 = pplane2 * pspim2;

            // COMPONENT 3
            // help variables, which are dependent on time, to write the full function
            double p0t3 = Math.sqrt(4 * pareq[8] * x + Math.pow(s, 2));
            double p1xt3 = ax + rx - pareq[2] * x;
            double p2xt3 = ax - rx + pareq[2] * x;
            double p3xt3 = rx - pareq[2] * x;
            double p1yt3 = ay + ry - pareq[3] * x;
            double p2yt3 = ay - ry + pareq[3] * x;
            double p3yt3 = ry - pareq[3] * x;
            double pexpxt3 = Math.exp(-Math.pow(p1xt3 / p0t3, 2)) + Math.exp(-Math.pow(p2xt3 / p0t3, 2))
                    - 2 * Math.exp(-Math.pow(p3xt3 / p0t3, 2));
            double perfxt3 = p1xt * Erf.erf(p1xt3 / p0t3) + p2xt3 * Erf.erf(p2xt3 / p0t3)
                    - 2 * p3xt3 * Erf.erf(p3xt3 / p0t3);
            double pexpyt3 = Math.exp(-Math.pow(p1yt3 / p0t3, 2)) + Math.exp(-Math.pow(p2yt3 / p0t3, 2))
                    - 2 * Math.exp(-Math.pow(p3yt3 / p0t3, 2));
            double perfyt3 = p1yt3 * Erf.erf(p1yt3 / p0t3) + p2yt3 * Erf.erf(p2yt3 / p0t3)
                    - 2 * p3yt3 * Erf.erf(p3yt3 / p0t3);

            double pplane3 = (p0t3 / sqrpi * pexpxt3 + perfxt3) * (p0t3 / sqrpi * pexpyt3 + perfyt3)
                    / (4 * Math.pow(ax * ay, 2) / fitobsvol);
            double pspim3 = 1 / Math.sqrt(1 + (4 * pareq[8] * x) / Math.pow(szeff, 2));
            double acf3 = pplane3 * pspim3;

            // TRIPLET
            double triplet = 1 + pareq[9] / (1 - pareq[9]) * Math.exp(-x / pareq[10]);

            return (1 / pareq[0])
                    * ((1 - pareq[5] - pareq[7]) * acf1 + Math.pow(q2, 2) * pareq[5] * acf2
                            + Math.pow(q3, 2) * pareq[7] * acf3)
                    / Math.pow(1 - pareq[5] - pareq[7] + q2 * pareq[5] + q3 * pareq[7], 2) * triplet + pareq[4];
        }
    }

    // generalized least square fit; uses FCS_3p, FCS_SX_3p, and FCCS_2p to tranform
    // the data via the covariance matrix
    class GLS_fitFunction implements ParametricUnivariateFunction {

        @Override
        public double[] gradient(double x, double[] params) {
            int num = params.length; // number of values in the gradient of FCS_3p, which depends on the number of
                                     // parameters to be fit
            double[][] gradtau = new double[num][chanum - 1];
            double[] tmpgrad = new double[num]; // temporary variable;
            double[] finalgrad = new double[num];
            int sol = 0; // get the index of the solution for the particular tau

            for (int g = 1; g < chanum; g++) { // determine which element is required by checking the lagtime
                if (lagtime[g] == x) {
                    sol = g - 1;
                }
            }

            if (x == lagtime[1]) {
                RealVector[] solution = new RealVector[num]; // vector containing solutions
                ParametricUnivariateFunction function;
                if ("FCS".equals(fitModel)) { // select the fit model to be used; extra fit models can be added here
                    function = new FCS_3p();
                } else {
                    function = new FCS_3p();
                }

                // else {
                // function = new FCCS_2p();
                // }
                for (int y = 1; y < chanum; y++) { // get all gradients for all lag times
                    tmpgrad = function.gradient(lagtime[y], params);
                    for (int z = 0; z < num; z++) {
                        gradtau[z][y - 1] = tmpgrad[z];
                    }
                }

                for (int z = 0; z < num; z++) { // solve for a new correlation vector with independent elements
                    DecompositionSolver solver = new LUDecomposition(lowerCholDecCovmats).getSolver();
                    RealVector constants = new ArrayRealVector(gradtau[z]);
                    solution[z] = solver.solve(constants);
                    transTheoreticalGradientACF[z] = solution[z].toArray(); // remember the transformed theoretical
                                                                            // function
                }

                for (int z = 0; z < num; z++) {
                    finalgrad[z] = transTheoreticalGradientACF[z][sol];
                }
            } else {
                for (int g = 2; g < chanum; g++) { // determine which element is required by checking the lagtime
                    if (lagtime[g] == x) {
                        sol = g - 1;
                    }
                }
                for (int z = 0; z < num; z++) {
                    finalgrad[z] = transTheoreticalGradientACF[z][sol];
                }
            }

            return finalgrad;
        }

        @Override
        public double value(double x, double[] params) {
            double[] valtau = new double[chanum - 1]; // array for theoretical ACF before transormation
            double retval; // return value
            int sol = 0; // index of the solution for a particular tau

            if (x == lagtime[1]) {
                ParametricUnivariateFunction function;
                if ("FCS".equals((String) fitModel)) { // select the fit model to be used; extra fit models can be added
                                                       // here
                    function = new FCS_3p();
                } else {
                    function = new FCS_3p();
                }

                // else {
                // function = new FCCS_2p();
                // }
                // calculate the correlation function for this particular set of parameters; do
                // not take the zero lagtime into account
                for (int y = 0; y < chanum - 1; y++) {
                    valtau[y] = function.value(lagtime[y + 1], params);
                }

                // use the regularized covariance matrix to transform the data
                DecompositionSolver solver = new LUDecomposition(lowerCholDecCovmats).getSolver(); // solve for a new
                                                                                                   // correlation vector
                                                                                                   // with independent
                                                                                                   // elements
                RealVector constants = new ArrayRealVector(valtau);
                RealVector solution = solver.solve(constants);

                for (int y = 0; y < chanum - 1; y++) { // remember the transformed theoretical function
                    transTheoreticalACF[y] = solution.getEntry(y);
                }

                retval = transTheoreticalACF[sol];

            } else {
                for (int g = 2; g < chanum; g++) {
                    if (lagtime[g] == x) {
                        sol = g - 1;
                    }
                }
                retval = transTheoreticalACF[sol];
            }

            return retval;
        }
    }

    /*
     * Miscellaneous functions
     * 
     * public double determinant(double A[][], int N): calculate determinant public
     * double[][] generateSubArray (double A[][], int N, int j1): generate a
     * subarray
     * 
     */
    // calculate determinant of a subarray; required for the Bayes model
    // probabilities
    // this was adapted from
    // http://stackoverflow.com/questions/16602350/calculating-matrix-determinant
    public double determinant(double A[][], int N) {
        double res;
        double[][] m;

        // 1x1 matrix
        if (N == 1) {
            res = A[0][0];
        } // 2x2 matrix
        else if (N == 2) {
            res = A[0][0] * A[1][1] - A[1][0] * A[0][1];
        } // NxN matrix
        else {
            res = 0;
            for (int j1 = 0; j1 < N; j1++) {
                m = generateSubArray(A, N, j1);
                res += Math.pow(-1.0, 1.0 + j1 + 1.0) * A[0][j1] * determinant(m, N - 1);
            }
        }
        return res;
    }

    public double[][] generateSubArray(double A[][], int N, int j1) {
        double[][] m = new double[N - 1][];
        for (int k = 0; k < (N - 1); k++) {
            m[k] = new double[N - 1];
        }
        for (int i = 1; i < N; i++) {
            int j2 = 0;
            for (int j = 0; j < N; j++) {
                if (j == j1) {
                    continue;
                }
                m[i - 1][j2] = A[i][j];
                j2++;
            }
        }
        return m;
    }

    // Dan
    public static class storeTempFit { // what is different from ImFCS impementation is it takes full double fit value
                                       // instead of 2 decimal place
        private static double[] tempfitres;

        public static void init(int a) {
            tempfitres = new double[a];
        }

        public static void storeTempFitres(int ind, double value) {
            tempfitres[ind] = value;
        }

        public static double getTempFitres(int ind) {
            return tempfitres[ind];
        }
    }

}
