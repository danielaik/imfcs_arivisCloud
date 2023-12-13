package com.apeer.modules;

import java.util.Formatter;
import java.util.HashMap;
import java.util.Map;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import java.io.File;
import java.io.FileWriter;
import ij.io.FileSaver;
import ij.process.FloatProcessor;
//import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collection;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.AbstractCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.linear.DiagonalMatrix;

/**
 * Bleach correction takes in 16-bit or 32-bit and export a 32-bit
 * 
 * @author Daniel Aik
 * @version 1.0
 * @date 19/07/2020
 */
public class PolynomialBleachCorrection {

    private double frametime; // acquisition time per frame
    private int roi1WidthX = 0; // dimensions of the ROI used for correlations
    private int roi1HeightY = 0;
    private int width; // width of loaded stack
    private int height; // height of loaded stack
    private int cfXDistance; // x,y distance between pixels to be correlated; in DC-FCCS this is the distance
                             // between corresponding pixels in red and green channels
    private int cfYDistance;
    private int binningX = 1;
    private int binningY = 1;
    private int roi1StartX = 0; // first pixel of the ROI used for correlation
    private int roi1StartY = 0;
    private int firstframe; // first frame to be used in the correlation
    private int lastframe; // last frame to be used in the correlation
    private int frames; // number of frames in stack
    private int polyOrder = 0; // polynomial order for bleach correction
    private float[][] filterArray; // filter array; intensity filters to select the pixels to be correlated
    private int pixbinX = 1; // this is a multiplication factor to determine positions in the image window;
                             // it is 1 for overlap and equal to binning for non-overlap
    private int pixbinY = 1;
    private String fitModel = "FCS";
    private String BleachCor = "Polynomial";
    private int nopit = 1; // points in the shortened intensity traces for plotting
    private double[] intTrace1; // intensity traces
    private double[] intTime; // time for intensity traces
    private float background; // background is determined from the smallest count value from the stack; this
                              // can be maually changed in the control panel
    private double[][] datac; // temporary array which stores the values along two pixel columns to be
                              // correlated through the stack
    private float impmin; // minimum value in the stack

    // bleach corrected array
    static ArrayList<Double> bc1Darray;
    ImageStack ims_bc;
    ImagePlus imp_bc;

    // input file
    public ImagePlus imp;

    // logger variable and methods
    public Formatter x;
    public FileWriter fw;
    public File nf;

    private boolean is16Bit;
    private final int maxWidth = 40;
    private final int maxHeight = 40;

    public void openFile() {
        try {
            x = new Formatter("logger.txt");
        } catch (Exception e) {
            System.out.println("openFile error");
        }
    }

    public void addRecords(String inputfilepath, int polydeg) {
        x.format("%s", "this is a log file update 2");
        x.format("\n%s \n%s", inputfilepath, polydeg);

    }

    public void closeFile() {
        x.close();
    }

    public Map<String, Object> log(String inputImagepath, int polydeg) {
        var out = new HashMap<String, Object>();
        nf = new File("test.csv");
        String b = nf.getAbsolutePath();
        String c = nf.getPath();
        openFile();
        addRecords(inputImagepath, polydeg);
        closeFile();

        out.put("output_file", "logger.txt");
        return out;

    }

    public Map<String, Object> run(String inputImagePath, int polydeg) {

        var outputs = new HashMap<String, Object>();

        // reading tiff file
        imp = new ImagePlus(inputImagePath);
        String $title;

        // get width, height, number of frames in stack, and magnification of the image
        // window
        width = imp.getWidth();
        height = imp.getHeight();
        frames = imp.getStackSize();
        firstframe = 1;
        lastframe = imp.getImageStackSize();
        frametime = 0.001;
        polyOrder = polydeg;
        impmin = 0;// backgdoun set to zero//minDetermination(imp); // calculate the minimum of the
                   // image stack; this will
                   // be used as the default
                   // background value
        background = impmin; // setting background to lowest value of all pixels before performing polynomial
                             // bleach correction

        // limiting dimension to 40 x 40
        if (width > maxWidth) {
            width = maxWidth;
        }
        if (height > maxHeight) {
            height = maxHeight;
        }

        if (imp.getBitDepth() == 16) {
            is16Bit = true;
            IJ.log("it is a 16-bit iamge");
        } else {
            is16Bit = false;
            IJ.log("it is a 32-bit image");
        }

        $title = imp.getTitle();
        $title = $title.substring(0, $title.lastIndexOf('.'));

        if (setParameters()) {

            // overlap = false
            roi1WidthX = (int) Math.floor((width - Math.abs(cfXDistance)) / binningX) * binningX;
            roi1HeightY = (int) Math.floor((height - Math.abs(cfYDistance)) / binningY) * binningY;

            Roi impRoi1 = new Roi(roi1StartX, roi1StartY, roi1WidthX, roi1HeightY);
            imp.setRoi(impRoi1);

            correlateROI(impRoi1);

            imp_bc = new ImagePlus("bc imageplus",
                    getBleachCorrectStack(firstframe, lastframe, width, height, bc1Darray));

        }

        String fileName = $title + "_BleachCorrected" + ".tiff";
        File outputFile = new File(fileName);

        // write tiff
        FileSaver fs = new FileSaver(imp_bc);
        if (imp_bc.getStackSize() > 1) {
            fs.saveAsTiffStack(outputFile.getAbsolutePath());
        } else {
            fs.saveAsTiff(outputFile.getAbsolutePath());
        }

        outputs.put("output_image", fileName);

        return outputs;

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

        return min;
    }

    public boolean setParameters() {

        // set common arrays and parameters according to user settings in the panel
        if ((lastframe - firstframe + 1) < 1000) { // use 1000 points for the intensity, except when less than 1000
                                                   // frames are present
            nopit = (lastframe - firstframe + 1);
        } else {
            nopit = 1000;
        }

        // initialize arrays required for calculations; they change with each new
        // paramter setting and are thus re-initialized
        intTrace1 = new double[nopit];
        intTime = new double[nopit];

        int totalpixel = frames * width * height;
        bc1Darray = new ArrayList<Double>(totalpixel);
        return true;

    }

    // Calculate correlations for a ROI
    public void correlateROI(Roi improi) {

        // CPU calculation

        /*
         * //TODO error importing java.awt.*; Rectangle imprect = improi.getBounds();
         * int startXmap = (int) Math.ceil(imprect.getX() / pixbinX); int startYmap =
         * (int) Math.ceil(imprect.getY() / pixbinY); int endXmap = (int)
         * Math.floor((imprect.getX() + imprect.getWidth() - binningX) / pixbinX); int
         * endYmap = (int) Math.floor((imprect.getY() + imprect.getHeight() - binningY)
         * / pixbinY);
         */

        int startXmap = (int) Math.ceil(0 / pixbinX);
        int startYmap = (int) Math.ceil(0 / pixbinY);
        int endXmap = (int) Math.floor((0 + width - binningX) / pixbinX);
        int endYmap = (int) Math.floor((0 + height - binningY) / pixbinY);
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
            for (int x = startXmap; x <= endXmap; x++) {
                for (int y = startYmap; y <= endYmap; y++) {
                    if (!Float.isNaN(filterArray[x * pixbinX][y * pixbinY])) {
                        calcIntensityTrace(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance,
                                y * pixbinY + cfYDistance, firstframe, lastframe);
                        correlate(imp, x * pixbinX, y * pixbinY, x * pixbinX + cfXDistance, y * pixbinY + cfYDistance,
                                0, firstframe, lastframe);
                    }
                    IJ.showProgress(x - startXmap, startXmap - endXmap);
                }
            }
        }

    }

    // calculate reduced intensity traces for plotting average and use not more than
    // 'nopit' (defined in setParameters()) points for a trace
    public void calcIntensityTrace(ImagePlus image, int ipx1, int ipy1, int ipx2, int ipy2, int initialframe,
            int finalframe) {

        // image: imp form which intensity will be taken
        // px1, py1, px2, py2: coordinates of pixels to be correlated
        // initialframe and finalframe provide the range of frames to be used
        int ave = (int) Math.floor((finalframe - initialframe + 1) / nopit); // calculate number of data points which
                                                                             // are averaged
        int sum1;
        int sum2;
        float bckg1 = background;
        float bckg2 = background; // needs to be adapted for background2 once available

        for (int x = 0; x < nopit; x++) {
            sum1 = 0; // initialize arrays with 0
            sum2 = 0;
            for (int i = 0; i < binningX; i++) {
                for (int k = 0; k < binningY; k++) {
                    for (int y = initialframe + x * ave; y <= initialframe + (x + 1) * ave - 1; y++) {

                        sum1 += image.getStack().getProcessor(y).getf(ipx1 + i, ipy1 + k) - bckg1;
                        sum2 += image.getStack().getProcessor(y).getf(ipx2 + i, ipy2 + k) - bckg2;
                    }
                }
            }

            intTime[x] = frametime * (x + 0.5) * ave;
            intTrace1[x] = sum1 / ave; // calculate average intensity for the 'ave' points
        }

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
        int pxm1; // pixel coordinates on the binned grid used to store the output and map it to
                  // the parameter map
        int pym1;

        pxm1 = (int) px1 / pixbinX; // map to the pixel on the pixel on the binned grid
        pym1 = (int) py1 / pixbinY;

        // Sliding window is not selected, correlate the full intensity trace
        datac = new double[2][num + 1]; // get the intensity data for the correlation
        datac[0] = getIntensity(image, px1, py1, 1, initialframe, finalframe); // getIntensity for first pixel; performs
                                                                               // a bleach correction if indicated in
                                                                               // the panel
        datac[1] = datac[0];
        fill1Dbc(bc1Darray, datac, frames);

    }

    // filling bleach corrected 1D array of double
    public void fill1Dbc(ArrayList<Double> fulllist, double[][] sublist, int noframes) {
        for (int f = 1; f <= noframes; f++) {
            fulllist.add(sublist[0][f]);
        }
    }

    // recovering bleach corrected 1Darray to imagestack
    public static ImageStack getBleachCorrectStack(int ff, int lf, int width, int height, ArrayList<Double> fulllist) {
        int num = (lf - ff + 1);
        ImageStack ims = new ImageStack();

        for (int f = 1; f <= num; f++) {
            FloatProcessor ip = new FloatProcessor(width, height);
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    int index = (f - 1) + (x * height * num) + (y * num);

                    ip.putPixelValue(x, y, fulllist.get(index));
                }
            }
            ims.addSlice(ip);
        }
        return ims;
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
                    intdat[x] += image.getStack().getProcessor(initialframe + x - 1).getf(px + i, py + k) - bckg;
                }
            }

        }

        String $bcmode = (String) BleachCor; // perform single or double exponential bleach corrections if selected

        if ("Polynomial".equals($bcmode)) { // fitting with a polynomial of selected order
            PolynomFit polfit = new PolynomFit();
            double corfunc;
            int maxord = polyOrder;

            // note that the bleach correction is performed on the averaged intensity traces
            // to make it faster
            // while you may have 20,000 intensity points, intTrace1 and 2 contain only
            // 1,000 points
            // see definition in setParameters()
            res = polfit.doFit(intTrace1);

            for (int x = 1; x <= num; x++) {
                corfunc = 0;
                for (int i = 0; i <= maxord; i++) {
                    corfunc += res[i] * Math.pow(frametime * (x + 0.5), i);
                }
                intdat[x] = intdat[x] / Math.sqrt(corfunc / res[0]) + res[0] * (1 - Math.sqrt(corfunc / res[0]));
            }
            if (mode == 1) {
                for (int x = 0; x < nopit; x++) {
                    corfunc = 0;
                    for (int i = 0; i <= maxord; i++) {
                        corfunc += res[i] * Math.pow(intTime[x], i);
                    }
                    intTrace1[x] = intTrace1[x] / Math.sqrt(corfunc / res[0])
                            + res[0] * (1 - Math.sqrt(corfunc / res[0]));
                }
            }

        }

        return (intdat);
    }

    // polynomial for bleach correction
    class Polynomial implements ParametricUnivariateFunction {

        @Override
        public double[] gradient(double x, double[] params) {
            int maxord = polyOrder;

            double[] grad = new double[maxord + 1];
            for (int i = 0; i <= maxord; i++) {
                grad[i] = Math.pow(x, i);
            }

            return grad;
        }

        @Override
        public double value(double x, double[] params) {
            int maxord = polyOrder;
            double[] A = new double[maxord + 1];
            for (int i = 0; i <= maxord; i++) {
                A[i] = params[i];
            }

            double val = 0;
            for (int i = 0; i <= maxord; i++) {
                val += A[i] * Math.pow(x, i);
            }

            return val;
        }
    }

    // Polynomial fit for bleach correction
    public class PolynomFit extends AbstractCurveFitter {

        @Override
        protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> points) {
            final int len = points.size();
            final double[] target = new double[len];
            final double[] weights = new double[len];
            final double[] initialGuess = new double[polyOrder + 1];

            int i = 0;
            for (WeightedObservedPoint point : points) {
                target[i] = point.getY();
                weights[i] = point.getWeight();
                i += 1;
            }

            // initial guesses
            initialGuess[0] = target[len - 1]; // use the last point as offset estimate
            for (int j = 1; j <= polyOrder; j++) { // use a straight line as the first estimate
                initialGuess[j] = 0;
            }

            ParametricUnivariateFunction function;
            function = new Polynomial();

            final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(
                    function, points);

            return new LeastSquaresBuilder().maxEvaluations(Integer.MAX_VALUE).maxIterations(Integer.MAX_VALUE)
                    .start(initialGuess).target(target).weight(new DiagonalMatrix(weights))
                    .model(model.getModelFunction(), model.getModelFunctionJacobian()).build();
        }

        public double[] doFit(double[] itrace) {
            PolynomFit fitter = new PolynomFit();
            ArrayList<WeightedObservedPoint> points = new ArrayList<>();
            int num = itrace.length;

            // Add points here
            for (int i = 0; i < num; i++) {
                WeightedObservedPoint point = new WeightedObservedPoint(1, intTime[i], itrace[i]);
                points.add(point);
            }

            double result[] = fitter.fit(points);
            return (result);
        }
    }

}