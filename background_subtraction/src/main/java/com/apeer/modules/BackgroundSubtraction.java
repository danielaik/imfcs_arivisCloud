package com.apeer.modules;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

/**
 * Background subtraction module
 * 
 * @author Daniel Aik
 * @version 1.0
 * @date 19/07/2020
 */

public class BackgroundSubtraction {

    // input file
    public ImagePlus imp;
    public ImagePlus impBGR = null;

    // background subtracted
    public ImageStack ims_bg;
    public ImagePlus imp_bg;

    ImageProcessor bgrCorip;
    double[][] bgrmean;
    int background;
    boolean isBGRavail = false;

    public Map<String, Object> run(String inputImagePath, int bgrnum, String inputBckgPath) {

        var outputs = new HashMap<String, Object>();

        // reading tiff file
        // Sampple file
        imp = new ImagePlus(inputImagePath);
        String $title;

        // background file
        impBGR = new ImagePlus(inputBckgPath);

        int size = imp.getImageStackSize();
        int width = imp.getStack().getWidth();
        int height = imp.getStack().getHeight();
        int sizeBGR = 0;
        int widthBGR = 0;
        int heightBGR = 0;

        background = bgrnum;

        $title = imp.getTitle();
        $title = $title.substring(0, $title.lastIndexOf('.'));

        // Find out if impBGR is imp? Will proceed to subtract from constant background
        // if impBGR == imp

        isBGRavail = !isSameImageStack(imp, impBGR);

        if (isBGRavail) {
            sizeBGR = impBGR.getStackSize();
            widthBGR = impBGR.getStack().getWidth();
            heightBGR = impBGR.getStack().getHeight();
        }

        if (background == 65536) {
            background = minDetermination(imp, size, width, height);
        }

        ims_bg = new ImageStack(width, height);
        if (!isBGRavail) {
            System.out.println("background file is NOT available, subtracting a constant background of " + background);
            for (int f = 1; f <= size; f++) {
                ShortProcessor ip = new ShortProcessor(width, height);
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        int value = imp.getStack().getProcessor(f).get(x, y) - background;
                        ip.putPixel(x, y, value);
                    }
                }
                ims_bg.addSlice(ip);
            }
        } else {
            System.out.println(
                    "background file is available, subracting each pixels with unique counts calculated from average background image");
            // calculate mean
            bgrmean = new double[width][height];
            for (int z = 1; z <= sizeBGR; z++) {
                bgrCorip = impBGR.getStack().getProcessor(z);
                for (int x = 0; x < widthBGR; x++) {
                    for (int y = 0; y < heightBGR; y++) {
                        bgrmean[x][y] += bgrCorip.getPixel(x, y);
                    }
                }
            }
            for (int x = 0; x < widthBGR; x++) { // caluclate mean and variance: E(x) = Sum(x)/n, E(x^2) = Sum(x^2)/n
                                                 // and Var = E(x^2) - E(x)^2
                for (int y = 0; y < heightBGR; y++) {
                    bgrmean[x][y] /= sizeBGR;
                }
            }

            // subtracting all pixels accross frame with mean bgr [][] from background file
            for (int f = 1; f <= size; f++) {
                ShortProcessor ip = new ShortProcessor(width, height);
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        background = (int) Math.round(bgrmean[x][y]);
                        int value = imp.getStack().getProcessor(f).get(x, y) - background;
                        ip.putPixel(x, y, value);
                    }
                }
                ims_bg.addSlice(ip);
            }
        }

        imp_bg = new ImagePlus("after subtrack", ims_bg);
        imp_bg.resetDisplayRange();

        // Saving bacgkround subtracted image stack
        String fileName = $title + "_BackgroundSubtracted" + ".tiff";
        File outputFile = new File(fileName);

        // write tiff
        FileSaver fs = new FileSaver(imp_bg);
        if (imp_bg.getStackSize() > 1) {
            fs.saveAsTiffStack(outputFile.getAbsolutePath());
        } else {
            fs.saveAsTiff(outputFile.getAbsolutePath());
        }

        outputs.put("output_image", fileName);

        return outputs;

    }

    private static int minDetermination(ImagePlus image, int frames, int width, int height) {
        int min;
        min = image.getStack().getProcessor(1).get(1, 1);
        for (int z = 1; z <= frames; z++) {
            for (int x = 1; x < width; x++) {
                for (int y = 1; y < height; y++) {
                    if (image.getStack().getProcessor(z).get(x, y) < min) {
                        min = image.getStack().getProcessor(z).get(x, y);
                    }
                }
            }
        }

        return min;
    }

    private static boolean isSameImageStack(ImagePlus image1, ImagePlus image2) {
        boolean isSame;

        ImagePlus imp_sample = image1;
        ImagePlus imp_background = image2;

        isSame = imp_sample.getStackSize() == imp_background.getStackSize();
        if (!isSame) {
            return false;
        }
        isSame = imp_sample.getImageStack().getHeight() == imp_sample.getImageStack().getHeight();
        if (!isSame) {
            return false;
        }
        isSame = imp_sample.getImageStack().getWidth() == imp_sample.getImageStack().getWidth();
        if (!isSame) {
            return false;
        }
        int size = imp_sample.getStackSize();
        int wid = imp_sample.getImageStack().getWidth();
        int hei = imp_sample.getImageStack().getHeight();
        int valuePixSample;
        int valuePixBackground;
        for (int z = 1; z <= size; z++) {
            for (int x = 0; x < wid; x++) {
                for (int y = 0; y < hei; y++) {
                    valuePixSample = imp_sample.getStack().getProcessor(z).get(x, y);
                    valuePixBackground = imp_background.getStack().getProcessor(z).get(x, y);
                    if (valuePixSample != valuePixBackground) {
                        isSame = false;
                        return isSame;
                    }
                }
            }
        }
        // retrun true if all condistion satisfied; true means two ImagePlus are
        // identical
        return true;

    }

}
