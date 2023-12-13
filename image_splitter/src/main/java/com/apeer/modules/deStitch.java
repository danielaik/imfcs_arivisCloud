package com.apeer.modules;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.process.ShortProcessor;

/**
 * Imagestack/Image splitter
 * 
 * @author Daniel Aik
 * @version 1.0
 * @date 19/07/2020
 */

public class deStitch {

    public Map<String, Object> run(String inputImagePath, boolean isSplitLeftRight) {

        // var outputs = new HashMap<String, ImagePlus>();
        Map<String, Object> outputs = new HashMap<String, Object>();

        // input
        ImagePlus impOptosplit = new ImagePlus(inputImagePath);
        ImagePlus impaft_1;
        ImagePlus impaft_2;
        ImageStack imsaft_1;
        ImageStack imsaft_2;
        String $title;

        String file1;
        String file2;

        int widthaft;
        int heightaft;
        int sizeaft;

        $title = impOptosplit.getTitle();
        $title = $title.substring(0, $title.lastIndexOf('.'));

        boolean isLR = isSplitLeftRight; // split left-right or top-bottom

        /*
         * // might not allply universally to other camera manufacturer // automatically
         * findout if left-right or top-bottom split is necessary if
         * (impOptosplit.getImageStack().getWidth() <=
         * impOptosplit.getImageStack().getHeight()){ isLR = false; } else { isLR =
         * true; }
         */

        if (isLR) {
            widthaft = impOptosplit.getImageStack().getWidth() / 2;
            heightaft = impOptosplit.getImageStack().getHeight();
            sizeaft = impOptosplit.getStackSize();

        } else {
            widthaft = impOptosplit.getImageStack().getWidth();
            heightaft = impOptosplit.getImageStack().getHeight() / 2;
            sizeaft = impOptosplit.getStackSize();
        }

        if (isLR) {
            // destitch
            // left
            imsaft_1 = new ImageStack(widthaft, heightaft);
            for (int f = 1; f <= sizeaft; f++) {
                ShortProcessor ip = new ShortProcessor(widthaft, heightaft);
                for (int y = 0; y < heightaft; y++) {
                    for (int x = 0; x < widthaft; x++) {

                        ip.putPixel(x, y, impOptosplit.getImageStack().getProcessor(f).get(x, y));
                    }
                }
                imsaft_1.addSlice(ip);
            }

            // right
            // bottom
            imsaft_2 = new ImageStack(widthaft, heightaft);
            for (int f = 1; f <= sizeaft; f++) {
                ShortProcessor ip = new ShortProcessor(widthaft, heightaft);
                for (int y = 0; y < heightaft; y++) {
                    for (int x = 0; x < widthaft; x++) {

                        ip.putPixel(x, y, impOptosplit.getImageStack().getProcessor(f).get(x + widthaft, y));
                    }
                }
                imsaft_2.addSlice(ip);
            }
        } else {
            // destitch
            // top
            imsaft_1 = new ImageStack(widthaft, heightaft);
            for (int f = 1; f <= sizeaft; f++) {
                ShortProcessor ip = new ShortProcessor(widthaft, heightaft);
                for (int y = 0; y < heightaft; y++) {
                    for (int x = 0; x < widthaft; x++) {

                        ip.putPixel(x, y, impOptosplit.getImageStack().getProcessor(f).get(x, y));
                    }
                }
                imsaft_1.addSlice(ip);
            }

            // bottom
            imsaft_2 = new ImageStack(widthaft, heightaft);
            for (int f = 1; f <= sizeaft; f++) {
                ShortProcessor ip = new ShortProcessor(widthaft, heightaft);
                for (int y = 0; y < heightaft; y++) {
                    for (int x = 0; x < widthaft; x++) {
                        ip.putPixel(x, y, impOptosplit.getImageStack().getProcessor(f).get(x, y + heightaft));
                        // ip = null;
                    }
                }
                imsaft_2.addSlice(ip);
            }
        }

        impaft_1 = new ImagePlus("top", imsaft_1);
        impaft_2 = new ImagePlus("bottom", imsaft_2);
        impaft_1.resetDisplayRange();
        impaft_2.resetDisplayRange();

        if (isLR) {
            file1 = $title + "_left" + ".tiff";
            file2 = $title + "_right" + ".tiff";
        } else {
            file1 = $title + "_top" + ".tiff";
            file2 = $title + "_bottom" + ".tiff";
        }

        File outputFile = new File(file1);
        FileSaver fs = new FileSaver(impaft_1);
        if (impaft_1.getStackSize() > 1) {
            fs.saveAsTiffStack(outputFile.getAbsolutePath());
        } else {
            fs.saveAsTiff(outputFile.getAbsolutePath());
        }
        outputs.put("1", file1);

        File outputFile2 = new File(file2);
        FileSaver fs2 = new FileSaver(impaft_2);
        if (impaft_2.getStackSize() > 1) {
            fs2.saveAsTiffStack(outputFile2.getAbsolutePath());
        } else {
            fs2.saveAsTiff(outputFile2.getAbsolutePath());
        }
        outputs.put("2", file2);

        return outputs;

    }
}