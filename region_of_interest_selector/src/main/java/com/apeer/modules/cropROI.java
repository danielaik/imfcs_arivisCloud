package com.apeer.modules;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.process.ShortProcessor;
import ij.process.FloatProcessor;

/**
 * Region of Interest Selector for imagestacks
 * 
 * @author Daniel Aik
 * @version 1.0
 * @date 19/07/2020
 */

// takes in 16-bit or 32-bit and export bithdepth equivalent to input image
public class cropROI {

    public Map<String, Object> run(String inputImagePath, int w, int h, int l, int t) {

        // var outputs = new HashMap<String, ImagePlus>();
        Map<String, Object> outputs = new HashMap<String, Object>();

        // setting max to 40 pixels. max pixel dimension for polynomial bleaach module
        // before exceed maximum heap space (REMOVED; restriction handled over to
        // polynomial bleach correction module)
        final int maxwidth;
        final int maxheight;
        String $title;

        // input
        ImagePlus impbef = new ImagePlus(inputImagePath);
        $title = impbef.getTitle();
        $title = $title.substring(0, $title.lastIndexOf('.'));
        ImagePlus impaft;
        int width_ori = impbef.getWidth();
        int height_ori = impbef.getHeight();
        int width;
        int height;
        int size = impbef.getImageStackSize();
        boolean is16Bit;

        if (impbef.getBitDepth() == 16) {
            is16Bit = true;
        } else {
            is16Bit = false;
        }

        // Attempting to make modules universal eg. not limited by 40 x 40 as required
        // by polynom bleacgh correction module
        // to limit replace width_ori and height_ori with the upper limit
        maxwidth = width_ori;
        maxheight = height_ori;

        if (width_ori <= maxwidth) {
            if (w > width_ori) {
                width = width_ori;
            } else {
                width = w;
            }
        } else {

            if (w > maxwidth) {
                width = maxwidth;
            } else {
                width = w;
            }
        }

        if (height_ori <= maxheight) {
            if (h > height_ori) {
                height = height_ori;
            } else {
                height = h;
            }

        } else {

            if (h > maxheight) {
                height = maxheight;
            } else {
                height = h;
            }

        }

        // shift
        int xshift;
        int yshift;

        if (l + width - 1 > width_ori) {
            xshift = width_ori - width;
        } else {
            xshift = l - 1;
        }

        if (t + height - 1 > height_ori) {
            yshift = height_ori - height;
        } else {
            yshift = t - 1;
        }

        ImageStack imsaft = new ImageStack(width, height);

        if (is16Bit) {
            for (int f = 1; f <= size; f++) {
                ShortProcessor ip = new ShortProcessor(width, height);
                for (int y = 0; y < height; y++) {
                    int yaft = y + yshift;
                    for (int x = 0; x < width; x++) {
                        int xaft = x + xshift;
                        ip.putPixel(x, y, impbef.getImageStack().getProcessor(f).get(xaft, yaft));
                    }
                }
                imsaft.addSlice(ip);
            }
        } else {
            for (int f = 1; f <= size; f++) {
                FloatProcessor ip = new FloatProcessor(width, height);
                for (int y = 0; y < height; y++) {
                    int yaft = y + yshift;
                    for (int x = 0; x < width; x++) {
                        int xaft = x + xshift;
                        ip.setf(x, y, impbef.getImageStack().getProcessor(f).getf(xaft, yaft));
                    }
                }
                imsaft.addSlice(ip);
            }

        }

        impaft = new ImagePlus("imp aft", imsaft);
        impaft.resetDisplayRange();

        String fileName = $title + "_croppedW" + width + "H" + height + "L" + l + "T" + t + ".tiff";

        File outputFile = new File(fileName);
        FileSaver fs = new FileSaver(impaft);
        if (impaft.getStackSize() > 1) {
            fs.saveAsTiffStack(outputFile.getAbsolutePath());
        } else {
            fs.saveAsTiff(outputFile.getAbsolutePath());
        }
        outputs.put("cropped", fileName);
        return outputs;

    }
}