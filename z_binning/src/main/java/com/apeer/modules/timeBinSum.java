package com.apeer.modules;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.process.ShortProcessor;

/**
 * Z binning/ Frame binning
 * 
 * @author Daniel Aik
 * @version 1.0
 * @date 19/07/2030
 */

public class timeBinSum {

    public Map<String, Object> run(String inputImagePath, int shrinkFactorZ) {

        // var outputs = new HashMap<String, ImagePlus>();
        Map<String, Object> outputs = new HashMap<String, Object>();

        // input
        ImagePlus impbef = new ImagePlus(inputImagePath);
        String $title;
        ImagePlus impaft;
        ImageStack imsaft;
        int width = impbef.getImageStack().getWidth();
        int height = impbef.getImageStack().getHeight();
        int size = impbef.getImageStackSize();
        int sfz = shrinkFactorZ;
        boolean oddTermination = (size % sfz != 0); // return true if eg. 49999 total frame with shrink factor = 2.
                                                    // return false if total frame = 50,000 with shrink factor = 2
        int terminationframe = size - (size % sfz - 1) - (sfz - 1); // last frame to bin + 1

        $title = impbef.getTitle();
        $title = $title.substring(0, $title.lastIndexOf('.'));

        // do nothing if shrink factor = 1
        if (sfz == 1) {
            String fileName = $title + ".tiff";
            File outputFile = new File(fileName);
            FileSaver fs = new FileSaver(impbef);
            if (impbef.getStackSize() > 1) {
                fs.saveAsTiffStack(outputFile.getAbsolutePath());
            } else {
                fs.saveAsTiff(outputFile.getAbsolutePath());
            }
            outputs.put("timebinned", fileName);
            return outputs;
        }

        // process image and shrink images if shrink factor other than 1 is selected
        imsaft = new ImageStack(width, height);
        for (int f = 1; f <= size; f++) {
            ShortProcessor ip = new ShortProcessor(width, height);

            if (oddTermination && (f == terminationframe)) {
                break;
            }

            if (f % sfz == 1) {
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        int sum = 0;
                        for (int i = 0; i < sfz; i++) {
                            sum += impbef.getImageStack().getProcessor(f + i).get(x, y);
                        }
                        ip.putPixel(x, y, sum);
                    }
                }
                imsaft.addSlice(ip);
            }
        }

        impaft = new ImagePlus("impaft " + sfz, imsaft);
        impaft.resetDisplayRange();

        String fileName = $title + "_timebinned" + sfz + ".tiff";
        File outputFile = new File(fileName);
        FileSaver fs = new FileSaver(impaft);
        if (impaft.getStackSize() > 1) {
            fs.saveAsTiffStack(outputFile.getAbsolutePath());
        } else {
            fs.saveAsTiff(outputFile.getAbsolutePath());
        }
        outputs.put("timebinned", fileName);
        return outputs;

    }
}