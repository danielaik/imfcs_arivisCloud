package com.apeer.modules;

import java.io.IOException;

import com.apeer.sdk.ApeerDevKit;
import com.apeer.sdk.ApeerEnvironmentException;
import com.apeer.sdk.ApeerInputException;
import com.apeer.sdk.ApeerOutputException;

/**
 * autocorrelations of image stacks and data fitting for Fluorescence
 * Correlation Spectroscopy (FCS)
 * 
 * @author Daniel Aik
 * @version 1.0
 * @date 19/07/2020
 */

public class ApeerMain {
    public static void main(String[] args) {
        try {
            var adk = new ApeerDevKit();
            var inputImagePath = adk.getInput("input_image", String.class);
            var inputBinX = adk.getInput("input_binx", Integer.class);
            var inputFrametime = adk.getInput("input_frametime", Double.class);
            var inputBinY = adk.getInput("input_biny", Integer.class);
            var inputp = adk.getInput("input_correlator_p", Integer.class);
            var inputq = adk.getInput("input_correlator_q", Integer.class);
            var inputisFit = adk.getInput("input_fit", Integer.class);
            var inputPixelSize = adk.getInput("input_pixelsize", Double.class);
            var inputMag = adk.getInput("input_magnification", Double.class);
            var inputNA = adk.getInput("input_na", Double.class);
            var inputEmLambda = adk.getInput("input_emlambda", Double.class);
            var inputisGLS = adk.getInput("input_gls", Integer.class);
            var inputisBayes = adk.getInput("input_bayes", Integer.class);
            var inputFilterLL = adk.getInput("filterll", Integer.class);
            var inputFilterUL = adk.getInput("filterul", Integer.class);

            MultiTauCor ImageObject1 = new MultiTauCor();
            var outputs = ImageObject1.run(inputImagePath, inputFrametime, inputBinX, inputBinY, inputp, inputq,
                    inputisFit, inputPixelSize, inputMag, inputNA, inputEmLambda, inputisGLS, inputisBayes,
                    inputFilterLL, inputFilterUL);
            // var outlogs = ImageObject1.log(inputImagePath, inputFrametime, inputBinX,
            // inputBinY, inputp, inputq);
            var outhtml = ImageObject1.getHtml();

            adk.setFileOutput("output_xlsx", (String) outputs.get("output_excel"));
            // adk.setFileOutput("output_log", (String) outlogs.get("output_file"));
            adk.setFileOutput("output_html", (String) outhtml.get("output_new"));

            adk.finalizeModule();

        } catch (ApeerEnvironmentException | ApeerInputException | ApeerOutputException | IOException e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }
    }
}
