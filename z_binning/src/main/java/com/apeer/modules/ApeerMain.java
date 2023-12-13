package com.apeer.modules;

import com.apeer.sdk.ApeerDevKit;
import com.apeer.sdk.ApeerEnvironmentException;
import com.apeer.sdk.ApeerInputException;
import com.apeer.sdk.ApeerOutputException;

/**
 * Z binning/ Frame binning
 * 
 * @author Daniel Aik
 * @version 1.0
 * @date 19/07/2030
 */

public class ApeerMain {
    public static void main(String[] args) {
        try {
            var adk = new ApeerDevKit();
            var inputImagePath = adk.getInput("input_image", String.class);
            var inputShrinkFactor = adk.getInput("z_shrink_factor", Integer.class);

            var outputs = new timeBinSum().run(inputImagePath, inputShrinkFactor);

            adk.setFileOutput("output_image", (String) outputs.get("timebinned"));
            adk.finalizeModule();

        } catch (ApeerEnvironmentException | ApeerInputException | ApeerOutputException e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }
    }
}
