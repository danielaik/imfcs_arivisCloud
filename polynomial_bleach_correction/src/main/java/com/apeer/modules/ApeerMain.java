package com.apeer.modules;

import com.apeer.sdk.ApeerDevKit;
import com.apeer.sdk.ApeerEnvironmentException;
import com.apeer.sdk.ApeerInputException;
import com.apeer.sdk.ApeerOutputException;

/**
 * Bleach correction takes in 16-bit or 32-bit and export a 32-bit
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
            var polydeg = adk.getInput("polynomial_degree", Integer.class);

            var outputs = new PolynomialBleachCorrection().run(inputImagePath, polydeg);
            // var outputlog = new PolynomialBleachCorrection().log(inputImagePath,
            // polydeg);

            adk.setFileOutput("output_image", (String) outputs.get("output_image"));
            // adk.setFileOutput("output_log", (String) outputlog.get("output_file"));
            adk.finalizeModule();

        } catch (ApeerEnvironmentException | ApeerInputException | ApeerOutputException e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }
    }
}
