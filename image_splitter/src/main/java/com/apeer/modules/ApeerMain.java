package com.apeer.modules;

import com.apeer.sdk.ApeerDevKit;
import com.apeer.sdk.ApeerEnvironmentException;
import com.apeer.sdk.ApeerInputException;
import com.apeer.sdk.ApeerOutputException;

/**
 * Imagestack/Image splitter
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
            var inputisLeftRight = adk.getInput("split_left-right", Boolean.class);
            var outputs = new deStitch().run(inputImagePath, inputisLeftRight);

            adk.setFileOutput("output_1", (String) outputs.get("1"));
            adk.setFileOutput("output_2", (String) outputs.get("2"));
            adk.finalizeModule();

        } catch (ApeerEnvironmentException | ApeerInputException | ApeerOutputException e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }
    }
}
