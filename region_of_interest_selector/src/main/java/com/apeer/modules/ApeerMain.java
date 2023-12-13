package com.apeer.modules;

import com.apeer.sdk.ApeerDevKit;
import com.apeer.sdk.ApeerEnvironmentException;
import com.apeer.sdk.ApeerInputException;
import com.apeer.sdk.ApeerOutputException;

/**
 * Region of Interest Selector for imagestacks
 * 
 * @author Daniel Aik
 * @version 1.0
 * @date 19/07/2020
 */

public class ApeerMain {
    public static void main(String[] args) {
        try {
            var adk = new ApeerDevKit();
            var inputWidth = adk.getInput("input_width", Integer.class);
            var inputHeight = adk.getInput("input_height", Integer.class);
            var inputLeft = adk.getInput("input_left", Integer.class);
            var inputTop = adk.getInput("input_top", Integer.class);
            var inputImagePath = adk.getInput("input_image", String.class);

            var outputs = new cropROI().run(inputImagePath, inputWidth, inputHeight, inputLeft, inputTop);

            adk.setFileOutput("output_image", (String) outputs.get("cropped"));
            adk.finalizeModule();

        } catch (ApeerEnvironmentException | ApeerInputException | ApeerOutputException e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }
    }
}
