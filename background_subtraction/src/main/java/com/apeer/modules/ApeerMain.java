package com.apeer.modules;

import com.apeer.sdk.ApeerDevKit;
import com.apeer.sdk.ApeerEnvironmentException;
import com.apeer.sdk.ApeerInputException;
import com.apeer.sdk.ApeerOutputException;

/**
 * Background subtraction module
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
            var bgrnum = adk.getInput("background_num", Integer.class);
            var inputBckgPath = adk.getInput("input_background", String.class);

            var outputs = new BackgroundSubtraction().run(inputImagePath, bgrnum, inputBckgPath);

            adk.setFileOutput("output_image", (String) outputs.get("output_image"));

            adk.finalizeModule();

        } catch (ApeerEnvironmentException | ApeerInputException | ApeerOutputException e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }
    }
}
