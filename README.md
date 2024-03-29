# microtubule width measurement from MINFLUX data

## Workflow Description:

This workflow takes image rendered from localization data as input. With the image open (and active) in Fiji, a [1st Groovy script](/1_getLineSegments_.groovy) will apply a set of filters to highlight curvilinear structures within the image, and then extract skeletons, and eventually the line segments. The line segments will be displayed as ImageJ line-type ROIs, and stored in ROI Manager of Fiji.  
It's not always neccesary, at this stage, the user can modify, or fine tuning the line segments ROIs. The same as editing ROIs with ROI Manager, the user can add (means manual drawing), remove, or modify on the ROIs.  

Once the line segments are approved by the user, a [2nd Groovy script](/2_get_intensity_profile_plot_.groovy) can be applied to generate the cumulative intensity profile plot of the image. It takes the input image and the genreated line segment ROIs stored currently in the ROI Manager, scanning through all or selected line segments, and plot cumulative signal intensity against signal's perpendicular distances to the line segments.

##

## Demo:
prerequisite: Need Fiji and MorphoLibJ plugin. To install the MorphoLibJ plugin, in Fiji, go to Menu "Help" -> "Update...". After self-checking, within "ImageJ Updater" window, click "Manage update sites". Down the list, locate "IJPB-plugins" and tick the checkbox if it's not yet ticked. Close the "Manage update sites" window, "Apply changes" in "ImageJ Updater" window, and restart Fiji.

-  [sample image](/sample_data_rendered_with_4nm_pixel_size.tif), rendered from MINFLUX localization data as density map with 4nm pixel size. The Brightness and Contrast was tuned, and "Red Hot" Lookup table was applied in Fiji, to make visible the microtubule structures.  
note: this preview image was coverted to PNG format to be visible on Webpage. To run the script, both TIFF and PNG format works but we always use TIFF format as our input.  
    <p align="center">
    <img src="/preview/sample_data_preview.png" width="800" height=auto>
    </p>

 <br />
 
- with the sample image open in Fiji, run the Groovy script [1_getLineSegments_.groovy](/1_getLineSegments_.groovy). A dialog window will pop up to ask for a few input parameters:  
    <p align="center">
    <img src="/preview/script_input_demo_1.png" width="550" height=auto>
    </p>
    
 <br />    
-  intermediate filter results:  
    <p align="center">
    <img src="/preview/sample_data_filtered_image_preview.png" width="1100" height=auto>
    </p>  
    
    1st row left to right: 
       original input image;                directional filtered image;           Tuneness filtered image  
    
    2nd row left to right: 
       skeleton of tubeness image;          skeleton dilated;                     extracted line segments (overlay with ROIs)

 <br />
-  extracted line segments overlay onto sample data, as line-type ROIs stored in ROI Manager.  
    <p align="center">
    <img src="/preview/sample_data_with_line_segments_preview.png" width="1300" height=auto>
    </p>
    
    
 <br />
 
- with the now populated ROI Manager, together with the sample input image being the active image in Fiji, run the Groovy script [2_get_intensity_profile_plot_.groovy](/2_get_intensity_profile_plot_.groovy). Another dialog window will appear to ask for input parameters:  
note: the script will run over selected ROI entries within the ROI Manager. If nothing is selected, it runs over all ROI entries. The user can multi-select ROIs by holding down "ctrl" while clicking onto the ROI entries within the ROI Manager.
    <p align="center">
    <img src="/preview/script_input_demo_2.png" width="550" height=auto>
    </p>

<br />
- result: cumulative intensity profile plot:  
    <p align="center">
    <img src="/preview/script_output_demo_3.png" width="550" height=auto>
    </p>
    
<br />
- result: statistics of the intensity profile reported to Fiji Log Window:  

    Statistics of data sample_data_rendered_with_4nm_pixel_size.tif , sum
    filament length: 17861 pixels
    data:    signal intensity max:0.808, at: 50 ; min:0.016
              halfMax = 0.404
              FWHM = 6.550 pixels
    RMS:     signal intensity max:0.808, at: 50 ; min:0.016
              halfMax = 0.404
              FWHM = 6.550 pixels
