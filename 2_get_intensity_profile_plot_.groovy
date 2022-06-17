//#@ ImagePlus (label = "ignore this") null
//#@ ImagePlus (label = "image to analyze") imp_input
#@ Integer	(label = "sampling distance to filament (px)", value = 150) maxDist
#@ Integer	(label = "RMS window size (one-sided)", value=5) rmsWinSize
#@ Boolean	(label = "show individual plot?", value = true) showEach
#@ Integer	(label = "shift speed (%)", value=95, min=0, max=100, style="slider") processSpeed
#@ Boolean	(label = "close all Plot window prior to processing?", value = true) closePlot

/*
 * 	Fiji groovy script to measure the width (Full-Width at Half Maximum, FWHM) of microtubule filaments
 * 	
 * 	It takes a rendered MINFLUX microtubule image, and line segments selection (in ROI Manager) as input.
 * 	1) The line segments selection can be generated with script: 1_getLineSegments_.groovy;
 * 	2) It plot the intensity profile along each line segments, against the distance to the central line,
 * 	similar to Fiji built-in "Plot Profile" function;
 * 	3) The scanning distance will be defined by the user, depending on the packing density of the microtubule filaments.
 * 	It's recommended to use a scanning distance that scans a region not including too much of other filament sigmals.
 * 	4) It plots also the RMS smoothed curve on top of the intensity profile plot. The RMS window size can be decided by the user.
 * 	5) The user can choose to display plot for each segment, or omit them. The scanning process can be visualized by choosing a 
 * 	processing speed less than 100. The less the processing speed, the slower the scanning process.
 * 	6) The measurement results will be reported to Fiji Log window.
 * 	
 *  author: Ziqiang Huang <ziqiang.huang@embl.de>
 *  date: 2022.06.17
 */

/*
 * generate intensity plot along line segment
 * 
 */

import ij.*;
import ij.gui.*;
import ij.plugin.*;
import ij.process.*;
import ij.plugin.filter.ThresholdToSelection;
import ij.plugin.frame.RoiManager;
import process3d.DistanceTransform3D;
import java.awt.Point;
import java.awt.Polygon;
	
		ImagePlus imp_input = IJ.getImage();
		if (imp_input.getWindow() instanceof PlotWindow) { IJ.error("Input image appear to be a Plot."); return; }
		String name = imp_input.getTitle();
		RoiManager rm = RoiManager.getInstance();
		if (null==rm || 0==rm.getCount()) { IJ.error("Need ROI Manager with microtuble segments as ROIs."); return; }
		Roi[] rois = rm.getSelectedRoisAsArray();
		ImagePlus imp = convertToGray8(imp_input.duplicate());
		if (closePlot) closeAllPlotWindow();
		
		
		double[] resultSum = new double[maxDist*2 + 1];
		int totalPoint = 0;
		
		for (int r=0; r<rois.length; r++) {
			int nRoiPoint = getRoiPointNum( rois[r] );
			totalPoint += nRoiPoint;
			double[] result_current = getDistanceShiftData(imp, rois[r], maxDist, processSpeed, imp_input);
			resultSum = elementWiseAdd(resultSum, result_current, nRoiPoint);
			if (showEach) {
				IJ.log("\n\n Statistics of data " + name + " , " +  rois[r].getName());
				IJ.log(" filament length: " + nRoiPoint + " pixels");
				//plotSegment(result_current, true, rois[r].getName());
				plotSegmentWithRMS(result_current, rmsWinSize, rois[r].getName());
			}
		}
	
		for (int i=-maxDist; i<=maxDist; i++) {
			resultSum[i+maxDist] /= totalPoint;
		}
		IJ.log("\n\n Statistics of data " + name + " , sum" );
		IJ.log(" filament length: " + totalPoint + " pixels");
		
		plotSegmentWithRMS(resultSum, rmsWinSize, name);
		
		// clean up 
		IJ.run("Collect Garbage", "");
		System.gc();


	// close all Plot Windows 
	public void closeAllPlotWindow () {
		int[] IDs = WindowManager.getIDList();
		for (int i=0; i<IDs.length; i++) {
			ImagePlus imp = WindowManager.getImage(IDs[i]);
			if (imp.getWindow() instanceof PlotWindow)
				imp.close();
		}
		IJ.run("Collect Garbage", "");
	}
	
	// check input image, convert to 8-bit
    public ImagePlus convertToGray8 (ImagePlus imp) {
		if (imp.getType() != ImagePlus.GRAY8)
			IJ.run(imp, "8-bit", "");
    	return imp;
    }
    
	// shift the ROI profile along the orthogonal of its long axis, and generate the distance shifted profile plot
	public double[] getDistanceShiftData (ImagePlus imp, Roi roi_in, int maxDist, int processSpeed, ImagePlus imp_show) {
		Roi roi = roi_in.clone();
		ImageProcessor ip = imp.getProcessor();
		
		double ori_x = roi.getXBase();
		double ori_y = roi.getYBase();
		Point[]	points = roi.getContainedPoints();
		double[] slope =  getInterceptSlope (points);
		
		double dx = slope[0];
		double dy = slope[1];

		double[] result = new double[maxDist*2 + 1];
		long waitTime = (long) Math.round(550 * (100-processSpeed)  / (points.length));
		for (int i=-maxDist; i<=maxDist; i++) {
			roi.setLocation(ori_x+ i*dx, ori_y+ i*dy);
			if (processSpeed<100) {
				imp_show.setRoi(roi);
				Thread.sleep(waitTime);
			}
			Point[]	points_new = roi.getContainedPoints();
			double sum = 0.0;
			for (int j=0; j<points_new.length; j++) {
				sum += ip.getPixel( (int)points_new[j].x, (int)points_new[j].y );
			}
			result[i+maxDist] = (sum / points_new.length);
		}
		roi.setLocation(ori_x, ori_y);
		if (processSpeed<100) imp.setRoi(roi);
		return result;
	}
    
    // interpolate line ROI, so instead of the two ends, it contains all coordinates in the line
	public Roi interpolate_line_roi (ImagePlus imp, Roi roi) { // make sure input is "Straight Line" ROI
		FloatPolygon poly = roi.getFloatPoints();
		double x1 = poly.xpoints[0];
		double x2 = poly.xpoints[1];
		double y1 = poly.ypoints[0];
		double y2 = poly.ypoints[1];
		double x_center = (x1 + x2) / 2;
		double y_center = (y1 + y2) / 2;
		double target_c = roi.getStrokeWidth()/2;
		double dx = x2 - x1; double dy = y2 - y1;
		double c = Math.sqrt(dx * dx + dy * dy);
		double delta_x = target_c * dy / c;
		double delta_y = target_c * dx / c;
		double new_x1 = x_center - delta_x;
		double new_y1 = y_center + delta_y;
		double new_x2 = x_center + delta_x;
		double new_y2 = y_center - delta_y;
		Roi new_roi = new Line(new_x1, new_y1, new_x2, new_y2);
		new_roi.setStrokeWidth(1.0);
		imp.setRoi( new_roi );
		ImagePlus mask = new ImagePlus ("mask", imp.createRoiMask());
		IJ.run(mask, "Create Selection", "");
		return mask.getRoi();
	}

	// element-wise sum of two array with the same size
	public double[] elementWiseAdd (double[] data1, double[] data2, int n_data2) {
		if (data1.length != data2.length) return null;
		int N = data1.length;
		double[] sum = new double[N];
		for (int i=0; i<N; i++) {
			sum[i] = data1[i] + data2[i]*n_data2;
		}
		return sum;
	}

	// get number of contained Points in ROI, could be less in case of straight line or polygon
	public int getRoiPointNum (Roi roi) {
		return roi.getContainedPoints().length;
	}
	
	// compute the intercept slope of line profile, integrate along the entire line
	public double[] getInterceptSlope (Point[]	points) {
		int N = points.length;
		double[] x = new double[N]; double[] y = new double[N];
		double x_mean = 0.0; double y_mean = 0.0;
		double x_sum = 0.0; double y_sum = 0.0;
		double x_sum2 = 0.0; double xy_sum = 0.0;
		for (int i=0; i<N; i++) {
			x[i] = points[i].x; y[i] = points[i].y;
			x_sum += x[i]; y_sum += y[i];
			x_sum2 += x[i] * x[i];
			xy_sum += x[i] * y[i];
		}
		x_mean = x_sum / N; y_mean = y_sum / N;
		double a = 	xy_sum - (y_mean * x_sum);	// dx 
		double b = 	x_sum2 - (x_mean * x_sum);	// dy
		double c = Math.sqrt( (a*a) + (b*b) );
		return [-a/c, b/c];
	}
	
	// plot segment with RMS smoothing curve
	public void plotSegmentWithRMS (double[] distData, int windowSize, String name) {
		// compute statistics, max value, index, and FWHM from distance-shifted data
		double[] stat = getStatFast(distData, false); // 0:mean, 1:std, 2:max, 3:min
		int max_index = getElementIndex(distData, stat[2]);
		double data_FWHM = getFWHMfromHist(distData);
		// report to Log
		IJ.log("  data:    signal intensity max:" + IJ.d2s(stat[2], 3) + ", at: " + max_index +" ; min:" + IJ.d2s(stat[3],3) );
		IJ.log("              halfMax = " + IJ.d2s((stat[2]/2), 3) );
		IJ.log("              FWHM = " + IJ.d2s(data_FWHM, 3) + " pixels");
		// plot distance-shifted data
		String xLabel = "Distance (pixel)";
        String yLabel = "signal Intensity (avg)"
        int n = distData.length;
        double[] xValues = new double[n];
		for (int i=0; i<n; i++) { xValues[i] = i; }
        Plot plot = new Plot("Profile Plot - " + name, xLabel, yLabel);
        plot.setLineWidth(2);
		plot.add("SEPARATED_BAR", xValues, distData);
		
		// compute RMS of data
		double[] data_RMS = generateRMS_patch(distData, windowSize);
		double[] stat_RMS = getStatFast(data_RMS, false); // 0:mean, 1:std, 2:max, 3:min
		int max_index_RMS = getElementIndex(data_RMS, stat_RMS[2]);
		double data_FWHM_RMS = getFWHMfromHist(data_RMS);
		// report to Log
		IJ.log("  RMS:   signal intensity max:" + IJ.d2s(stat_RMS[2], 3) + ", at: " + max_index_RMS +" ; min:" + IJ.d2s(stat_RMS[3], 3) );
		IJ.log("              halfMax = " + IJ.d2s((stat_RMS[2]/2), 3) );
		IJ.log("              FWHM = " + IJ.d2s(data_FWHM_RMS, 3) + " pixels");
		// plot RMS smoothed data
		plot.setLineWidth(1); plot.setColor(java.awt.Color.RED);
		plot.add("line", xValues, data_RMS);	
        plot.show();
	}


	// get 1st appearance of element index in data array, that matching specific value
	public int getElementIndex (double[] data, double value) {
		for (int i=0; i<data.length; i++) {
			if (value == data[i]) return i;
		}
		return -1;
	}
	
	// compute mean, standard deviation, max and min value from 1D data
    public double[] getStatFast (double[] data, boolean correctSampleSize) {
		if (data == null) return null;
		int N = data.length;	// get size of data
		if (N == 0) return null;
		
		double value = 0; double sum = 0; double sum2 = 0; 
		double max = -Double.MAX_VALUE ; double min = Double.MAX_VALUE ;
		for (int i=0; i<N; i++) {
			value = data[i];
			sum += value;
			sum2 += value * value;
			max = Math.max(max, value);
			min = Math.min(min, value);
		}
		double mean = sum / N;
		if (N==1) return [mean, 0, max, min];
		double var = sum2 / N - mean * mean;
		if (correctSampleSize) var *= N / (N-1);
		return [mean, Math.sqrt(var), max, min];
	}
	
	// generate patched data to compute RMS
	public double[] generateRMS_patch(double[] data, int windowSize) {
	 	double[] data_patched = patchData(data, windowSize, "const");
	 	return generateRMS (data_patched, windowSize);
	}
	public double[] generateRMS (double[] data, int windowSize) {
		int N = data.length;
		int N_new = N - 2 * windowSize;
		int N_window = 2*windowSize + 1;
		double[] RMS = new double[N_new];
		double[] data2 = new double[N];
		for (int i=0; i<N; i++) { data2[i] = data[i] * data[i]; }
		double[] sum2 = new double[N_new]; 
		for (int i=0; i<N_new; i++) {
			for (int j=0; j<N_window; j++) {
				sum2[i] += data2[i+j];
			}
		}
		for (int i=0; i<N_new; i++) {
			RMS[i] = Math.sqrt( sum2[i] / N_window );
		}
		return RMS;
	}
	
	// patch data: mirror, zero, constant
	public double[] patchData(double[] data, int windowSize, String type) {
		int N = data.length;
		int size = N + 2 * windowSize;
		double[] data_patched = new double[size];
		for (int i=0; i<windowSize; i++) {
			switch (type) { 
				case "mirror":
					data_patched[i] = data[windowSize - 1 - i];
					data_patched[i + N + windowSize] = data[N - 1 - i];
				case "zero":
					data_patched[i] = 0;
					data_patched[i + N + windowSize] = 0;
				case "const":
					data_patched[i] = data[0];
					data_patched[i + N + windowSize] = data[N-1];
			}
		}
		for (int i=0; i<N; i++) {
			data_patched[i + windowSize] = data[i];
		}
		return data_patched;
	}
	
	// recover photon count data from histogram
	public double[] getStatFromHist (double[] hist, double rejectValue) {
		int N = hist.length;
		double sum = 0.0;
		double totalNum = 0.0;
		// calculate sum of data, and count of data point
		for (int i=0; i<N; i++) { 
			if (hist[i] < rejectValue) continue;
			sum += hist[i]*i;
			totalNum += hist[i];
		}
		// calculate mean of data
		double mean = sum / totalNum;
		// calculate variance of data
		double variance = 0.0;
		for (int i=0; i<N; i++) {
			if (hist[i] < rejectValue) continue;
			variance += (hist[i] * (i-mean) * (i-mean));
		}
		double std = Math.sqrt(variance/totalNum);
		//var = var / totalNum;
		return [mean, std, totalNum];
	}

	// compute FWHM directly from distance shift measurement data
	public double getFWHMfromHist(double[] data) {
		int N = data.length;
		// find max value and index
		double max = 0.0; int max_index = 0;
		for (int i=0; i<N; i++) {
			if (data[i] > max) {
				max = data[i];
				max_index = i;
			}
		}
		double halfMax = max / 2;
		// using half of max value to locate left and right bound
		double leftBound = rightBound = 0;
		for (int i=max_index; i>0; i--) {
			if (data[i] <= halfMax) {	// left bound located
				leftBound = (double) i + ( (halfMax - data[i]) / (data[i+1] - data[i]) );
				break;
			}
		}
		for (int i=max_index; i<N; i++) {
			if (data[i] <= halfMax) {	// left bound located
				rightBound = (double) i - ( (halfMax - data[i]) / (data[i-1] - data[i]) );
				break;
			}
		}
		return (rightBound - leftBound);
	}

	