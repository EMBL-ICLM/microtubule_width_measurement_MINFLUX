#@ Integer		(label = "line length (pixel)") lineLength 
#@ Integer		(label = "line width (pixel)", min=1, max=100) lineWidth
#@ Integer		(label = "angle resolution", value=18, min=2, max=360, style="slider") nLine
#@ Boolean 		(label = "Gaussian intensity", value=true) doGaussian
#@ Integer		(label = "minimum profile length (pixel)", value=550) minLength
#@ Boolean		(label = "show kernel", value=false) showKernel
#@ Boolean		(label = "show directional image", value=false) showConv
#@ Boolean		(label = "show tubeness image", value=false) showTube
#@ Boolean		(label = "show skeleton", value=false) showSkel
#@ Boolean		(label = "show line segments", value=false) showLineSeg

/*
 * 	Fiji groovy script to segment and create microtubule filament central line selections
 * 	need 3rd party plugin "MorphoLibJ", which can be obtained by enabling plugin update site: IJPB-plugins
 * 
 * 	It takes a rendered MINFLUX microtubule image as input.
 * 	1) then it apply 1st a diretional filter as user defined: 
 * 	which should be roughly the strucutring element of the microtubule filaments;
 * 	2) the curvilinear structures within the directional fitlered imgaes will be enhanced by "Tubeness" filter, 
 * 	which is similar to Frangi filter, based on computation with Hessian eigenvalues of the image;
 * 	3) the resulting image will then be thresholded with the Fiji default method to create the segmentation,
 * 	and further skeletonized into central lines of the microtubule ;
 * 	4) intersection part of the lines will be removed by morphological operations, to create single line segments.
 * 	ROIs delineate each line segments will be added to Fiji ROI Manager for further analysis or storage.
 * 	
 *  author: Ziqiang Huang <ziqiang.huang@embl.de>
 *  date: 2022.06.17
 */

import java.awt.*;
import java.awt.geom.Rectangle2D.Double;

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.*;
import ij.plugin.filter.ThresholdToSelection;
import ij.plugin.frame.RoiManager;

import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.algorithm.fft2.FFTConvolution;

import inra.ijpb.morphology.Morphology;
import inra.ijpb.morphology.Strel;

import features.TubenessProcessor;

	//	create a copy of the input image (within ROI if there exist) and conver to 8-bit
	ImagePlus imp_input = IJ.getImage();
	String name = imp_input.getTitle();
	double xBound = 0.0; double yBound = 0.0;
	if (null != imp_input.getRoi()) {
		xBound = imp_input.getRoi().getFloatBounds().getX();
		yBound = imp_input.getRoi().getFloatBounds().getY();
	}
	ImagePlus imp = convertToGray8(imp_input.crop());
	// dilate/external gradient original image to enlarge signal clusters
	Strel strel = Strel.Shape.DISK.fromRadius((int)(lineLength/100));	// create structuring element (cube of radius 'radius')
	imp = new ImagePlus(name+" dilated", Morphology.dilation(imp.getProcessor(), strel));

	// create line kernels for convolution
	ImagePlus imp_kernel = makeLineKernel (lineLength, lineWidth, nLine, doGaussian, showKernel);

	// compute line-filtered image
	ImageStack stack = new ImageStack(imp.getWidth(), imp.getHeight());
	for (int i=0; i<nLine; i++) {
		Img<FloatType> image = ImageJFunctions.convertFloat(imp);
		Img<FloatType> kernel = ImageJFunctions.convertFloat( new ImagePlus(""+i, imp_kernel.getStack().getProcessor(i+1)) );
		new FFTConvolution<>( image, kernel ).convolve();
		stack.addSlice(""+(1+i), image.getImagePlus().getProcessor() );
	}
	ImagePlus stk_conv = new ImagePlus(name + " convolution stack", stack);
	ImagePlus imp_conv = ZProjector.run(stk_conv, "max all" );
	if (showConv) {
		imp_conv.setTitle(name + " directional");
		stk_conv.show(); imp_conv.show();
	} else {
		stk_conv.close();
		IJ.run("Collect Garbage", "");
	}
	// tubeness of line-filtered image
	ImagePlus imp_tube = new TubenessProcessor(lineWidth, false).generateImage(imp_conv);
	if (showTube) {
		imp_tube.setTitle(name + "tubeness");
		imp_tube.show();
	}
	// generate skeleton from tubeness image 
	ImageProcessor skel = getSkeleton(imp_tube, "Default", name, showSkel);

	// trim intersections with morphological filter, and generate line segments as Roi[]
	Roi[] segments = getLineSegments(skel, 5, 15, xBound, yBound, name, showLineSeg);
	segments = filterSegments(segments, minLength);
	addToManager(segments);

	// clean up 
	imp.close();
	if(!showConv) imp_conv.close();
	if(!showTube) imp_tube.close();
	System.gc();
	IJ.run("Collect Garbage", "");
		
	// check input image, convert to 8-bit
    public ImagePlus convertToGray8 (ImagePlus imp) {
		if (imp.getType() != ImagePlus.GRAY8)
			IJ.run(imp, "8-bit", "");
    	return imp;
    }

	// check if input image is 3D 
	public boolean is3D (ImagePlus imp) {
		return imp.getNSlices()>1;
	}
	
    // create line profiles, with 180/nLine degree precision
    public ImagePlus makeLineKernel (int length, int width, int nLine, boolean doGaussian, boolean display) {
    	double angle = 180.0 / (double) nLine;
    	ImagePlus kernel = IJ.createImage("directionl kernel", "8-bit black", length, length, nLine);
    	ImageProcessor ip = makeRadialKernelProcessor(length, doGaussian);
    	double x_center = y_center = r = (double)length / 2;
    	for (int i=0; i<nLine; i++) {
    		double theta = angle * i;
    		double x_diff = r * Math.cos(theta * Math.PI / 180);
    		double y_diff = r * Math.sin(theta * Math.PI / 180);
    		double x1 = x_center + x_diff;	double y1 = y_center - y_diff;
    		double x2 = x_center - x_diff; 	double y2 = y_center + y_diff;
    		RotatedRectRoi roi = new RotatedRectRoi(x1, y1, x2, y2, (double)width);
    		ImageProcessor ip2 = ip.duplicate();
    		ip2.setValue(0);
    		ip2.fillOutside(roi);
    		kernel.getStack().setProcessor(ip2, i+1);
    	}
		if (display) kernel.show();
    	return kernel;
    }

	// make base radial symmetrical kernel image 
	public ImageProcessor makeRadialKernelProcessor(int diameter, boolean doGaussian) {
		ImageProcessor ip = NewImage.createByteImage("", diameter, diameter, 1, NewImage.FILL_BLACK).getProcessor();
		int Radius = (int) (diameter/2);
		if (doGaussian) {
			double sigma = (double) diameter / 6;
			double pxValueFactor = 255 / pdf(1, 0, sigma);
			for (int r=Radius; r>=1; r--) {
				ip.setValue( pxValueFactor * pdf(r, 0, sigma) );
				ip.fill( new OvalRoi(Radius-r, Radius-r, r*2, r*2) );
			}
		} else {
			ip.setValue(255);
			ip.fill( new OvalRoi(0, 0, diameter, diameter) );
		}
		return ip;
	}
	
    // return pdf(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
    public double pdf(double x, double mu, double sigma) {
        return pdf((x - mu) / sigma) / sigma;
    }
    	// return pdf(x) = standard Gaussian pdf
    public double pdf(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }


    // get skeleton from gray scale image 
    public ImageProcessor getSkeleton (ImagePlus imp, String method, String name, boolean display) {
    	imp.getProcessor().setAutoThreshold(method, true, ImageProcessor.NO_LUT_UPDATE);
		BinaryProcessor skel = new BinaryProcessor( imp.getProcessor().createMask() );
		skel.skeletonize(255);
		if (display)
			new ImagePlus(name + " skeleton", skel.duplicate()).show();
		return (ImageProcessor) skel;
    }

	// get line segments from skeleton image, as Roi[]
	public Roi[] getLineSegments(ImageProcessor ip, int blackTopDiam, int dilationDiam, double xShift, double yShift, String name, boolean display) {
		if (null == ip) return null;
		// enhance intersection regions from blackTopHat, and dilation operation // strel shape SQUARE
		Strel strel = Strel.Shape.DISK.fromRadius( blackTopDiam );
		ImageProcessor blackTop = Morphology.blackTopHat(ip, strel);
		if (0.0 != blackTop.getStats().max) {	// skip this step if there's no intersection in skeleton 
			Strel strel2 = Strel.Shape.DISK.fromDiameter( dilationDiam );
			ImageProcessor dilation = Morphology.dilation(blackTop, strel2);
			dilation.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE);
			// remove intersection regions from skeleton image
			ip.setValue(0);
			ip.fill( (Roi) new ThresholdToSelection().convert(dilation) );
			ip.setValue(255);
		}
		// create ROI selections from now separated skeleton image
		ip.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE);
		Roi roi = new ThresholdToSelection().convert(ip);
		if (null == roi) return null;
		roi.setLocation(roi.getXBase()+xShift, roi.getYBase()+yShift);
		if (display)
			new ImagePlus(name + " line segments", ip).show();
		if (Roi.COMPOSITE == roi.getType()) {	
			return ((ShapeRoi)roi).getRois();
		}
		return [roi];
	}

	// filter line segments by feret diameter 
	public Roi[] filterSegments (Roi[] rois, int minLength) {
		if (null == rois || 0 == rois.length) return null;
		int nROI = rois.length;
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i=0; i<nROI; i++) {
			//if (rois[i].getFeretsDiameter() >= minLength)
			//println("" + i + " : " + rois[i].getContainedPoints().length);
			if (rois[i].getContainedPoints().length >= minLength)
				list.add(i);
		}
		Roi[] selectedRois = new Roi[list.size()];
		for (int i=0; i<list.size(); i++) {
			selectedRois[i] = rois[list.get(i)];
		}
		return selectedRois;
	}
	
	// add ROIs to ROI Manager 
	public void addToManager(Roi[] rois) {
		if (null == rois || 0 == rois.length) return;
		RoiManager rm = RoiManager.getInstance();
		if (null == rm) rm = new RoiManager();
		else rm.reset();
		for (int i=0; i<rois.length; i++) {
			rm.addRoi(rois[i]);
			rm.rename(i, "segment "+(i+1));
		}
		return;
	}
