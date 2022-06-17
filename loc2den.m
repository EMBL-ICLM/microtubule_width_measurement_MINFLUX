function loc2den(MINFLUX_mat_file_path, auto_separation, num_channel, sigmaCutoff, marginRatio, pixelSize)
    % MATLAB script to export MINFLUX data as density map, saved to tiff
    % multi-page image format
    %   The rendering is done same as in the python script (mfx2png.py) in paper:
    %       Schmidt, R., et al. MINFLUX nanometer-scale 3D imaging and microsecond-range tracking on a common fluorescence microscope. 
    %       Nat Commun 12, 1478 (2021). https://doi.org/10.1038/s41467-021-21652-z

    %% load MINFLUX raw data in MATLAB data format
    [filepath,name, ~] = fileparts(MINFLUX_mat_file_path);
    data = load(MINFLUX_mat_file_path);
    vld = data.vld;
    
    %
    if (numel(size(data.loc)) == 2) % in case only 1 iteration in data
        xy = data.loc(vld, :);
    else
        xy = squeeze(data.loc(vld, end, 1:2)); % select last iteration of MINFLUX data
    end
    x = xy(:, 1); y = xy(:, 2);
    % calculate range for each dimension:
    xMin = min(x); xMax = max(x); xRange = xMax-xMin; xMid = (xMin+xMax)/2;
    yMin = min(y); yMax = max(y); yRange = yMax-yMin; yMid = (yMin+yMax)/2;
    ratio = 0.5 * (marginRatio/100 + 1);
    xMin = xMid - ratio * xRange; xMax = xMid + ratio * xRange;
    yMin = yMid - ratio * yRange; yMax = yMid + ratio * yRange;
    % pixel / voxel resolution in nm
    resolution = pixelSize * 1e-9;              
    numX = length(xMin: resolution :xMax);    % X dimension
    numY = length(yMin: resolution :yMax);    % Y dimension
    
    % if data is two channel, separate them by dcr attribute
    if (num_channel == 2) 
        dcr = data.dcr(vld, end);
        cut1 = .35; cut2 = .4;
        if (auto_separation)
            start_val = ones (numel(dcr), 1);
            start_val(dcr>0.4, :) = 2;
            model = fitgmdist(dcr, num_channel, 'Start', start_val);
            mu1 = model.mu(1);  sigma1 = sqrt(model.Sigma(:,:,1));
            mu2 = model.mu(2);  sigma2 = sqrt(model.Sigma(:,:,2));
            if (mu1 < mu2)
                cut1 = mu1 + sigmaCutoff*sigma1;
                cut2 = mu2 - sigmaCutoff*sigma2;
            else
                cut1 = mu2 + sigmaCutoff*sigma2;
                cut2 = mu1 - sigmaCutoff*sigma1;
            end  
        end
        
        fprintf("cut1: %d\ncut2: %d\n", cut1, cut2);
        ch1 = dcr>=0 & dcr<=cut1;
        ch2 = dcr<=1 & dcr>=cut2;
 
        denMap_1 = histcounts2(x(ch1), y(ch1), [xMin: resolution :xMax], [yMin: resolution :yMax]); %#ok<NBRAK> 
        denMap_1 = uint8(denMap_1);
        
        denMap_2 = histcounts2(x(ch2), y(ch2), [xMin: resolution :xMax], [yMin: resolution :yMax]); %#ok<NBRAK> 
        denMap_2 = uint8(denMap_2);

        % create save file path
        filename_tif1 = strcat(name, '_C1_render(', num2str(pixelSize), 'nm).tiff');
        filename_tif2 = strcat(name, '_C2_render(', num2str(pixelSize), 'nm).tiff');
        pathSave1 = fullfile(filepath, filename_tif1);
        pathSave2 = fullfile(filepath, filename_tif2);
        t1 = Tiff(pathSave1, 'w');
        t2 = Tiff(pathSave2, 'w');

    elseif (num_channel == 1)
        denMap = histcounts2(x, y, [xMin: resolution :xMax], [yMin: resolution :yMax]); %#ok<NBRAK> 
        denMap = uint8(denMap);
        % create save file path
        filename_tif = strcat(name, '_render(', num2str(pixelSize), 'nm).tiff');
        pathSave = fullfile(filepath, filename_tif);
        t = Tiff(pathSave, 'w');
    else
        error(" Channel number wrong! ");
    end

    %% populate TIFF property structure, same for 2D/3D image
    tagstruct.ImageLength = numX;
    tagstruct.ImageWidth = numY;
    %tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Centimeter;
    %tagstruct.XResolution = pixelSize * 1e-7;
    %tagstruct.YResolution = pixelSize * 1e-7;
    %tagstruct.ZResolution = pixelSize * 1e-7;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 8;
    tagstruct.SampleFormat = 1;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    
    %% Create a progress bar to monitor the exporting process
    if (num_channel == 2)
        setTag(t1, tagstruct);
        setTag(t2, tagstruct);
        write(t1, denMap_1(:, :));
        write(t2, denMap_2(:, :));
        writeDirectory(t1);
        writeDirectory(t2);
    elseif (num_channel == 1)
        setTag(t, tagstruct);
        write(t, denMap(:, :));
        writeDirectory(t);
    else

    end

    
    if (num_channel == 2)
        close(t1);
        close(t2);
    elseif (num_channel ==1)
        close(t);
    else

    end

    clear;
 
end

