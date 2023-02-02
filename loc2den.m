
    %% MATLAB script to export MINFLUX data as density map, saved to tiff multi-page image format
    %   The rendering is done same as in the python script (mfx2png.py) in paper:
    %       Schmidt, R., et al. MINFLUX nanometer-scale 3D imaging and microsecond-range tracking on a common fluorescence microscope. 
    %       Nat Commun 12, 1478 (2021). https://doi.org/10.1038/s41467-021-21652-z
    % last modified 2023.02.01
    % Ziqiang Huang <ziqiang.huang@embl.de>
    %%

    %% load data file
    [filename, filepath] = uigetfile({'*.mat'}, 'MINFLUX raw data file');
    if isequal(filename, 0)
        return;
    end

    %% get user input: render pixel size, margin ratio, channel settings
    prompt = {'pixel size (nm):',...
        'margin ratio (%):',...
        'channel number: (e.g.: 1 or 2)',...
        'auto channel separation: (e.g.: yes or no)',... 
        'sigma cut off: (e.g.: 0.5)'};
    dlgtitle = 'Input';
    dims = [1, 55];
    definput = {'4','5', '1', 'no', '0.5'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    if isempty(answer)
        return;
    end
    pixelSize = str2double(answer{1});
    marginRatio = str2double(answer{2});
    numChannel = uint8(str2double(answer{3}));
    autoChannelSep = ~isequal('n', lower(answer{4}(1))); % whether to perform auto channel separation: if first letter input is n (or N), then No.
    sigmaCutoff = str2double(answer{5});

    %% load MINFLUX raw data
    data = load(fullfile(filepath, filename));
    vld = data.vld;
    % get localization coordinates
    if (numel(size(data.loc)) == 2) % in case only 1 iteration in data
        xy = data.loc(vld, :);
    else
        xy = squeeze(data.loc(vld, end, 1:2)); % select last iteration of MINFLUX data
    end
    x = xy(:, 2); y = xy(:, 1); % transpose 1st and 2nd dimension of matrix
    % calculate range for each dimension:
    xMin = min(x); xMax = max(x); xRange = xMax-xMin; xMid = (xMin+xMax)/2;
    yMin = min(y); yMax = max(y); yRange = yMax-yMin; yMid = (yMin+yMax)/2;
    halfSize = 0.5 * (1 + marginRatio/100);
    xMin = xMid - halfSize * xRange; xMax = xMid + halfSize * xRange;
    yMin = yMid - halfSize * yRange; yMax = yMid + halfSize * yRange;
    % pixel/voxel resolution in nm
    resolution = pixelSize * 1e-9;   
    xGrid = xMin : resolution : xMax;   % X grid
    yGrid = yMin : resolution : yMax;   % Y grid
    
    % if data is two channel, separate them by dcr attribute
    if (numChannel == 2)
        dcr = data.dcr(vld, end);
        [cut1, cut2] = channelSeparation (dcr, autoChannelSep, sigmaCutoff);
        ch1 = dcr>=0 & dcr<=cut1;
        ch2 = dcr<=1 & dcr>=cut2;
        % create save file path
        filename_tif1 = strcat(filename, '_C1_render(', num2str(pixelSize), 'nm).tiff');
        filename_tif2 = strcat(filename, '_C2_render(', num2str(pixelSize), 'nm).tiff');
        saveAsTiff (x(ch1), y(ch1), xGrid, yGrid, fullfile(filepath, filename_tif1));
        saveAsTiff (x(ch2), y(ch2), xGrid, yGrid, fullfile(filepath, filename_tif2));
    elseif (numChannel == 1)
        filename_tif = strcat(filename, '_render(', num2str(pixelSize), 'nm).tiff');
        saveAsTiff (x, y, xGrid, yGrid, fullfile(filepath, filename_tif));
    else
        error(" channel number wrong! ");
    end
 
    


    function [cut1, cut2] = channelSeparation (dcr, autoSep, sigmaCutoff)
        cut1 = .35; cut2 = .4;
        if (autoSep)
            start_val = ones (numel(dcr), 1);
            start_val(dcr>0.4, :) = 2;
            model = fitgmdist(dcr, 2, 'Start', start_val);
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
        fprintf(" channel (dcr) cutoff at\n\t\t1: %d\n\t\t2: %d\n", cut1, cut2);
    end
    
    function saveAsTiff (x, y, xGrid, yGrid, pathSave)
        densityMap = histcounts2(x, y, xGrid, yGrid);
        densityMap = uint8(densityMap);
        t = Tiff(pathSave, 'w');
        % populate TIFF property structure, same for 2D/3D image
        tagstruct.ImageLength = size(densityMap, 1);
        tagstruct.ImageWidth = size(densityMap, 2);
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 8;
        tagstruct.SampleFormat = 1;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.Software = 'MATLAB';
        setTag(t, tagstruct);
        write(t, densityMap(:, :));
        writeDirectory(t);
        close(t);
        clear;
    end

