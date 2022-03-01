diffuserDirectory = 'C:\Users\admin\OneDrive - Faurecia\Desktop\dataset - v9\diffuser\local dimming\lpf';
origDirectory = 'C:\Users\admin\OneDrive - Faurecia\Desktop\dataset - v9\input image\imageFit';

% filePattern = fullfile(origDirectory, '*.JPG');
filePattern = fullfile(origDirectory, '*.png');

imageNames = dir(filePattern);
 
for index = 1 : 1 % number of images
    imageName = getfield(imageNames(index),'name')
%     newStr = erase(imageName,".JPG")
    newStr = erase(imageName,".png")

%     newStr = erase(newStr,".jpg")
    newStr = erase(newStr,".png")


    imageDirectoryDiff = strcat( diffuserDirectory, '\', newStr, ' SimulationBacklightDiffusionModule_local dimming.png');
%     imageName = strcat('\', newStr, '.jpg');
    imageName = strcat('\', newStr, '.png');

    imageDirectoryOrig = strcat( origDirectory, imageName);

    IDiff = imread(imageDirectoryDiff);
    IOrig = imread(imageDirectoryOrig);
    IOrig = imresize(IOrig,[720 1920]);

    % Read json TC
    T = loadjson( join(['YinAll2.json'], '') );

    bit = 255;
    level = 10;
    thresh = multithresh(IDiff, level);
    seg_I = imquantize(IDiff, thresh);
    RGB = label2rgb(seg_I); 	 
    figure(1);
    imshow(RGB);
    axis off;
    title('RGB Segmented Image');
    file_name = strcat('C:\Users\admin\OneDrive - Faurecia\Desktop\dataset - v9\threshold level\LPF\', newStr, '.png');
    imwrite(RGB, file_name);
    % NewImage = I .* (I >= 0 & I <= 11);
    thresh = [0, thresh, 255];
    images = cell(length(thresh)-1,1);
    compensated = cell(length(thresh)-1,1);
    for index = 1 : length(thresh)-1
        newImage = IDiff >= thresh(index) & IDiff < thresh(index + 1);
        mask = uint8(cat(3, newImage, newImage, newImage));
        images{index} = mask .* IOrig;
    end

    for index = 1 : length(thresh) - 1
        %figure(index + 1);
        %imshow(images{index});
        %axis off;
        %title('threshold Segmented Image');
        (double(thresh(index)) ./ 255) .* 1000
        compensated{index} = lrtLocalDimming(double(images{index}) ./255, 1000, (double(thresh(index+1)) ./ 255).*1000, (double(thresh(index+1)) ./ 255).*1000, T);

    end

    imageWholeComp = zeros(size(IOrig),'like', IOrig);
    for index = 1 : length(thresh) - 1
     %   figure(index+27);
        imageInt = uint8(compensated{index} .* 255);
      %  imshow(imageInt);
        imageWholeComp = imageWholeComp + imageInt;
       % axis off;
       % title('threshold Segmented Image');
    end

    %figure (299)
    %imshow(imageWholeComp)
    file_name = strcat('C:\Users\admin\OneDrive - Faurecia\Desktop\dataset - v9\input image\local dimming\pdp_source500\', newStr, '.png');
    imwrite(imageWholeComp, file_name);
end