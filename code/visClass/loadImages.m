function [imData, txtData, imDims] = loadImages(directory, sz, preserveAR, checkForText, filterTextRegions)

txtData = [];

nX = sz(1);
nY = sz(2);
nC = sz(3);

filenames = getAllFiles(directory);

if (checkForText)
    textFilenames = getAllFiles(fullfile(directory,'\text\'));
    nText = length(textFilenames);
    haveTxt = zeros(nText,1);
    
    for i=1:nText
        haveTxt(i) = filenameToIndex(textFilenames{i});
    end
    filteredFilenames = [];
    for i=1:length(filenames)
        if (any(filenameToIndex(filenames{i}) == haveTxt))
            filteredFilenames = [filteredFilenames; filenames(i)];
        end
    end
    filenames = filteredFilenames;
    
    txtData = TextFeatureSet.empty(nText,0);
    for i = 1:length(filenames)
        txtData(i) = TextFeatureSet(textFilenames{i});
    end
end

imData = zeros(length(filenames),nX*nY*nC);
imDims = zeros(length(filenames),2);

for i = 1:length(filenames)
    I = imread(filenames{i});

    if (length(size(I)) == 3)
        I = rgb2gray(I);
    end
    imDims(i,:) = size(I);
    
    if (filterTextRegions)
        I = txtData(i).filterOutTextAreasSmooth(I);
    end
    if preserveAR
        I = imresize_constantAR(I,[nX nY]);
    else
        I = imresize(I, [nX nY]);
    end
    imData(i,:) = reshape(I,[1,nX*nY*nC]);
end

% meanDim = mean(imDims,1);
% disp(meanDim);

end

function fileList = getAllFiles(dirName)

dirData = dir(dirName);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
        fileList,'UniformOutput',false);
end

end
