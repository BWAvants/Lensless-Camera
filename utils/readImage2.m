function outImage = readImage2(dataPath,dataName,cameraSize,numberOfFrames)

nameImg = [dataPath dataName '_1.bmp'];
cameraSize = size((imread(nameImg)));
outImage = zeros(cameraSize);
for i = 1:numberOfFrames
    outImage = outImage + double(imread([dataPath dataName '_' num2str(i) '.bmp']));
end
outImage = mean(outImage/numberOfFrames,3);
