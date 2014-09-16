function testX = loadImageData(data, Nimages, CIFAR_DIM)

corpusdir = data.basedir;
catdirs = data.catdirs;

fprintf('Loading data...\n');
tic;
Ncategories = length(catdirs);
testX = zeros(Nimages*Ncategories, prod(CIFAR_DIM));

for i=1:Ncategories
    catdir = catdirs{i};
    catImages = loadImages(strcat(corpusdir,catdir), CIFAR_DIM, false);
    rng = ((i-1)*Nimages + 1) : (i*Nimages);
    testX(rng,:) = catImages(1:Nimages,:);
end

testX = double(testX);
toc

end