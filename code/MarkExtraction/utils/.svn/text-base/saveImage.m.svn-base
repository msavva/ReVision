function saveImage(imHandle, filename, suffix, outpath)
    extInd = strfind(filename, '.');
    outFile = strcat(filename(1:extInd-1),suffix);
    saveas(imHandle, fullfile(outpath,outFile));
end