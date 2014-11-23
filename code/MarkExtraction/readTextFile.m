function [imsize textlocs] = readTextFile(file)
    tf = fopen(file);
    
    % Get rid of the image size header
    fgetl(tf);
    imsize = fscanf(tf, '%d,%d\n');
    
    % Get second header (text locations)
    headerline = fgetl(tf);
    headers = textscan(headerline,'%q','Delimiter',',');
    data = textscan(tf, '%f%f%f%f%f%q%q','Delimiter',',');
    
    textlocs = struct();
    for k = 1:length(headers{:})
        eval(['textlocs.' headers{1}{k} '= data{k};']);
    end
    
    fclose(tf);
end