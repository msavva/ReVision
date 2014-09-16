function n = filenameToIndex(fname)
[~, base, ~] = fileparts(fname);
[~, number] = strtok(base,'_');
n = str2double(number(2:end));
end