function out = strnum(str, n)
out = strcat(str,'_',num2str((1:n)','%-d'));
out = strtrim(cellstr(out));
end