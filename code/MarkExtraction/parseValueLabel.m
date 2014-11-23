function lab = parseValueLabel(curtext)
    lab = {};
    lab.text = curtext;
    
    % Remove some common characters from the label

    % Remove percentages
    if ~isempty(findstr(curtext, '%'))
        curtext = curtext(1:findstr(curtext,'%')-1);
    end
    % Remove $ sign
    if ~isempty(findstr(curtext, '$'))
        curtext = curtext(findstr(curtext,'$')+1:end);
    end
    % Remove commas
    if ~isempty(findstr(curtext, ','))
        commas = findstr(curtext,',');
        nocommas = '';
        for cc=1:length(commas)
            if cc == 1
                nocommas = strcat(nocommas, curtext(1:commas(cc)-1));
            elseif cc == length(commas)
                nocommas = strcat(nocommas, curtext(commas(cc)+1:end));
            else
                nocommas = strcat(nocommas, curtext(commas(cc)+1:commas(cc+1)-1));
            end
        end
        curtext = nocommas;
    end

    lab.value = str2num(curtext);
end