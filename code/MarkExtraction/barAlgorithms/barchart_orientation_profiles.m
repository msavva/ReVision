% Return the orientation of a bar chart based on its projection profiles.
%
% INPUT
%
% profiles: structure containing x and y projection profiles
%
% OUTPUT
%
% o: one of 'h' or 'v'
function o = barchart_orientation_profiles(profiles)
    % Find local maxima that are greater than one stdev away from the average
    % in the projection profiles
    [xPeaks, xLocs] = findpeaks(profiles.x, 'MINPEAKHEIGHT', ...
        sum(profiles.x)/length(profiles.x) + std(profiles.x));
    [yPeaks, yLocs] = findpeaks(profiles.y, 'MINPEAKHEIGHT', ...
        sum(profiles.y)/length(profiles.y) + std(profiles.y));

    if mean(xPeaks) > mean(yPeaks)
        o = 'h';
    else
        o = 'v';
    end
end