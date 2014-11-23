function ret = adjacentsegs(seg1, seg2)
    ret = (seg1(2) <= (seg2(1)-1)) || (seg1(1) >= (seg2(2)+1));
    ret = ~ret;
end