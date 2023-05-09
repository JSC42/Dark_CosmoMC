function [r]=Find_Index(x_axis,x)
% Return index of x_axis element closest to x
% This script is 2000 times faster than built-in find function
% Warning: x must be asymptotically ascending and has no repeated elements
id1=1;
id3=length(x_axis);
xmin = x_axis(1);
xmax = x_axis(end);
nx = length(x_axis);

if x_axis(id1)>x_axis(id3)
    error('x must be asymptotically ascending!\n')
end
if xmin<0
    error('x_axis must all be positive\n')
end

if x<xmin
    if x>xmin*0.999
        r = 1;
    else
        error('overflow detected')
    end
elseif x>xmax
    if x<xmax*1.001
        r = id3;
    else
        error('overflow detected')
    end
else
    % Proceed if in range
    PROCEED=1;
    while PROCEED
        dx=floor((id3-id1)/2);
        id2=id1+dx;
        if x_axis(id1)<=x && x<x_axis(id2)
            id3=id2;
        else
            id1=id2;
        end
        
        if id3-id1<2
            PROCEED=0;
        end
    end
    r = id1;
    if r==nx
        r = nx-1;
    end
    
end
end