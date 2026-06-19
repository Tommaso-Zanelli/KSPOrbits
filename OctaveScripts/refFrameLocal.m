function [T] = refFrameLocal(x, v)
%
%   function [T] = refFrameLocal(x, v)
%
%   Local frame of reference (as in use in KSP)
%

    ir = x / norm(x, 2);
    i1 = v / norm(v, 2);
    i2 = cross(ir, i1);
    i2 = i2 / norm(i2, 2);
    i3 = cross(i1, i2);
    T = [i1, i2, i3];

end