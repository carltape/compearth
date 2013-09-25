function q = convertq(q)
%CONVERTQ ensure that quaternion has the convention q = (w,x,y,z) with w > 0
%   See TapeTape2012 "Angle between principal axis triples"

[a,n] = size(q);
if a~=4, error('q must be 4 x n'); end
if ~isreal(q)
    %q
    disp('WARNING: q must be real-valued');
end

q = repmat(sign(q(1,:)),4,1) .* q;