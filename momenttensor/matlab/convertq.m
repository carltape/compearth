function q = convertq(q)
%CONVERTQ ensure that quaternion has the convention q = (w,x,y,z) with w >= 0
%   See TapeTape2012 "Angle between principal axis triples"

[a,n] = size(q);
if a~=4, error('q must be 4 x n'); end
if ~isreal(q)
    %q
    disp('WARNING: q must be real-valued');
end

% if w=0, do not multiply by 0; otherwise you zero-out the rotation axis (x,y,z)
qfac = sign(q(1,:));
qfac(qfac==0) = 1;
q = repmat(qfac,4,1) .* q;