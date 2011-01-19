%
% function atanx = arctan(x)
% CARL TAPE, 27-Oct-2005
% printed xxx
%
% Modified arc-tangent function.
% Output is in the range [0, 180], not [-90, 90].
%
% calls xxx
% called by xxx
%

function atanx = arctan(x)

atanx = atan(x);

ineg = find(x < 0);
atanx(ineg) = atanx(ineg) + pi;

% n = length(x);
% atanx = zeros(1,n);
% 
% for i=1:n
%     if x(i) <= 0,
%         atanx(i) = pi + atan(x(i));
%     else
%         atanx(i) = atan(x(i));
%     end
% end

%========================================
