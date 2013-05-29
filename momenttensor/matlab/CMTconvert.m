function Mout = CMTconvert(M,i1,i2)
%CMTCONVERT convert moment tensor matrices among different bases
%
% This program converts between different moment tensor conventions.
% All conventions are associated with a local coordinate system.
%
% M = [M11 M22 M33 M12 M13 M23]
%
% Convention 1: Harvard CMT
%   1: up (r), 2: south (theta), 3: east (phi)
%
% Convention 2: Aki and Richards
%   1: north, 2: east, 3: down
%
% Convention 3: Kanamori
%   1: north, 2: west, 3: up
% 
% Convention 4:
%   1: east, 2: north, 3: up
% 
% Convention 5:
%   1: south, 2: east, 3: up
%
% NOTE:
% The transformations are best thought of a change-of-basis operations such
% as M = T * Mp * T', where T is the transformation matrix, for example
% T = [0 0 1 ; 1 0 0 ; 0 1 0] to go from Harvard to C5.
% 
% Carl Tape, 11/2010
%

if i1==i2, error('i1 must differ from i2'); end

% make sure M is 6 x n
[M,n] = Mdim(M);

Mout = zeros(6,n);

stlabs = {'GCMT (up-south-east)',...
    'Aki (north-east-down)',...
    'Kanamori (north-west-up)',...
    'C4 (east-north-up)',...
    'C5 (south-east-up)'};
disp(sprintf('CMTconvert.n: %s to %s',stlabs{i1},stlabs{i2}));

if i1==1
    if i2==2        % Harvard to Aki (AR, 1980, p. 118)
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(3,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = -M(6,:);
        Mout(5,:) = M(4,:);
        Mout(6,:) = -M(5,:);
    elseif i2==3    % Harvard to Kanamori (/opt/seismo-util/bin/faultpar2cmtsol.pl)
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(3,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = M(6,:);
        Mout(5,:) = -M(4,:);
        Mout(6,:) = -M(5,:);
    elseif i2==4    % Harvard to C4
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = -M(6,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(4,:);
    elseif i2==5    % Harvard to C5
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(3,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = M(6,:);
        Mout(5,:) = M(4,:);
        Mout(6,:) = M(5,:);  
    end
    
elseif i1==2
    if i2==1        % Aki to Harvard (AR, 1980, p. 118)
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(2,:);
        Mout(4,:) = M(5,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = -M(4,:);
    elseif i2==3    % Aki to Kanamori
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(6,:);   
    elseif i2==4
        error('not yet implemented to convert from Aki to C4');  
    elseif i2==5
        error('not yet implemented to convert from Aki to C5');     
    end
    
elseif i1==3
    if i2==1        % Kanamori to Harvard
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(2,:);
        Mout(4,:) = -M(5,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = M(4,:);
    elseif i2==2    % Kanamori to Aki
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(6,:); 
    elseif i2==4
        error('not yet implemented to convert from Kanamori to C4'); 
    elseif i2==5
        error('not yet implemented to convert from Kanamori to C5');     
    end
    
elseif i1==4
    if i2==1        % C4 to Harvard
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = -M(6,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(4,:);
    elseif i2==2
        error('not yet implemented to convert from C4 to Aki');
    elseif i2==3
        error('not yet implemented to convert from C4 to Kanamori');
    elseif i2==5
        error('not yet implemented to convert from C4 to C5');    
    end
    
elseif i1==5        % C5 to Harvard
    if i2==1
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(2,:);
        Mout(4,:) = M(5,:);
        Mout(5,:) = M(6,:);
        Mout(6,:) = M(4,:);
    elseif i2==2
        error('not yet implemented to convert from C5 to Aki');
    elseif i2==3
        error('not yet implemented to convert from C5 to Kanamori');
    elseif i2==4
        error('not yet implemented to convert from C5 to C4');    
    end
    
end

%==========================================================================
% EXAMPLES

if 0==1
    A = [1:6]'
    M1 = CMTconvert(A,1,5)
    M2 = CMTconvert(M1,5,1)
end

%==========================================================================

