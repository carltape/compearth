function [Mout,T] = CMTconvert(M,i1,i2)
%CMTCONVERT convert moment tensor matrices among different bases
%
% This program converts between different moment tensor conventions.
% All conventions are associated with a local coordinate system.
%
% M = [M11 M22 M33 M12 M13 M23]
%
% INPUT
%   M       6 x n set of moment tensors, M = [M11 M22 M33 M12 M13 M23]
%   i1      index of input moment tensor basis
%   i2      index of output moment tensor basis
%
% OUTPUT
%   Mout    6 x n set of moment tensors in basis of i2
%   P       transformation matrix to change basis of M from i1 to i2: Mout = T*M*T'
%
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
% Convention 5 (Tape and Tape, 2012):
%   1: south, 2: east, 3: up
%
% Carl Tape, 11/2010
%

NTYPE = 5;  % number of right-handed bases to consider

% permutation matrices
Tall = cell(NTYPE,NTYPE);
Tall{1,1} = eye(3);
Tall{2,2} = eye(3);
Tall{3,3} = eye(3);
Tall{4,4} = eye(3);
Tall{5,5} = eye(3);
% from i1 to i2
Tall{1,2} = [0 -1 0 ; 0 0 1 ; -1 0 0];
Tall{1,3} = [0 -1 0 ; 0 0 -1 ; 1 0 0];
Tall{1,4} = [0 0 1 ; 0 -1 0 ; 1 0 0];
Tall{1,5} = [0 1 0 ; 0 0 1 ; 1 0 0];
Tall{2,3} = [1 0 0 ; 0 -1 0 ; 0 0 -1];
Tall{2,4} = [0 1 0 ; 1 0 0 ; 0 0 -1];
Tall{2,5} = [-1 0 0 ; 0 1 0 ; 0 0 -1];
Tall{3,4} = [0 -1 0 ; 1 0 0 ; 0 0 1];
Tall{3,5} = [-1 0 0 ; 0 -1 0 ; 0 0 1];
Tall{4,5} = [0 -1 0 ; 1 0 0 ; 0 0 1];
% from i2 to i1
Tall{2,1} = Tall{1,2}';
Tall{3,1} = Tall{1,3}';
Tall{3,2} = Tall{2,3}';
Tall{4,1} = Tall{1,4}';
Tall{4,2} = Tall{2,4}';
Tall{4,3} = Tall{3,4}';
Tall{5,1} = Tall{1,5}';
Tall{5,2} = Tall{2,5}';
Tall{5,3} = Tall{3,5}';
Tall{5,4} = Tall{4,5}';
% transformation matrix
T = Tall{i1,i2};

if i1==i2
    %error('i1 must differ from i2');
    disp('warning: i1 = i2, so no change');
    Mout = M;
    return
end

% make sure M is 6 x n
[M,n] = Mdim(M);

Mout = zeros(6,n);

stlabs = {'GCMT (up-south-east)',...
    'Aki (north-east-down)',...
    'Kanamori (north-west-up)',...
    'C4 (east-north-up)',...
    'C5 (south-east-up)'};
disp(sprintf('CMTconvert.n: %s to %s',stlabs{i1},stlabs{i2}));

Mout = [];

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
        Mout(5,:) = -M(5,:);
        Mout(6,:) = M(6,:);   
    elseif i2==4    % Aki to C4
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = M(4,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = -M(5,:);
    elseif i2==5    % Aki to C5
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(6,:);   
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
        Mout(5,:) = -M(5,:);
        Mout(6,:) = M(6,:); 
    elseif i2==4    % Kanamori to C4
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = M(5,:); 
    elseif i2==5    % Kanamori to C5
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = M(4,:);
        Mout(5,:) = -M(5,:);
        Mout(6,:) = -M(6,:); 
    end
    
elseif i1==4
    if i2==1        % C4 to Harvard
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = -M(6,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(4,:);
    elseif i2==2    % C4 to Aki
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = M(4,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = -M(5,:);
    elseif i2==3    % C4 to Kanamori
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(6,:);
        Mout(6,:) = -M(5,:); 
    elseif i2==5    % C4 to C5
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = M(5,:); 
    end
    
elseif i1==5        % C5 to Harvard
    if i2==1
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(2,:);
        Mout(4,:) = M(5,:);
        Mout(5,:) = M(6,:);
        Mout(6,:) = M(4,:);
    elseif i2==2    % C5 to Aki
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(6,:);
    elseif i2==3    % C5 to Kanamori
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = M(4,:);
        Mout(5,:) = -M(5,:);
        Mout(6,:) = -M(6,:);
    elseif i2==4    % C5 to C4
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(6,:);
        Mout(6,:) = -M(5,:); 
    end
    
end

%==========================================================================
% EXAMPLES

if 0==1
    A = [1:6]'
    M1 = CMTconvert(A,1,5)
    M2 = CMTconvert(M1,5,1)
    
    %% checking the transformation matrix
    i1 = 1; i2 = 5;
    M1 = rand(6,1);
    [M2,T] = CMTconvert(M1,i1,i2);
    Mvec2Mmat(M1,1)                 % up-south-east
    Mcheck = T*Mvec2Mmat(M1,1)*T'   % south-east-up
    Mvec2Mmat(M2,1)                 % (check)
    v1 = rand(3,1)                  % up-south-east
    v2 = T*v1                       % south-east-up
    X1 = randi(10,3); X1=X1'*X1,  X2=T*X*T'
    
    %% check all possible transofrmations
    M1 = rand(6,1);
    NTYPE = 5;
    for i1=1:NTYPE
        for i2=2:NTYPE
            if i2==i1, continue; end
            [M2,T] = CMTconvert(M1,i1,i2);
            if and(~isempty(T),~isempty(M2))
                Mcheck = T*Mvec2Mmat(M1,1)*T';
                % display info if the check fails
                if ncheck > 1e-6
                   disp(sprintf('from %i to %i',i1,i2));
                   Mvec2Mmat(M2,1),Mcheck
                   v1 = rand(3,1) 
                   v2 = T*v1
                   error('check') 
                end
            end
        end
    end
end

%==========================================================================
