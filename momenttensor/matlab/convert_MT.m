function [Mout,T] = convert_MT(i1,i2,M,boption)
%CONVERT_MT convert moment tensor matrices among different bases
%
% This program converts between different moment tensor conventions.
% All conventions are associated with a local coordinate system.
%
% M = [M11 M22 M33 M12 M13 M23]
%
% INPUT
%   i1      index of input moment tensor basis (see convert_getbasis.m)
%   i2      index of output moment tensor basis (see convert_getbasis.m)
%   M       6 x n set of moment tensors, M = [M11 M22 M33 M12 M13 M23]
%   boption choice of how to calculate the change of basis
%
% OUTPUT
%   Mout    6 x n set of moment tensors in basis of i2
%   T       transformation matrix to change basis of M from i1 to i2: Mout = T*M*T'
%
% See convert_getbasis.m for details
% Convention 1: up-south-east (GCMT)
% Convention 2: north-east-down (Aki and Richards, 1980)
% Convention 3: north-west-down (Stein and Wysession (2003)
% Convention 4: east-north-up
% Convention 5: south-east-down
%
% See also convertv.m, the vector version of this function.
%
% calls convert_getbasis.m
%
% Carl Tape, 11/2010
%

% two options to implement the transformation
%   true:   perform the linear algebra operation Mout = U*Min*U'
%   false:  swap entries and flip signs to achieve the correct transformation
% The 'true' option is mathematically cleanest but also has some
% possibility of introducing numerical errors.
if nargin < 4, boption = false; end

T = convert_getbasis(i1,i2);
if nargin==2
   disp('returning transformation matrix only, as requested');
   Mout = T;
   return
end

% make sure M is 6 x n
[M,n] = Mdim(M);

if i1==i2
    %error('i1 must differ from i2');
    disp('warning: i1 = i2, so no change');
    Mout = M;
    return
end

Mout = [];  % initialize

if boption   % transformation matrix T is from convert_getbasis.m
    % convert from 6 x n to 3 x 3 x n
    Mmat = Mvec2Mmat(M,1);
    % loop over each matrix and apply Mout = T*Min*T'
    Moutmat = NaN(size(Mmat));
    for ii=1:n
        Moutmat(:,:,ii) = T*Mmat(:,:,ii)*T';
    end
    % convert back to 6 x n (note: this assumes that Mout is symmetric)
    Mout = Mvec2Mmat(Moutmat,0);
    
else            % WARNING: THESE FORMULAS ARE TIED TO SPECIFIC BASES
    if i1==1
        if i2==2        % up-south-east (GCMT) to north-east-down (AkiRichards) (AR, 1980, p. 118)
            Mout(1,:) = M(2,:);
            Mout(2,:) = M(3,:);
            Mout(3,:) = M(1,:);
            Mout(4,:) = -M(6,:);
            Mout(5,:) = M(4,:);
            Mout(6,:) = -M(5,:);
        elseif i2==3    % up-south-east (GCMT) to north-west-up (/opt/seismo-util/bin/faultpar2cmtsol.pl)
            Mout(1,:) = M(2,:);
            Mout(2,:) = M(3,:);
            Mout(3,:) = M(1,:);
            Mout(4,:) = M(6,:);
            Mout(5,:) = -M(4,:);
            Mout(6,:) = -M(5,:);
        elseif i2==4    % up-south-east (GCMT) to east-north-up
            Mout(1,:) = M(3,:);
            Mout(2,:) = M(2,:);
            Mout(3,:) = M(1,:);
            Mout(4,:) = -M(6,:);
            Mout(5,:) = M(5,:);
            Mout(6,:) = -M(4,:);
        elseif i2==5    % up-south-east (GCMT) to south-east-up
            Mout(1,:) = M(2,:);
            Mout(2,:) = M(3,:);
            Mout(3,:) = M(1,:);
            Mout(4,:) = M(6,:);
            Mout(5,:) = M(4,:);
            Mout(6,:) = M(5,:);  
        end

    elseif i1==2
        if i2==1        % north-east-down (AkiRichards) to up-south-east (GCMT) (AR, 1980, p. 118)
            Mout(1,:) = M(3,:);
            Mout(2,:) = M(1,:);
            Mout(3,:) = M(2,:);
            Mout(4,:) = M(5,:);
            Mout(5,:) = -M(6,:);
            Mout(6,:) = -M(4,:);
        elseif i2==3    % north-east-down (AkiRichards) to north-west-up
            Mout(1,:) = M(1,:);
            Mout(2,:) = M(2,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = -M(4,:);
            Mout(5,:) = -M(5,:);
            Mout(6,:) = M(6,:);   
        elseif i2==4    % north-east-down (AkiRichards) to east-north-up
            Mout(1,:) = M(2,:);
            Mout(2,:) = M(1,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = M(4,:);
            Mout(5,:) = -M(6,:);
            Mout(6,:) = -M(5,:);
        elseif i2==5    % north-east-down (AkiRichards) to south-east-up
            Mout(1,:) = M(1,:);
            Mout(2,:) = M(2,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = -M(4,:);
            Mout(5,:) = M(5,:);
            Mout(6,:) = -M(6,:);   
        end

    elseif i1==3
        if i2==1        % north-west-up to up-south-east (GCMT)
            Mout(1,:) = M(3,:);
            Mout(2,:) = M(1,:);
            Mout(3,:) = M(2,:);
            Mout(4,:) = -M(5,:);
            Mout(5,:) = -M(6,:);
            Mout(6,:) = M(4,:);
        elseif i2==2    % north-west-up to north-east-down (AkiRichards)
            Mout(1,:) = M(1,:);
            Mout(2,:) = M(2,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = -M(4,:);
            Mout(5,:) = -M(5,:);
            Mout(6,:) = M(6,:); 
        elseif i2==4    % north-west-up to east-north-up
            Mout(1,:) = M(2,:);
            Mout(2,:) = M(1,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = -M(4,:);
            Mout(5,:) = -M(6,:);
            Mout(6,:) = M(5,:); 
        elseif i2==5    % north-west-up to south-east-up
            Mout(1,:) = M(1,:);
            Mout(2,:) = M(2,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = M(4,:);
            Mout(5,:) = -M(5,:);
            Mout(6,:) = -M(6,:); 
        end

    elseif i1==4
        if i2==1        % east-north-up to up-south-east (GCMT)
            Mout(1,:) = M(3,:);
            Mout(2,:) = M(2,:);
            Mout(3,:) = M(1,:);
            Mout(4,:) = -M(6,:);
            Mout(5,:) = M(5,:);
            Mout(6,:) = -M(4,:);
        elseif i2==2    % east-north-up to north-east-down (AkiRichards)
            Mout(1,:) = M(2,:);
            Mout(2,:) = M(1,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = M(4,:);
            Mout(5,:) = -M(6,:);
            Mout(6,:) = -M(5,:);
        elseif i2==3    % east-north-up to north-west-up
            Mout(1,:) = M(2,:);
            Mout(2,:) = M(1,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = -M(4,:);
            Mout(5,:) = M(6,:);
            Mout(6,:) = -M(5,:); 
        elseif i2==5    % east-north-up to south-east-up
            Mout(1,:) = M(2,:);
            Mout(2,:) = M(1,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = -M(4,:);
            Mout(5,:) = -M(6,:);
            Mout(6,:) = M(5,:); 
        end

    elseif i1==5        % south-east-up to up-south-east (GCMT)
        if i2==1
            Mout(1,:) = M(3,:);
            Mout(2,:) = M(1,:);
            Mout(3,:) = M(2,:);
            Mout(4,:) = M(5,:);
            Mout(5,:) = M(6,:);
            Mout(6,:) = M(4,:);
        elseif i2==2    % south-east-up to north-east-down (AkiRichards)
            Mout(1,:) = M(1,:);
            Mout(2,:) = M(2,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = -M(4,:);
            Mout(5,:) = M(5,:);
            Mout(6,:) = -M(6,:);
        elseif i2==3    % south-east-up to north-west-up
            Mout(1,:) = M(1,:);
            Mout(2,:) = M(2,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = M(4,:);
            Mout(5,:) = -M(5,:);
            Mout(6,:) = -M(6,:);
        elseif i2==4    % south-east-up to east-north-up
            Mout(1,:) = M(2,:);
            Mout(2,:) = M(1,:);
            Mout(3,:) = M(3,:);
            Mout(4,:) = -M(4,:);
            Mout(5,:) = M(6,:);
            Mout(6,:) = -M(5,:); 
        end
    end
end  % boption
    
%==========================================================================
% EXAMPLES

if 0==1
    % transformation matrix only
    i1 = 1; i2 = 2;
    T = convert_MT(i1,i2)
    
    % simple example
    i1 = 1; i2 = 2;
    A = [1:6]'
    M1 = convert_MT(i1,i2,A)    % convert from i1 to i2
    M2 = convert_MT(i2,i1,M1)   % convert back from i2 to i1
    
    % checking the transformation matrix
    i1 = 1; i2 = 5;
    M1 = rand(6,1);
    [M2,T] = convert_MT(i1,i2,M1);
    Mvec2Mmat(M1,1)                 % up-south-east
    Mcheck = T*Mvec2Mmat(M1,1)*T'   % south-east-up
    Mvec2Mmat(M2,1)                 % (check)
    % example vector v
    v1 = rand(3,1)                  % up-south-east
    v2 = T*v1                       % south-east-up
    % example symmetric matrix X1
    X = randi(10,3); X1=X'*X, X2=T*X1*T'
    
    % check the two different implementations (boption = true/false) for
    % all possible change of bases
    M1 = rand(6,1);     % single symmetric matrix with random entries
    NTYPE = 5;
    for i1=1:NTYPE
        for i2=2:NTYPE
            if i2==i1, continue; end
            M2true  = convert_MT(i1,i2,M1,true);
            M2false = convert_MT(i1,i2,M1,false);
            [M2true M2false M2true-M2false]
            if norm(M2true-M2false) > 1e-6
                error('from %i to %i',i1,i2);
            end
        end
    end
    
%     % check all possible transformations
%     M1 = rand(6,1);     % single symmetric matrix with random entries
%     NTYPE = 5;
%     for i1=1:NTYPE
%         for i2=2:NTYPE
%             if i2==i1, continue; end
%             [M2,T] = convert_MT(i1,i2,M1);
%             if and(~isempty(T),~isempty(M2))
%                 Mcheck = T*Mvec2Mmat(M1,1)*T';
%                 ncheck = norm(Mvec2Mmat(M2,1)-Mcheck);
%                 % display info if the numerical check fails
%                 if ncheck > 1e-6
%                    disp(sprintf('from %i to %i',i1,i2));
%                    Mvec2Mmat(M2,1),Mcheck
%                    error('check') 
%                 end
%             end
%         end
    end
    
end

%==========================================================================

