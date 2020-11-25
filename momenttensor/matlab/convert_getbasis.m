function T = convert_getbasis(i1,i2)
%CONVERT_GETBASIS get transformation matrix to go between moment tensor bases
%
% INPUT
%   i1  index of input moment tensor basis
%   i2  index of output moment tensor basis
%
% OUTPUT
%   T   transformation matrix to change basis of M from i1 to i2: Mout = T*M*T'
%
% Convention 1: up-south-east (GCMT) (www.globalcmt.org)
%   1: up (r), 2: south (theta), 3: east (phi)
%
% Convention 2: Aki and Richards (1980, p. 114-115, 118)
%   also Jost and Herrman (1989, Fig. 1)
%   1: north, 2: east, 3: down
%
% Convention 3: Stein and Wysession (2003, p. 218)
%   also TapeTape2012a "A geometric setting for moment tensors" (p. 478)
%   also several Kanamori codes
%   1: north, 2: west, 3: up
% 
% Convention 4: 
%   1: east, 2: north, 3: up
% 
% Convention 5: TapeTape2013 "The classical model for moment tensors" (p. 1704)
%   1: south, 2: east, 3: up
%
% EXAMPLE: T = convert_getbasis(2,1)
%
% called by convert_MT.m, convertv.m 
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

stlabs = {'up-south-east (GCMT)',...
    'north-east-down (AkiRichards)',...
    'north-west-up',...
    'east-north-up',...
    'south-east-up'};
disp(sprintf('convert_getbasis.m: %s to %s',stlabs{i1},stlabs{i2}));

%==========================================================================

