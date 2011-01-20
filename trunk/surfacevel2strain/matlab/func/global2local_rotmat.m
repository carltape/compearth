%
% function T = global2local_rotmat(th,ph)
% Carl Tape, 03-Sept-2005
%
% This function converts from a global basis (x,y,z) to a local basis in
% (r, th, ph) or, equivalently, (up, south, east).  In other words, the
% components of the vector R will depend on the choice of basis: this
% program generates a matrix that converts between two possible bases.
%
% EXAMPLES in test_global2local.m.
%
% After Cox and Hart (1986), p. 155.
% Convention is from Dahlen and Tromp (1998), p. 832-833.
%
% format for conversion matrix (note ordering: r --> th --> ph):
%
%     |  r . x    r . y    r . z |   |  sin(th)cos(ph)  sin(th)sin(ph)   cos(th)  |
%     |                          |   |                                            |
% T = | th . x   th . y   th . z |   |  cos(th)cos(ph)  cos(th)sin(ph)   -sin(th) |
%     |                          |   |                                            |
%     | ph . x   ph . y   ph . z | = |  -sin(ph)        cos(ph)          0        |
%
% calls xxx 
% called by global2local.m, test_global2local.m
%

function T = global2local_rotmat(th,ph)

% convention: r, th, ph
T = zeros(3,3);
T(1,1) = sin(th)*cos(ph);
T(1,2) = sin(th)*sin(ph);
T(1,3) = cos(th);
T(2,1) = cos(th)*cos(ph);
T(2,2) = cos(th)*sin(ph);
T(2,3) = -sin(th);
T(3,1) = -sin(ph);
T(3,2) = cos(ph);
T(3,3) = 0;

% convention: th, ph, r
% T(1,1) = cos(th)*cos(ph);
% T(1,2) = cos(th)*sin(ph);
% T(1,3) = -sin(th);
% T(2,1) = -sin(ph);
% T(2,2) = cos(ph);
% T(2,3) = 0;
% T(3,1) = sin(th)*cos(ph);
% T(3,2) = sin(th)*sin(ph);
% T(3,3) = cos(th);

%==============================================================
