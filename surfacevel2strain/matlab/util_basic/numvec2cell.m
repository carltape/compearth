function z = numvec2cell(numvec)
%NUMVEC2CELL convert a vector of integers into a vector of strings
%
% This is useful for plotting text labels on figures.
%
% EXAMPLE: z = numvec2cell([1:12])
%

z = strtrim(cellstr(num2str(numvec(:))));
