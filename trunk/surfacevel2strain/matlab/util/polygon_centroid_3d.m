function centroid = polygon_centroid_3d( n, v )

%% POLYGON_CENTROID_3D computes the centroid of a polygon in 3D.
%
%  Discussion:
%
%    The polygon is described by its vertices.  In many applications,
%    these vertices will lie in a common plane, and the polygon will
%    be "flat".  However, that is not required for this formula.
%
%    This formula triangulates the polygon, computes the area of
%    each triangle and its centroid, and then computes the centroid
%    of the polygon as the weight-averaged sum of the triangle centroids.
%
%    Thanks to Sam Hemmings for pointing a missing transpose sign,
%    05 April 2005.
%
%  Modified:
%
%    05 April 2005
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Adrian Bowyer and John Woodwark,
%    A Programmer's Geometry,
%    Butterworths, 1983.
%
%  Parameters:
%
%    Input, integer N, the number of vertices of the polygon.
%
%    Input, real V(3,N), the coordinates of the vertices.
%
%    Output, real CENTROID(3), the coordinates of the centroid.
%
  dim_num = 3;
  
  area = 0.0;
  centroid(1:dim_num,1) = 0.0;

  for i = 1 : n - 2

    t(1:dim_num,1:3) = [ v(1:dim_num,i)'; v(1:dim_num,i+1)'; v(1:dim_num,n)' ]';

    area_triangle = triangle_area_3d ( t );

    area = area + area_triangle;

    centroid(1:dim_num,1) = centroid(1:dim_num,1) + area_triangle ...
      * ( v(1:dim_num,i) + v(1:dim_num,i+1) + v(1:dim_num,n) ) / 3.0;

  end

  if ( area == 0.0 )
    centroid(1:dim_num,1) = v(1:dim_num,1);
  else
    centroid(1:dim_num,1) = centroid(1:dim_num) / area;
  end
