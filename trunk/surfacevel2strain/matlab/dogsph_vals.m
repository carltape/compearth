%
% function ff = dogsph_vals(clon,lat,q,lon_vec,lat_vec,opts)
% Pablo Muse, April 3, 2007
%
% Given a gridpoint (clon,clat), return the value of the DOG wavelet 
% of scale 1/2^(q-2) centered at this point, and its first derivative 
% evaluated at all the input datapoints.
%
% INPUT:
%   clon, clat, q   = these describe the DOG wavelet. q is the folding
%                     order of the triangular mesh (must be q>2)
%   opts            = how many columns of ff you want returned (derivatives)
%   lon_vec,lat_vec = datapoints at which you want the wavelet evaluated
%
%

function ff = dogsph_vals(clon, clat, q, lon_vec, lat_vec, opts)

% convert to theta-phi
deg    = 180/pi;
ph     = clon/deg;
th     = (90-clat)/deg;
ph_vec = lon_vec/deg;
th_vec = (90-lat_vec)/deg;

% options and parameters -- q controls the scale (width) of the spline
ncol  = opts{1};      % number of columns of ff to return
ndata = length(lon_vec);

alpha = 1.25;
a0 = 1;
aj = a0./ 2.^q;
a_alpha = aj*alpha;
 
% data points
cos_th_vec = cos(th_vec);
sin_th_vec = sin(th_vec);
cos_ph_vec = cos(ph_vec);
sin_ph_vec = sin(ph_vec);
X = sin_th_vec.*cos_ph_vec;
Y = sin_th_vec.*sin_ph_vec;
Z = cos_th_vec;
% wavelet center
cX = sin(th).*cos(ph);
cY = sin(th).*sin(ph);
cZ = cos(th);

%norm_cst_j = 1/feval('dogsph',0,0,0,0,0,0,aj,0,waveopts{:}); % l_infty normalization
norm_cst_j = 1; % l2 normalization

if (ncol ~= 1 && ncol ~= 3 &&  ncol ~= 6) 
    display('Something may be wrong: the maximum number of columns in dogsph_vals is 3, and was exceeded');
end


% columns of ff :
% (1) f, function value
% (2) df/dph
% (3) df/dth
% (4) surf_del2 -- depends only on delta
% (5) |del f|   -- depends only on delta

ff = zeros(ndata,ncol);
if (ncol >= 1)
    sqr_dist = (X - cX).*(X - cX) + (Y - cY).*(Y - cY) + (Z - cZ).*(Z - cZ);
    tan2_halfangle = sqr_dist./(4 - sqr_dist);
    sqrt_lambda_a = 2*aj./(2*aj.*aj + 0.5*(1-aj.*aj).*sqr_dist);
    sqrt_lambda_a_alpha = 2*a_alpha./(2*a_alpha.*a_alpha + 0.5*(1-a_alpha.*a_alpha).*sqr_dist);
    ff(:,1) = sqrt_lambda_a.*exp(-(1./aj.^2).*tan2_halfangle) - ...
        (1/alpha)*sqrt_lambda_a_alpha.*exp(-(1./a_alpha.^2).*tan2_halfangle);
    % function value  
    %ff(:,1) = feval('dogsph',X,Y,Z,cX,cY,cZ,aj,0,waveopts{:})*norm_cst_j;
end

if (ncol >= 2)
   % derivatives (to do: a mexfile)
    derivative_tan2_halfangle = 4./((4-sqr_dist).^2);
    tmp = 2*aj.*aj + 0.5*(1-aj.*aj).*sqr_dist;
    derivative_sqrt_lambda_a = -aj*(1-aj.*aj)./tmp.^2;
    tmp_alpha = 2*a_alpha.*a_alpha + 0.5*(1-a_alpha.*a_alpha).*sqr_dist;
    derivative_sqrt_lambda_a_alpha = -a_alpha.*(1-a_alpha.*a_alpha)./...
        tmp_alpha.^2;

    % derivative of the pair of dilated gaussian wrt sqr_dist
    derivative_ga = (derivative_sqrt_lambda_a - (1./aj.^2).*sqrt_lambda_a.*derivative_tan2_halfangle)...
        .*exp(-(1./aj.^2).*tan2_halfangle);
    derivative_ga_alpha = (derivative_sqrt_lambda_a_alpha - (1./a_alpha.^2).*sqrt_lambda_a_alpha.*derivative_tan2_halfangle)...
        .*exp(-(1./a_alpha.^2).*tan2_halfangle);

    % derivative of the whole DOG wavelet wrt sqr_dist (the distance of a
    % point to the center of the wavelet)
    derivative_gaussians = norm_cst_j*(derivative_ga - (1/alpha)*derivative_ga_alpha);

    % derivative of sqr_dist wrt colatitude theta
    derivative_sqr_dist_dth = -2*(cX*cos_ph_vec + cY*sin_ph_vec).*cos_th_vec + 2*cZ*sin_th_vec;
       
    % derivatives of  sqr_dist wrt longitude phi
    derivative_sqr_dist_dph = 2*(cX*sin_ph_vec - cY*cos_ph_vec).*sin_th_vec;
    
    % final result by chain rule
    ff(:,2) = derivative_gaussians.*derivative_sqr_dist_dph; % df/dph
    ff(:,3) = derivative_gaussians.*derivative_sqr_dist_dth; % df/dth

end

if (ncol >= 4)
   % compute second derivatives
   snd_derivative_ga = (-8*exp(sqr_dist./(aj^2*(-4 + sqr_dist))).*(-2*aj^6*(-8 + sqr_dist).*(-4 + sqr_dist).^3 + ...
       aj^8*(-4 + sqr_dist).^4 + 8*aj^2*(-4 + sqr_dist).^2.*sqr_dist + 8*sqr_dist.^2 + ...
       aj^4*(-4 + sqr_dist).^2.*(40 - 24*sqr_dist + sqr_dist.^2)))./((-4 + sqr_dist).^4.*(aj^3*(-4 + sqr_dist) ...
       - aj*sqr_dist).^3);
   snd_derivative_ga_alpha = (-8*exp(sqr_dist./(a_alpha^2*(-4 + sqr_dist))).*(-2*a_alpha^6*(-8 + sqr_dist).*...
       (-4 + sqr_dist).^3 + a_alpha^8*(-4 + sqr_dist).^4 + 8*a_alpha^2*(-4 + sqr_dist).^2.*sqr_dist ...
       + 8*sqr_dist.^2 + a_alpha^4*(-4 + sqr_dist).^2.*(40 - 24*sqr_dist + sqr_dist.^2)))./...
       ((-4 + sqr_dist).^4.*(a_alpha^3*(-4 + sqr_dist) - a_alpha*sqr_dist).^3);
   
   snd_derivative_gaussians = norm_cst_j*(snd_derivative_ga - (1/alpha)*snd_derivative_ga_alpha);
   
   % second derivative of sqr_dist wrt colatitude theta
   snd_derivative_sqr_dist_dthdth = 2*(cZ*cos_th_vec + (cX*cos_ph_vec + cY*sin_ph_vec).*sin_th_vec);
   
   % second derivative of sqr_dist wrt longitude phi
   snd_derivative_sqr_dist_dphdph = 2*(cX*cos_ph_vec + cY*sin_ph_vec).*sin_th_vec;
   
   % mixed derivative of sqr_dist wrt longitude phi and colatitude theta
   snd_derivative_sqr_dist_dthdph = 2*(cX*sin_ph_vec - cY*cos_ph_vec).*cos_th_vec;
   
   % snd derivative: (d^2/dph^2)f
   ff(:,4) = derivative_gaussians.*snd_derivative_sqr_dist_dphdph + ...
       snd_derivative_gaussians.*derivative_sqr_dist_dph.^2;
   
   %snd derivative: (d^2/dth^2)f
   ff(:,5) = derivative_gaussians.*snd_derivative_sqr_dist_dthdth + ...
       snd_derivative_gaussians.*derivative_sqr_dist_dth.^2;
   
   %snd derivative: (d^2/dphdth)f
   ff(:,6) = derivative_gaussians.*snd_derivative_sqr_dist_dthdph + ...
       snd_derivative_gaussians.*derivative_sqr_dist_dph.*derivative_sqr_dist_dth;
   
end

% safety for a point on the opposite side of the sphere
ind = find(abs(sqr_dist - 4) < 1e-5);
ff(ind,:) = 0;

%===================================================
