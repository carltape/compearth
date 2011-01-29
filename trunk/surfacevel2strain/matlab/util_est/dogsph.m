%
% function ff = dogsph(a,theta,phi)
% Pablo Muse, April 3, 2007
%
% computes and plots a dogsh wavelet on a longitude-colatitude grid
%

function ff = dogsph(a,theta,phi)

alpha = 1.25;
a_alpha = a*alpha;

cos_th = cos(theta);
sin_th = sin(theta);
cos_ph = cos(phi);
sin_ph = sin(phi);

X = sin_th.*cos_ph;
Y = sin_th.*sin_ph;
Z = cos_th;

cX = 0; cY = 0; cZ = 1;

sqr_dist = (X - cX).*(X - cX) + (Y - cY).*(Y - cY) + (Z - cZ).*(Z - cZ);
tan2_halfangle = sqr_dist./(4 - sqr_dist);
sqrt_lambda_a = 2*a./(2*a*a + 0.5*(1-a*a)*sqr_dist);
sqrt_lambda_a_alpha = 2*a_alpha./(2*a_alpha*a_alpha + 0.5*(1-a_alpha*a_alpha)*sqr_dist);

ff = sqrt_lambda_a.*exp(-(1./a.^2).*tan2_halfangle) - ...
    (1/alpha)*sqrt_lambda_a_alpha.*exp(-(1./a_alpha.^2).*tan2_halfangle);
