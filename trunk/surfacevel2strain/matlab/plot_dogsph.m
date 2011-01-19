%
% function [ff,ff2] = plot_dogsph(a0,q,ntheta,nphi)
% Pablo Muse, April 3, 2007
%
% computes and plots a dogsh wavelet on a longitude-colatitude grid
%

function ff = plot_dogsph(a0,q,ntheta,nphi)

alpha = 1.25;

%a0 = 1;
aj = a0./ 2.^q;
a_alpha = alpha*aj;


norm_cst_j = 1; % l2 normalization

thetas = linspace(0, pi, ntheta);
phis = linspace(-pi, pi, nphi);
%ff = zeros(ntheta,nphi);
ff = zeros(1,ntheta);
cX = 0; cY = 0; cZ = 1;
for i= 1:ntheta
    
    cos_th = cos(thetas(i));
    sin_th = sin(thetas(i));

    sin_ph = 0; cos_ph = 1;
    X = sin_th*cos_ph;
    Y = sin_th*sin_ph;
    Z = cos_th;

    sqr_dist = (X - cX).*(X - cX) + (Y - cY).*(Y - cY) + (Z - cZ).*(Z - cZ);
    tan2_halfangle = sqr_dist./(4 - sqr_dist);
    sqrt_lambda_a = 2*aj./(2*aj*aj + 0.5*(1-aj*aj)*sqr_dist);
    sqrt_lambda_a_alpha = 2*a_alpha./(2*a_alpha*a_alpha + 0.5*(1-a_alpha*a_alpha)*sqr_dist);
    ff(i) = sqrt_lambda_a.*exp(-(1./aj.^2).*tan2_halfangle) - ...
        (1/alpha)*sqrt_lambda_a_alpha.*exp(-(1./a_alpha.^2).*tan2_halfangle);
 
end


%mesh(thetas,phis,ff)

subplot(1,2,1); plot(thetas,ff,'r'); grid
subplot(1,2,2); polar(pi/2 - thetas, 1 + ff,'r');