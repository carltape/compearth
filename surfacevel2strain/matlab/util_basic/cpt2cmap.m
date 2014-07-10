function outmat = cpt2cmap(inmat)
%CPT2CMAP input a nx3 matrix of RGB colors and output a matrix suitable for colormap

len = length(inmat);
temp = inmat/max(max(inmat));  % normalize to 0-1 range
rv = temp(:,1);
gv = temp(:,2);
bv = temp(:,3);

% there is probably a better way to do this -- this interpolates to numc points.
numc = 65;
outmat = [  interp1(linspace(1,numc,len), rv', [1:numc])' ...
            interp1(linspace(1,numc,len), gv', [1:numc])' ...
            interp1(linspace(1,numc,len), bv', [1:numc])'];
