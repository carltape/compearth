function map = sar(m)
% SAR, a wrapped around sar colormap.

if nargin + nargout == 0
   h = get(gca,'child');
   m = length(h);
elseif nargin == 0
   m = size(get(gcf,'colormap'),1);
end

R=zeros(16,3);
R(:,1)=ones(16,1);
for i=1:8;
	R(i,2)=1-i*1/8;
end;
for i=9:16;
	R(i,2)=(i-8)*1/8;
end;
clear i;

e = ones(ceil(m/16),1);
R = kron(e,R);
R = R(1:m,:);
R(1,:)=[.5 .5 .5];
R(length(R),:)=[.5 .5 .5];

if nargin + nargout == 0
   % Apply to lines in current axes.
   for k = 1:m
      if strcmp(get(h(k),'type'),'line')
         set(h(k),'color',R(k,:))
      end
   end
else
   map = R;
end
