function im = ConstructPredMaxProj2(preddata,predshape,t,varargin)

[xlim,ylim,zlim,thresh] = myparse(varargin,'xlim',[],'ylim',[],'zlim',[],'thresh',[]);

% preddata is x,y,z,t
idx = preddata(:,4) == t;
if numel(zlim) >= 1 && ~isnan(zlim(1)),
  idx = idx & preddata(:,3) >= zlim(1);
else
  zlim(1) = 1;
end
if numel(zlim) >= 2 && ~isnan(zlim(2)),
  idx = idx & preddata(:,3) <= zlim(2);
else
  zlim(2) = predshape(3);
end

if numel(xlim) >= 1 && ~isnan(xlim(1)),
  idx = idx & preddata(:,1) >= xlim(1);
else
  xlim(1) = 1;
end
if numel(xlim) >= 2 && ~isnan(xlim(2)),
  idx = idx & preddata(:,1) <= xlim(2);
else
  xlim(2) = predshape(1);
end

if numel(ylim) >= 1 && ~isnan(ylim(1)),
  idx = idx & preddata(:,2) >= ylim(1);
else
  ylim(1) = 1;
end
if numel(ylim) >= 2 && ~isnan(ylim(2)),
  idx = idx & preddata(:,2) <= ylim(2);
else
  ylim(2) = predshape(2);
end


if ~isempty(thresh),
  idx = idx & preddata(:,5) >= thresh;
end

im = zeros([xlim(2)-xlim(1)+1,ylim(2)-ylim(1)+1]);
for i = find(idx'),
  x = preddata(i,1) - xlim(1) + 1;
  y = preddata(i,2) - ylim(1) + 1;
  im(x,y) = max(im(x,y),preddata(i,5));
end