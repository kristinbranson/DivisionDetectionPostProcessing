function im = ConstructPredMaxProj(preddata,predshape,z0,z1,thresh)

% preddata is z,y,x
idx = true(1,size(preddata,2));
if nargin >= 3 && ~isempty(z0),
  idx = idx & preddata(1,:) >= z0;
end
if nargin >= 4 && ~isempty(z1),
  idx = idx & preddata(1,:) <= z0;
end
if nargin >= 5 && ~isempty(thresh),
  idx = idx & preddata(4,:) >= thresh;
end

im = zeros(predshape([3,2]));
for i = find(idx),
  x = preddata(3,i);
  y = preddata(2,i);
  im(x,y) = max(im(x,y),preddata(4,i));
end