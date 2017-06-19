function im = ConstructPredIm(z,preddata,predshape)

% preddata is z,y,x
im = zeros(predshape([3,2]));
idx = find(round(preddata(1,:))==z);
im(sub2ind(predshape([3,2]),preddata(3,idx),preddata(2,idx))) = preddata(4,idx);