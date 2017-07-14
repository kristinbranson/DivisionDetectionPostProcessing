function [allpredidx,densedatadate] = SaveSparsePredictions(infile,outmatfile,predtimestamp,predthresh,ncores)

preddatasetname = '/predictions';

if isdeployed,
  if ischar(predthresh),
    predthresh = str2double(predthresh);
  end
  if ischar(predtimestamp),
    predtimestamp = str2double(predtimestamp);
  end
  
  if exist('ncores','var'),
    if ischar(ncores),
      ncores = str2double(ncores);
    end
    maxNumCompThreads(ncores);
  end
end

tmp = dir(infile);
densedatadate = tmp.datenum;
preddata = h5read(infile,preddatasetname);

idxthresh = find(preddata>=predthresh);
szcurr = size(preddata);
locscurr = ind2subv(szcurr,idxthresh);
scorescurr = preddata(idxthresh);
fprintf('%s, nthresh = %d / %d\n',infile,numel(idxthresh),numel(preddata));
allpredidx = cat(2,locscurr,repmat(predtimestamp,[numel(idxthresh),1]),scorescurr);

save(outmatfile,'allpredidx','predtimestamp','infile','densedatadate');