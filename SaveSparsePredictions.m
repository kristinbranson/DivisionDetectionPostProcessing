% [allpredidx,densedatadate] = SaveSparsePredictions(infile,outmatfile,predtimestamp,predthresh,ncores)
% Reads in an hdf5 file of detection maps (with a score for each voxel
% location in the volume), thresholds, and saves the sparse list of
% locations of voxels above threshold and their scores to a file.
%
% Inputs: 
% infile: name of hdf5 file of detection maps (with a score for each voxel
% location in the volume)
% outmatfile: name of .mat file to save sparse detection list to
% predtimestamp: which timestamp this volume corresponds to
% predthresh: threshold above which we will save the detections. Too low a
% value will result in bigger files, too high will result in missed
% detections. We used .01. 
% ncores: when run in compiled mode, number of cores to restrict the
% process to. not sure MATLAB actually listens to this for everything...
%
% Outputs:
% allpredidx: sparse list of spatiotemporal locations and scores of
% predictions above threshold. this will be ndetections x 5, and columns
% correspond to x, y, z, t, score
% densedatadate: date stored in the hdf5 metadata

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