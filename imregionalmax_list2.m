% [newlist,isdone] = imregionalmax_list(list,sz,...)
% 
% Performs filtering and non-maximal suppression in x, y, z, and t. To
% avoid reading the entire spatiotemporal volume into memory
% simultaneously, this function operates on sub-blocks of volume, then
% combines the results for the middle regions of these blocks. If
% startrunoncluster==true or fixrunoncluster==true, this function will
% divide up the entire volume into chunks and call a compiled version of
% this function on these chunks on the Janelia cluster. If
% finishrunoncluster==true, this function will collect the results of
% runs on chunks and combine them.
% 
% On a given spatiotemporal chunk, this function first filters with a
% 4-d filter (either box or Gaussian, depending on parameters), as true
% divisions tended to correspond to high detection values in a region,
% not just a single voxel. Then, it calls imregionalmax to find local
% maxima.
%
% Inputs:
% list: n x 5 list of detection outputs from CNN (and sparsified using
% SaveSparsePredictions). This can also be a file name, in which case it is
% loaded and overrides any other parameters input. 
% sz: 1 x 4 array indicating the size of the volume. This can also be a
% file name, in which case it is the location to save results to. This
% should only be used when list is an input file that specifies sz.
% 
% Outputs:
% newlist: If startrunoncluster==true or fixrunoncluster==true, this is a
% list of the output mat files where per-chunk results will be saved.
% Otherwise, this is an m x 5 list of detections after filtering and
% non-maximal suppression.
% isdone: if finishrunsoncluster==true, whether all jobs were done or not.
% Otherwise, this is just an empty matrix.
%
% Optional inputs:
% 
% chunksize: 1 x 4 array indicating the size of each chunk to operate on.
% This should be much bigger than the expected size of a blob of detections
% we want to reduce to a single detection via non-maximal suppression. The
% bigger this is, the more exact the results, and the less overhead in
% computation, but the more RAM is required. We used chunksize =
% [100,100,50,20]. If this is not specified, it defaults to the min of 50
% and the volume size. 
% 
% chunkstep: 1 x 4 array indicating the step in between chunks. This should
% be smaller than the chunk size, as we want some overlap between adjacent
% chunks so that we only need to use results on the interior of the chunk.
% If not specified, this will default to chunkstep =
% max(1,round(chunksize/2)). 
%
% thresh: threshold on the filtered detection scores. If not specified (or
% empty), no thresholding will be done. This is the default behavior.
%
% filrad: 1 x 4 vector indicating the radius of the filter. Default:
% filsig (either filsig or filrad must be specified for filtering to
% happen). 
%
% filsig: 1 x 4 vector indicating the standard deviation of the filter. If
% this is specified (or non-empty), Gaussian filtering is used, otherwise
% box filtering is used. Default is []. 
%
% chunki0, chunki1: Range of chunks to run filtering and non-max
% suppression on. Default: 1, inf. 
% 
% startrunoncluster: If this flag is true, then this will split up the
% volume into chunks and run the chunks on the cluster (nchunksperjob at a
% time). Default = false
% 
% fixrunoncluster: If this flag is true, then this will split up the
% volume into chunks and run the chunks on the cluster (nchunksperjob at a
% time). This differs from startrunoncluster in that it will only start
% jobs that do not appear to be completed previously. Default = false.
%
% finishrunoncluster: If this flag is true, then this will collect and
% combine results of jobs run on the cluster. Default = false
%
% nchunksperjob: Number of chunks to run in each job (see
% startrunoncluster). Default=100.
%
% inmatfile: Required for starting and finishing runs on the cluster.
% Name of file to store information about cluster runs. Default=''
%
% outdir: Required for starting runs on the cluster. Directory to save job
% results to. Default = ''. 
% 
% outmatfilestr, scriptfilestr, logfilestr: File roots for temporary files
% created for each job. Default = 'RegionalMax'
% 
% tmprootdir: where to store the MCR cache. Only needed for running jobs on
% the cluster. default = /scratch/<username>
% 
% script: location of compiled version of this code
% 
% mcr: location of MCR
%
% ncores: max number of cores to use when deployed. 
function [newlist,isdone] = imregionalmax_list2(list,sz,varargin)

newlist = [];
isdone = [];

% version that runs on the cluster
if ischar(list),
  
  inmatfile = list;
  outmatfile = sz;
  load(inmatfile);

  chunki0 = varargin{1};
  chunki1 = varargin{2};
  if isdeployed,
    if ischar(chunki0),
      chunki0 = str2double(chunki0);
    end
    if ischar(chunki1),
      chunki1 = str2double(chunki1);
    end
  end
    
  isjob = true;
  
  if isdeployed,
    maxNumCompThreads(ncores);
  end
  

else  
  
  isjob = false;

  username = GetUserName();
  [chunksize,chunkstep,thresh,minidx,maxidx,filrad,filsig,chunki0,chunki1,...
    startrunoncluster,finishrunoncluster,fixrunoncluster,nchunksperjob,inmatfile,outdir,...
    outmatfilestr,scriptfilestr,logfilestr,...
    TMP_ROOT_DIR,MCR_CACHE_ROOT,SCRIPT,MCR,ncores] = ...
    myparse(varargin,'chunksize',[],'chunkstep',[],'thresh',[],'minidx',ones(size(sz)),'maxidx',sz,...
    'filrad',[],'filsig',[],...
    'chunki0',1,'chunki1',inf,...
    'startrunoncluster',false,'finishrunoncluster',false,...
    'fixrunoncluster',false,...
    'nchunksperjob',100,...
    'inmatfile','','outdir','',...
    'outmatfilestr','RegionalMax','scriptfilestr','RegionalMax','logfilestr','RegionalMax',...
    'tmprootdir',fullfile('/scratch',username),...
    'mcrcacheroot','',...
    'script','/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/imregionalmax_list/for_redistribution_files_only/run_imregionalmax_list.sh',...
    'mcr','/groups/branson/bransonlab/share/MCR/v91',...
    'ncores',1);
  
  if isempty(MCR_CACHE_ROOT),
    MCR_CACHE_ROOT = fullfile(TMP_ROOT_DIR,'mcr_cache_root');
  end

  if finishrunoncluster,
    load(inmatfile);
    
    jobisdone = cellfun(@(x) exist(x,'file'), outmatfiles);
    if ~all(jobisdone>0),
      fprintf('%d / %d jobs complete\n',nnz(jobisdone),njobs);
      for i = 1:njobs,
        fprintf('%s: %d\n',outmatfiles{i},jobisdone(i)>0);
      end
      return;
    end
    
    newlist = zeros([0,5],class(list));
    isdone = false(n,1);
    for i = 1:njobs,
      res = load(outmatfiles{i});
      isnew = ~ismember(res.newlist,newlist,'rows');
      nnew = nnz(isnew);
      newlist(end+1:end+nnew,:) = res.newlist(isnew,:);
      isdone = isdone | res.isdone;
    end
    
    return;
  end
  
  nd = numel(sz);
  n = size(list,1);

  if isempty(chunksize),
    chunksize = min(sz,repmat(50,[1,nd]));
  end
  if isempty(chunkstep),
    chunkstep = max(1,round(chunksize/2));
  end
  assert(all(chunkstep<=chunksize));
  
  dofil = ~isempty(filsig) || ~isempty(filrad);
  if dofil,
    if isempty(filrad),
      filrad = filsig;
    end
    chunksizefil = chunksize+filrad*2;
    imchunkfil = zeros(chunksizefil,class(list));
    if isempty(filsig),
      fil = ones(filrad*2+1,class(list))/prod(2*filrad+1);
    else
      [xgrid,ygrid,zgrid,tgrid] = ndgrid(-filrad(1):filrad(1),-filrad(2):filrad(2),-filrad(3):filrad(3),-filrad(4):filrad(4));
      fil = mvnpdf([xgrid(:),ygrid(:),zgrid(:),tgrid(:)],zeros(1,4),diag(filsig));
      fil = reshape(fil,size(xgrid));
      fil = cast(fil,class(list));
    end
  end
  dothresh = ~isempty(thresh);

  nidx = maxidx - minidx + 1;

  nchunks = round((nidx-1)./chunkstep)+1;
  %nchunkstotal = prod(nchunks);
  
  ismiddlechunk = false(chunksize);
  s = cell(1,nd);
  for d = 1:nd,
    r = ceil(chunkstep(d)/2);
    c = round(chunksize(d)/2);
    s{d} = max(1,c-r):min(chunksize(d),c+r);
  end
  ismiddlechunk(s{:}) = true;
  
  % s = cell(1,nd);
  % [s{:}] = deal(-1:1);
  % [s{:}] = ndgrid(s{:});
  
  chunkv0 = round((list(:,1:4)-minidx)./chunkstep)+1;
  chunkidx = sub2indv(nchunks,chunkv0);
  
  uniquechunkidx = unique(chunkidx);
  nuniquechunks = numel(uniquechunkidx);
  
  if startrunoncluster || fixrunoncluster,
    
    assert(~isempty(inmatfile));
    if ~exist(outdir,'dir'),
      mkdir(outdir);
    end
    assert(exist(outdir,'dir')>0);
    
    njobs = ceil(nuniquechunks/nchunksperjob);
    outmatfiles = cell(1,njobs);
    for i = 1:njobs,
      chunki0 = (i-1)*nchunksperjob + 1;
      chunki1 = min(nuniquechunks,i*nchunksperjob);
      outmatfiles{i} = fullfile(outdir,[outmatfilestr,sprintf('_chunk%dto%d.mat',chunki0,chunki1)]);
    end
    clear varargin;
    save(inmatfile);
    
    for i = 1:njobs,
      chunki0 = (i-1)*nchunksperjob + 1;
      chunki1 = min(nuniquechunks,i*nchunksperjob);
      outmatfile = outmatfiles{i};
      if fixrunoncluster && exist(outmatfile,'file'),
        continue;
      end
      scriptfile = fullfile(outdir,[scriptfilestr,sprintf('_chunk%dto%d.sh',chunki0,chunki1)]);
      logfile = fullfile(outdir,[logfilestr,sprintf('_chunk%dto%d.log',chunki0,chunki1)]);
      jobid = sprintf('RegionalMax_chunk%dto%d_%d',chunki0,chunki1,i);
      fid = fopen(scriptfile,'w');
      fprintf(fid,'if [ -d %s ]\n',TMP_ROOT_DIR);
      fprintf(fid,'  then export MCR_CACHE_ROOT=%s.%s\n',MCR_CACHE_ROOT,jobid);
      fprintf(fid,'fi\n');
      fprintf(fid,'%s %s %s %s %d %d\n',...
        SCRIPT,MCR,inmatfile,outmatfile,chunki0,chunki1);
      fclose(fid);
      
      unix(sprintf('chmod u+x %s',scriptfile));
      cmd = sprintf('ssh login1 ''source /etc/profile; bsub -n %d -J %s -o ''%s'' ''\"%s\"''''',...
        ncores,jobid,logfile,scriptfile);
      unix(cmd);

    end
    
    newlist = outmatfiles;
    
    return;
  end
  
end

chunki1 = min(chunki1,nuniquechunks);


imchunk = zeros(chunksize,class(list));
newlist = zeros([0,5],class(list));

isdone = false(n,1);

for chunkii = chunki0:chunki1,
  
  chunki = uniquechunkidx(chunkii);  
  chunkv = ind2subv(nchunks,chunki);
  v0 = (chunkv-2).*chunkstep+minidx;
  v1 = v0+chunksize-1;
  
  idxcurr0 = true(n,1);
  for d = 1:nd,
    idxcurr0 = idxcurr0 & v0(d) <= list(:,d) & v1(d) >= list(:,d);
  end

  fprintf('Chunk %d: %s to %s, %d non-zero (%d / %d)\n',...
    chunki,mat2str(v0),mat2str(v1),nnz(idxcurr0),chunkii,nuniquechunks);
  
  if ~any(idxcurr0),
    continue;
  end
  
  idxcurr1 = sub2indv(chunksize,list(idxcurr0,1:nd)-v0+1);
  isdone(idxcurr0) = isdone(idxcurr0) | ismiddlechunk(idxcurr1);
  
  if dofil,
    v0fil = v0 - filrad;
    v1fil = v1 + filrad;
    idxcurr0fil = true(n,1);
    for d = 1:nd,
      idxcurr0fil = idxcurr0fil & v0fil(d) <= list(:,d) & v1fil(d) >= list(:,d);
    end
    imchunkfil(:) = 0;
    idxcurr1fil = sub2indv(chunksizefil,list(idxcurr0fil,1:nd)-v0fil+1);
    imchunkfil(idxcurr1fil) = list(idxcurr0fil,end);
    imchunk = convn(imchunkfil,fil,'valid');

  end
      
  if ~dofil,
    idxcurr1 = sub2indv(chunksize,list(idxcurr0,1:nd)-v0+1);
    imchunk(:) = 0;
    imchunk(idxcurr1) = list(idxcurr0,end);
  end
  
  if dothresh,
    imchunk(imchunk<thresh) = 0;
  end
  
  ischunk = ismiddlechunk & imregionalmax(imchunk);
  idxcurr2 = ind2subv(chunksize,find(ischunk)) + v0 - 1;
  ncurr = size(idxcurr2,1);
  newlist(end+1:end+ncurr,:) = cat(2,idxcurr2,imchunk(ischunk));

  fprintf('%d / %d non-zero pixels \n',ncurr,nnz(imchunk(ismiddlechunk)));

end

if isjob,
  save(outmatfile,'newlist','isdone');
end
  

