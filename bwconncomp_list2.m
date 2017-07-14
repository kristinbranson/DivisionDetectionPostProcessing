function [newlist,isdone] = bwconncomp_list2(list,sz,varargin)

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
  
  [chunksize,chunkstep,minidx,maxidx,...
    chunki0,chunki1,...
    startrunoncluster,finishrunoncluster,nchunksperjob,inmatfile,outdir,...
    outmatfilestr,scriptfilestr,logfilestr,...
    TMP_ROOT_DIR,MCR_CACHE_ROOT,SCRIPT,MCR,ncores] = ...
    myparse(varargin,'chunksize',[],'chunkstep',[],'minidx',ones(size(sz)),'maxidx',sz,...
    'chunki0',1,'chunki1',inf,...
    'startrunoncluster',false,'finishrunoncluster',false,...
    'nchunksperjob',100,...
    'inmatfile','','outdir','',...
    'outmatfilestr','ConnComp','scriptfilestr','ConnComp','logfilestr','ConnComp',...
    'tmprootdir',fullfile('/scratch',username),...
    'mcrcacheroot','',...
    'script','/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/bwconncomp_list/for_redistribution_files_only/run_bwconncomp_list.sh',...
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
  
  nidx = maxidx - minidx + 1;
  
  nchunks = ceil(nidx./chunkstep);
  %nchunkstotal = prod(nchunks);
  
  ischunk = false(chunksize);
  
  chunkv0 = max(min(round((list(:,1:4)-minidx)./chunkstep)+1,nchunks),1);
  chunkidx = sub2indv(nchunks,chunkv0);
  
  uniquechunkidx = unique(chunkidx);
  nuniquechunks = numel(uniquechunkidx);
  
  if startrunoncluster,
    
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
      scriptfile = fullfile(outdir,[scriptfilestr,sprintf('_chunk%dto%d.sh',chunki0,chunki1)]);
      logfile = fullfile(outdir,[logfilestr,sprintf('_chunk%dto%d.log',chunki0,chunki1)]);
      jobid = sprintf('ConnComp_chunk%dto%d_%d',chunki0,chunki1,i);
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
newlist = zeros([0,5],class(list));
isdone = false(n,1);

tic;
for chunkii = chunki0:chunki1,
  
  if toc > 10,
    fprintf('Chunk %d / %d\n',chunkii,nuniquechunks);
    tic;
  end
  
  chunki = uniquechunkidx(chunkii);  
  chunkv = ind2subv(nchunks,chunki);
  v0 = (chunkv-2).*chunkstep+minidx;
  v1 = v0+chunksize-1;
  
  idxcurr0 = ~isdone;
  for d = 1:nd,
    idxcurr0 = idxcurr0 & v0(d) <= list(:,d) & v1(d) >= list(:,d);
  end

%   fprintf('Chunk %d: %s to %s, %d non-zero (%d / %d)\n',...
%     chunki,mat2str(v0),mat2str(v1),nnz(idxcurr0),chunkii,nuniquechunks);
  
  if ~any(idxcurr0),
    continue;
  end
    
  idxcurr1 = sub2indv(chunksize,list(idxcurr0,1:nd)-v0+1);
  ischunk(:) = 0;
  ischunk(idxcurr1) = 1;

  idxcurr0 = find(idxcurr0);
  
  cccurr = bwconncomp(ischunk);

  islast = v1 >= maxidx;
  isfirst = v0 <= minidx;  
  
  for i = 1:cccurr.NumObjects,
    
    vs = ind2subv(chunksize,cccurr.PixelIdxList{i});
    isborder = any(any(~islast & ( (~isfirst & vs==1) | vs == chunksize)));
    if isborder,
      %fprintf('CC %d on border, ignoring in this chunk\n',i);
      continue;
    end
    ismcurr = ismember(idxcurr1,cccurr.PixelIdxList{i});
    isdone(idxcurr0(ismcurr)) = true;
    scorecurr = max(list(idxcurr0(ismcurr),end));
    vsglobal = vs + v0 - 1;
    ctr = mean(vsglobal,1);
    %fprintf('Adding CC %d at %s\n',i,mat2str(ctr));
    newlist(end+1,:) = [ctr,scorecurr];
%     imagesc(max(max(ischunk,[],3),[],4)); hold on; plot(ctr(2)-v0(2)+1,ctr(1)-v0(1)+1,'rx'); axis image
%     drawnow;
  end
end

if isjob,
  save(outmatfile,'newlist','isdone');
end
  


%keyboard;