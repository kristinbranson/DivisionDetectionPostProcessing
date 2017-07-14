function [newlist,isdone] = bwconncomp_list(list,sz,varargin)

[chunksize,chunkstep,minidx,maxidx] = ...
  myparse(varargin,'chunksize',[],'chunkstep',[],'minidx',ones(size(sz)),'maxidx',sz);
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
newlist = zeros([0,5],class(list));

chunkv0 = max(min(round((list(:,1:4)-minidx)./chunkstep)+1,nchunks),1);
chunkidx = sub2indv(nchunks,chunkv0);

uniquechunkidx = unique(chunkidx);
nuniquechunks = numel(uniquechunkidx);
isdone = false(n,1);

tic;
for chunkii = 1:nuniquechunks,
  
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
      vsglobal = vs + v0 - 1;
      ctr = mean(vsglobal,1);
      fprintf('CC %d (%s) on border, mind = %f, ignoring in this chunk\n',i,num2str(ctr),...
        min(sum(abs(vsglobal-vcenter),2)));
      continue;
    end
    ismcurr = ismember(idxcurr1,cccurr.PixelIdxList{i});
    isdone(idxcurr0(ismcurr)) = true;
    scorecurr = max(list(idxcurr0(ismcurr),end));
    vsglobal = vs + v0 - 1;
    ctr = mean(vsglobal,1);
    fprintf('Adding CC %d at %s, mind = %f\n',i,mat2str(ctr),min(sum(abs(vsglobal-vcenter),2)));
    newlist(end+1,:) = [ctr,scorecurr];
%     imagesc(max(max(ischunk,[],3),[],4)); hold on; plot(ctr(2)-v0(2)+1,ctr(1)-v0(1)+1,'rx'); axis image
%     drawnow;
  end
end

%keyboard;