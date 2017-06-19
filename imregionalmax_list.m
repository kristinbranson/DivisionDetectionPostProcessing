function newlist = imregionalmax_list(list,sz,varargin)

[chunksize,chunkstep,thresh,minidx,maxidx,filrad] = ...
  myparse(varargin,'chunksize',[],'chunkstep',[],'thresh',[],'minidx',ones(size(sz)),'maxidx',sz,'filrad',[]);
nd = numel(sz);
n = size(list,1);

if isempty(chunksize),
  chunksize = min(sz,repmat(50,[1,nd]));
end
if isempty(chunkstep),
  chunkstep = max(1,round(chunksize/2));
end
assert(all(chunkstep<=chunksize));

dofil = ~isempty(filrad);
if dofil,
  chunksizefil = chunksize+filrad*2;
  imchunkfil = zeros(chunksizefil,class(list));
  fil = ones(filrad*2+1,class(list))/prod(2*filrad+1);
end
dothresh = ~isempty(thresh);

nidx = maxidx - minidx + 1;

nchunks = ceil(nidx./chunkstep);
%nchunkstotal = prod(nchunks);

imchunk = zeros(chunksize,class(list));
newlist = zeros([0,5],class(list));

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
isdone = false(n,1);

for chunkii = 1:nuniquechunks,
  
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

