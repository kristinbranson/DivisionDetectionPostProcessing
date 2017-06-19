function im = imregionalmax_sparse(im,sz,varargin)

[chunksize,chunkstep,thresh,minidx,maxidx] = ...
  myparse(varargin,'chunksize',[],'chunkstep',[],'thresh',[],'minidx',ones(size(sz)),'maxidx',sz);
nd = numel(sz);

if isempty(chunksize),
  chunksize = min(sz,repmat(50,[1,nd]));
end
if isempty(chunkstep),
  chunkstep = round(chunksize/2);
end
assert(all(chunkstep<=chunksize));

nidx = maxidx - minidx + 1;

nchunks = ceil(nidx./chunkstep);
nchunkstotal = prod(nchunks);

for chunki = 1:nchunkstotal,
  
  chunkv = ind2subv(nchunks,chunki);
  
  v0 = (chunkv-1).*chunkstep+minidx;
  v1 = min(v0+chunksize-1,maxidx);
  ncurr = v1-v0+1;
  c = cell(1,nd);
  for j = 1:nd,
    c{j} = v0(j):v1(j);
  end
  idxcurr = sub2indc(sz,c);
  imchunk = im(idxcurr);
  if ~isempty(thresh),
    imchunk(imchunk < thresh) = 0;
  end
  n0 = nnz(imchunk);
  if ~any(imchunk>0),
    continue;
  end
  fprintf('Chunk %d / %d: %s / %s\n',chunki,nchunkstotal,mat2str(chunkv),mat2str(nchunks));
  imchunk = reshape(full(imchunk),ncurr);
  ischunk = imregionalmax(imchunk);
  imchunk(~ischunk) = 0;
  n1 = nnz(imchunk);
  fprintf('%d -> %d non-zero pixels\n',n0,n1);
  im(idxcurr) = imchunk;
  
end


%%

% for d1 = 1:nd,
%   for d2 = d1+1:nd,
% 
%     if isinit && d1 == 1 && d2 == 2,
%       continue;
%     end
%     
%     dcurr = 1:nd;
%     dcurr([d1,d2]) = [];
%     szcurr = sz(dcurr);
%     nelcurr = prod(szcurr);
%     idxadd = [];
%     vadd = [];
%     [cgrid,rgrid] = meshgrid(1:sz(d2),1:sz(d1));
%     rgrid = rgrid(:);
%     cgrid = cgrid(:);
%     subvcurr = nan([sz(d1)*sz(d2),nd]);
%     subvcurr(:,d1) = rgrid;
%     subvcurr(:,d2) = cgrid;
%     for i = 1:nelcurr,
%       if mod(i,100) == 0,
%         fprintf('d1 = %d, d2 = %d, i = %d\n',d1,d2,i);
%       end
%       tmp = ind2subv(sz(dcurr),i);
%       for j = 1:numel(tmp),
%         subvcurr(:,dcurr(j)) = tmp(j);
%       end
%       idxcurr = sub2indv(sz,subvcurr);
%       if ~any(im(idxcurr)),
%         continue;
%       end
%       imslice = full(reshape(im(idxcurr),sz([d1,d2])));
%       ismaxslice = imregionalmax(imslice);
%       % first time
%       if ~(d1 == 1 && d2 == 2),
%         ismaxslice = ismaxslice & vmaxmat(idxcurr)>0;
%       end
%       ncurr = nnz(ismaxslice);
%       idxadd(end+1:end+ncurr) = idxcurr(ismaxslice);
%       vadd(end+1:end+ncurr) = imslice(ismaxslice);
% 
%     end
%     
%     assert(numel(idxadd) == numel(unique(idxadd)));
%     vmaxmat(idxadd,1) = vadd;
%     
%   end
% end

