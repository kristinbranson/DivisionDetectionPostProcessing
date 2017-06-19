function [peaks_resize,scores,scores0,peaks] = selectLocalMaximaND(predim,predsz,varargin)

nd = numel(predsz);

[thresh,npadxyz,outsz,algorithm,chunksize,chunkstep,scorerad,minidx,maxidx] = myparse(varargin,'thresh',150,...
  'npadxyz',[16,16,37],'outsz',[2169,2048],...
  'algorithm','regionalmax','chunksize',[],'chunkstep',[],...
  'scorerad',5+zeros(1,nd),'minidx',[],'maxidx',[]);

if isempty(minidx) || isempty(maxidx),
  
  inds = ind2subv(predsz,find(predim));
  if isempty(minidx),
    minidx = min(inds,[],1);
  end
  if isempty(maxidx),
    maxidx = max(inds,[],1);
  end
  
end

switch lower(algorithm),
  
  case 'regionalmax',

    impeak = imregionalmax_sparse(predim,predsz,'chunksize',chunksize,'chunkstep',chunkstep,'thresh',thresh,...
      'minidx',minidx,'maxidx',maxidx);
    [peaks,scores0] = bwconncomp_sparse(impeak,predsz,'chunksize',chunksize,'chunkstep',chunkstep,...
      'minidx',minidx,'maxidx',maxidx);
    scores = nan(size(scores0));
    for i = 1:numel(scores0),
      
      loc = round(peaks(i,:));
      v0 = max(1,loc-scorerad);
      v1 = min(predsz,loc+scorerad);
      assert(all(v1>=v0));
      c = cell(1,nd);
      for j = 1:nd,
        c{j} = v0(j):v1(j);
      end
      idxcurr = sub2indc(predsz,c);
      scores(i) = sum(predim(idxcurr));
    end

  otherwise,

    error('Algorithm %s not implemented',algorithm);
    
end

peaks_resize = peaks;
peaks_resize(:,1:2) = (peaks(:,1:2)+npadxyz(1:2)).*outsz(2:-1:1)./(predsz(2:-1:1)+2*npadxyz(1:2));
peaks_resize(:,3) = peaks_resize(:,3) + npadxyz(3);