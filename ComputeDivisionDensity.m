% [maxprojdiv,zdiv] = ComputeDivisionDensity(list,t,sz,...)
%
% Inputs a list of divisions, creates a 4-d volume around a given time
% point, convolves with a 4-d Gaussian filter to compute density of
% divisions at each voxel, and then computes a max projection in a
% specified dimension.
%
% Inputs:
%
% list: n x 5 list of detection outputs from bwconncomp_list2. This can
% also be a file name, in which case it is loaded and overrides any other
% parameters input. 
% t: time point to compute division density for. 
% sz: 1 x 4 array indicating the size of the volume. This can also be a
% file name, in which case it is the location to save results to. This
% should only be used when list is an input file that specifies sz.
% 
% Outputs:
% 
% maxprojdiv: 2-d maximum-intensity projection of division density. 
% zdiv: 3-d division density before maximum intensity projection. 
%
% Optional inputs:
%
% filsig: Standard deviation os Gaussian filter. Default = [50,50,50,1].
% filrad: Radius of Gaussian filter. Default = ceil(filsig*2.5).
% thresh: Threshold on division score. 
% stepsz: Gaussian blurring and maximum intensity projection is done in
% chunks, stepsz specifies the size of chunks. Default = [110,110,50,3].
% dim: Which dimension to compute the maximum intensity projection over. 
%
function [maxprojdiv,zdiv] = ComputeDivisionDensity(list,t,sz,varargin)

isjob = ischar(list);
if isjob,
  inmatfile = list;
  realt = t;
  if ischar(realt),
    realt = str2double(realt);
  end
  outmatfile = sz;
  
  if isempty(varargin),
    out3dfile = '';
  else
    out3dfile = varargin{1};
  end
  load(inmatfile);
  t = realt;
  
  fprintf('Arguments:\n');
  fprintf('inmatfile: %s\n',inmatfile);
  fprintf('outmatfile: %s\n',outmatfile);
  fprintf('t: %d\n',t);
  fprintf('out3dfile: %s\n',out3dfile);
    
else

  [filsig,filrad,thresh,stepsz,ncores,dim,out3dfile] = ...
    myparse(varargin,'filsig',[50,50,50,1],'filrad',[],'thresh',0.1465,...
    'stepsz',[110,110,50,3],'ncores',[],'dim',3,'out3dfile','');
  
end

if isdeployed,
  maxNumCompThreads(ncores);
end

if isempty(filrad),
  filrad = ceil(filsig*2.5);
end

fils = cell(1,4);
for d = 1:4,
  fils{d} = normpdf((-filrad(d):filrad(d))',0,filsig(d));
  fils{d} = fils{d} / sum(fils{d}(:));
%   tmpsz = ones(1,4);
%   tmpsz(d) = numel(fils{d});
%   fils{d} = reshape(fils{d},tmpsz);
end

% [xgrid,ygrid,zgrid,tgrid] = ...
%   ndgrid(-filrad(1):filrad(1),...
%   -filrad(2):filrad(2),...
%   -filrad(3):filrad(3),...
%   -filrad(4):filrad(4));
% fil = mvnpdf([xgrid(:),ygrid(:),zgrid(:),tgrid(:)],zeros(1,4),diag(filsig));
% fil = reshape(fil,size(xgrid));


list(:,1:4) = round(list(:,1:4));

szcurr = [sz(1:3),filrad(4)*2+1];
nsteps = ceil(szcurr./stepsz);
nstepstotal = prod(nsteps);

t0 = t-filrad(4);
t1 = t+filrad(4);
idxcurr = list(:,4)>=t0 & list(:,4)<=t1 & list(:,5) >= thresh;
list = list(idxcurr,:);

tic;
otherdims = setdiff(1:3,dim);
maxprojdiv = zeros(sz(otherdims));
zdiv = zeros(sz(1:2));
for stepi = 1:nstepstotal,
  if toc > 10,
    fprintf('Step %d / %d\n',stepi,nstepstotal);
%     imagesc(maxprojdiv);
%     drawnow;
    tic;
  end
  stepsub = ind2subv(nsteps,stepi);
  v0 = (stepsub-1).*stepsz + 1;
  v0(4) = v0(4) + t0 - 1;
  v00 = v0;
  v0 = max(1,v0-filrad-1);
  v1 = min(stepsub.*stepsz,szcurr);
  v10 = v1;
  v1 = min(szcurr,v1+filrad+1);
  v1(4) = v1(4) + t0 - 1;
  idxcurr1 = all(v0 <= list(:,1:4),2) & ...
    all(v1 >= list(:,1:4),2);
  if ~any(idxcurr1),
    continue;
  end
%   divcurr = zeros(v1-v0+1);
%   idx = sub2indv(v1-v0+1,list(idxcurr1,1:4)-v0+1);
%   divcurr(idx) = 1; 
%   tic;
%   for d = 4:-1:1,
%     divcurr = convn(divcurr,fils{d},'same');
%   end
%   toc;
  
  div3 = zeros(v1-v0+1);
  idx = sub2indv(v1-v0+1,list(idxcurr1,1:4)-v0+1);
  div3(idx) = 1;
  for d = 1:4,
    order = [d,setdiff(1:4,d)];
    [~,unorder] = sort(order);
    div3 = permute(div3,order);
    for ti = 1:size(div3,4),
      for z = 1:size(div3,3),
        div3(:,:,z,ti) = conv2(div3(:,:,z,ti),fils{d},'same');
      end
    end
    div3 = permute(div3,unorder);
  end
  
%   divcurr = zeros(v1-v0+1);
%   idx = sub2indv(v1-v0+1,list(idxcurr1,1:4)-v0+1);
%   divcurr(idx) = 1;
%   tic;
%   divcurr = convn(divcurr,fil,'same');
%   toc;

  [maxprojdivcurr,zcurr] = max(div3(:,:,:,filrad(4)+1),[],dim);
  maxprojdivcurr = permute(maxprojdivcurr,[otherdims,dim]);
  zcurr = permute(zcurr,[otherdims,dim]);
  tmp = maxprojdiv(v0(otherdims(1)):v1(otherdims(1)),v0(otherdims(2)):v1(otherdims(2)));
  ztmp = zdiv(v0(otherdims(1)):v1(otherdims(1)),v0(otherdims(2)):v1(otherdims(2)));
  isbigger = maxprojdivcurr>=tmp;
  tmp(isbigger) = maxprojdivcurr(isbigger);
  ztmp(isbigger) = zcurr(isbigger)+v0(3)-1;
  maxprojdiv(v0(otherdims(1)):v1(otherdims(1)),v0(otherdims(2)):v1(otherdims(2))) = tmp;
  zdiv(v0(otherdims(1)):v1(otherdims(1)),v0(otherdims(2)):v1(otherdims(2))) = ztmp;
  
  if ~isempty(out3dfile),
    sz = v1-v0+1;
    i0 = min(max(1,v00-v0+1),sz);
    i1 = min(max(1,v10-v0+1),sz);
    ncurr = prod(i1(1:3)-i0(1:3)+1);
    maxd = max(reshape(div3(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),filrad(4)+1),[ncurr,1]));
    if maxd > 0,
      div3_uint16 = uint16(div3(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),filrad(4)+1)/maxd*(2^16-1));
      sz = size(div3_uint16);
      [p,n,ext] = fileparts(out3dfile);
      out3dfile2 = fullfile(p,sprintf('%s_%d%s',n,stepi,ext));
      switch ext,
        case '.mat'
          save(out3dfile2,'div3_uint16','maxd','v0','v1');
        case '.klb',
          if isempty(ncores),
            numThreads = -2;
          else
            numThreads = ncores;
          end
          pixelSize = -1;
          blockSize = [sz(1),sz(2),1];
          compressionType = 1;
          KLB_METADATA_SIZE = 256;
          metadata = sprintf('maxd:%e,v0:%s,v1:%s',maxd,mat2str(v00(1:3)),mat2str(v10(1:3)));
          metadata = [metadata,repmat(' ',[1,KLB_METADATA_SIZE-numel(metadata)])]; %#ok<AGROW>
          writeKLBstack(div3_uint16,out3dfile2,numThreads,pixelSize,blockSize,compressionType,metadata);
      end
    end
  end
  
end

if isjob,
  save(outmatfile,'maxprojdiv','zdiv');
end
