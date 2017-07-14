function [maxprojdiv,zdiv] = ComputeDivisionDensity(list,t,sz,varargin)

isjob = ischar(list);
if isjob,
  inmatfile = list;
  realt = t;
  if ischar(realt),
    realt = str2double(realt);
  end
  outmatfile = sz;
  load(inmatfile);
  t = realt;
    
else

  [filsig,filrad,thresh,stepsz,ncores,dim] = myparse(varargin,'filsig',[50,50,50,1],'filrad',[],'thresh',0.1465,...
    'stepsz',[110,110,50,3],'ncores',[],'dim',3);
  
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
  v0 = max(1,v0-filrad-1);
  v1 = stepsub.*stepsz;
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
  
  divcurr = zeros(v1-v0+1);
  idx = sub2indv(v1-v0+1,list(idxcurr1,1:4)-v0+1);
  divcurr(idx) = 1;
  for d = 1:4,
    order = [d,setdiff(1:4,d)];
    [~,unorder] = sort(order);
    divcurr = permute(divcurr,order);
    for ti = 1:size(divcurr,4),
      for z = 1:size(divcurr,3),
        divcurr(:,:,z,ti) = conv2(divcurr(:,:,z,ti),fils{d},'same');
      end
    end
    divcurr = permute(divcurr,unorder);
  end
  
%   divcurr = zeros(v1-v0+1);
%   idx = sub2indv(v1-v0+1,list(idxcurr1,1:4)-v0+1);
%   divcurr(idx) = 1;
%   tic;
%   divcurr = convn(divcurr,fil,'same');
%   toc;
  
  [maxprojdivcurr,zcurr] = max(divcurr(:,:,:,filrad(4)+1),[],dim);
  maxprojdivcurr = permute(maxprojdivcurr,[otherdims,dim]);
  zcurr = permute(zcurr,[otherdims,dim]);
  tmp = maxprojdiv(v0(otherdims(1)):v1(otherdims(1)),v0(otherdims(2)):v1(otherdims(2)));
  ztmp = zdiv(v0(otherdims(1)):v1(otherdims(1)),v0(otherdims(2)):v1(otherdims(2)));
  isbigger = maxprojdivcurr>=tmp;
  tmp(isbigger) = maxprojdivcurr(isbigger);
  ztmp(isbigger) = zcurr(isbigger)+v0(3)-1;
  maxprojdiv(v0(otherdims(1)):v1(otherdims(1)),v0(otherdims(2)):v1(otherdims(2))) = tmp;
  zdiv(v0(otherdims(1)):v1(otherdims(1)),v0(otherdims(2)):v1(otherdims(2))) = ztmp;
end

if isjob,
  save(outmatfile,'maxprojdiv','zdiv');
end
