%% paths
if ispc,

  addpath C:/Code/Jdetect/misc;
  addpath C:/Code/Jdetect/filehandling;
  addpath C:/Code/KLB;
  datadir = 'P:\SV1\14-05-21\DivisionDetection\Predictions_0427';
  
else

  addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;
  addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;
  addpath /groups/branson/home/bransonk/codepacks/keller-lab-block-filetype/matlabWrapper/
  %datadir = '/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/DivisionDetectionData_14-05-21/Predictions_0427';
  %rawdatadir = '/nrs/turaga/bergera/division_detection/full_h5_timeseries';
  rawdatadir = '/media/KellerS8/SV1/14-05-21/Mmu_E1_CAGTAG1.corrected/Results/TimeFused.Corrected';
  %preddatafile = '/nrs/turaga/bergera/division_detection/prediction_outbox/sharpen_only_no_augment_small.h5';
  preddatadir = '/nrs/turaga/bergera/division_detection/prediction_outbox/mk4_large_balanced_augmented_bn/sparse';
  
  preddatafile = '/nrs/turaga/bergera/division_detection/prediction_outbox/lr_large_balanced_bn.h5';
  %gtdatafiles = {'/groups/turaga/home/bergera/data/div_detect/annotations/partial/divisionAnnotations.mat'};
  gtdatafiles = {'/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/divisionAnnotations.mat'};
end

% sudo mount //Keller-S8/Processing -t cifs -o uid=990313,gid=93319,username=SiMView /media/KellerS8/

%% parameters

% datafiletype = 'klb';
% npadxyzt = [16,16,37,0];
predthresh = .6;

predfiletype = 'sparseh5';
preddatasetname = '/coo';
predshapename = '/shape';

rawfiletype = 'klb';
% predfiletype = 'h5';
% preddatasetname = '/predictions';
%npadxyzt = [8,8,8,3];
npadxyzt = [0,0,0,0];
%rawfilestr = '*00.klb';
rawisregexp = true;
rawisrecursive = true;
rawmaxdepth = 1;
rawfilestr = 'SPM00_TM(\d+)_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';
poslabels = [1,2,3,4,5,103];


%% find data files
switch rawfiletype,
  case 'klb'
    if rawisregexp,
      inputdatafiles = mydir(rawdatadir,'name',rawfilestr,'recursive',rawisrecursive,'maxdepth',rawmaxdepth);
    else
      inputdatafiles = mydir(fullfile(rawdatadir,rawfilestr),'recursive',rawisrecursive,'maxdepth',rawmaxdepth);
    end
  case 'h5',
    inputdatafiles = mydir(fullfile(rawdatadir,'*h5'));    
  otherwise,
    error('Not implemented');
end
timestamp = regexp(inputdatafiles,'TM(\d+)_','once','tokens');
timestamp = str2double([timestamp{:}]);
[timestamp,order] = sort(timestamp);
inputdatafiles = inputdatafiles(order);
assert(all(timestamp == (0:numel(timestamp)-1)));

switch predfiletype,
  case 'klb',
    preddatafiles = cell(size(inputdatafiles));
    ismissing = false(1,numel(inputdatafiles));
    for i = 1:numel(inputdatafiles),
      preddatafile = [inputdatafiles{i}(1:end-5),'1.klb'];
      if ~exist(preddatafile,'file'),
        ismissing(i) = true;
        [~,n] = myfileparts(inputdatafiles{i});
        fprintf('Could not find predictions for %s\n',n);
      else
        preddatafiles{i} = preddatafile;
      end
    end
    
    inputdatafiles(ismissing) = [];
    preddatafiles(ismissing) = [];
  case 'h5',
  case 'sparseh5',
    preddatafiles = mydir(fullfile(preddatadir,'*.h5'));
    m = regexp(preddatafiles,'/(\d+)\.h5','tokens','once');
    m = [m{:}];
    predtimestamps = str2double(m);
    [predtimestamps,order] = sort(predtimestamps);
    preddatafiles = preddatafiles(order);
    [ism,idx] = ismember(predtimestamps,timestamp);
    assert(all(ism));
    inputdatafiles = inputdatafiles(idx);
    timestamp = timestamp(idx);
    
end

labellocs = nan(0,4);

for i = 1:numel(gtdatafiles),
  
  gtdatafile = gtdatafiles{i};
  [~,~,ext] = fileparts(gtdatafile);
  
  switch ext,
    case '.csv',
      fid = fopen(gtdatafile,'r');
      s = fgetl(fid);
      ss = regexp(s,',','split');
      
      xidx = find(strcmpi(ss,'X'));
      yidx = find(strcmpi(ss,'Y'));
      zidx = find(strcmpi(ss,'Z'));
      tidx = find(ismember(ss,{'T','Time point'}));
      lidx = find(strcmpi(ss,'Annotation label'));
      
      while true,
        s = fgetl(fid);
        if ~ischar(s),
          break;
        end
        s = strtrim(s);
        ss = regexp(s,',','split');
        n = str2double(ss);
        if ~isempty(lidx) && ~ismember(n(lidx),poslabels),
          continue;
        end
        labellocs(end+1,:) = n([xidx,yidx,zidx,tidx]);
      end
      fclose(fid);
    case '.mat'
      ll = load(gtdatafile);
      ispos = ismember(ll.divisionAnnotations(:,2),[1,2,3,4,5,103]);
      labellocs(end+1:end+nnz(ispos),:) = ll.divisionAnnotations(ispos,[3,4,5,1]);
  end
end

%% load one

filenum = 6;
[rawreadframe,rawnframes,rawfid,rawheaderinfo] = get_readframe_fcn(inputdatafiles{filenum});

switch predfiletype,
  case 'klb',    
    [predreadframe,prednframes,predfid,predheaderinfo] = get_readframe_fcn(preddatafiles{filenum});
  case 'h5',
    [predreadframe,prednframes,predfid,predheaderinfo] = ...
      get_readframe_fcn(preddatafile,'datasetname',preddatasetname,'slicedims',3,'fixdims',4,'fixidx',filenum);
  case 'sparseh5'
    preddata = h5read(preddatafiles{filenum},preddatasetname);
    preddatashape = h5read(preddatafiles{filenum},predshapename);
    predreadframe = @(z) ConstructPredIm(z,preddata,preddatashape');
  otherwise
    error('Not implemented');
end

im = rawreadframe(1);
rawsz = [rawheaderinfo.xyzct(1:3),numel(inputdatafiles)];
rawclass = class(im);

im = predreadframe(1);
predsz = size(im);
predclass = class(im);


%% look at max projections
maxv = zeros(rawsz(1:2));
maxi = zeros(rawsz(1:2));
for i = 1:rawnframes,
  im = rawreadframe(i);
  idx = im > maxv;
  maxv(idx) = im(idx); 
  maxi(idx) = i;
  fprintf('i = %d/%d\n',i,rawnframes);
end

predmaxv = zeros(predsz);
predmaxi = zeros(predsz);
for i = 1:prednframes,
  im = predreadframe(i);
  idx = im > predmaxv;
  predmaxv(idx) = im(idx); 
  predmaxi(idx) = i;
  fprintf('i = %d/%d\n',i,prednframes);
end
predmaxi = predmaxi + npadxyzt(3);

miniall = min(prctile(maxi(maxi>0),1),prctile(predmaxi(predmaxv>predthresh),1));
maxiall = max(prctile(maxi(maxi>0),99),prctile(predmaxi(predmaxv>predthresh),99));

cmdepth = jet(maxiall-miniall+2);
preddepthim = colormap_image(predmaxi,cmdepth,[miniall,maxiall]).*predmaxv/max(predmaxv(:));
preddepthim_resize = pred2rawTransform(preddepthim,'outsz',rawsz(1:2),'npadxy',npadxyzt(1:2),'interp','bilinear');

rawdepthim = colormap_image(maxi,cmdepth,[miniall,maxiall]).*maxv/max(maxv(:));

hfig = 1;
figure(hfig);
clf;
hax = [];
hax(1) = subplot(1,2,1);
image(rawdepthim);
axis image;

hax(2) = subplot(1,2,2);
image(preddepthim_resize);
axis image;
impixelinfo;
linkaxes(hax);
colormap(cmdepth);
set(gca,'CLim',[miniall,maxiall]);
hcb = colorbar('Location','East');
set(hcb,'Color','w');

%% grab parts of the image stack according to max depth for prediction map
% predmaxi_resize = pred2rawTransform(predmaxi,'outsz',rawsz,'npadxy',npadxyzt(1:2),'interp','nearest');
% zs = unique(predmaxi_resize);
% rawnearpredmaxi = zeros(rawsz);
% for z = zs(:)',
%    im = rawreadframe(z);
%    idx = predmaxi_resize == z;
%    rawnearpredmaxi(idx) = im(idx);
% end
% 
% rawnearpredmaxi_depthim = colormap_image(predmaxi_resize,cmdepth,[miniall,maxiall]).*rawnearpredmaxi/max(rawnearpredmaxi(:));
% 
% hfig = 2;
% figure(hfig);
% clf;
% hax = [];
% hax(1) = subplot(1,2,1);
% image(rawnearpredmaxi_depthim);
% axis image;
% 
% hax(2) = subplot(1,2,2);
% image(preddepthim_resize);
% axis image;
% impixelinfo;
% linkaxes(hax);
% colormap(cmdepth);
% set(gca,'CLim',[miniall,maxiall]);
% hcb = colorbar('Location','East');
% set(hcb,'Color','w');
% 
% % this is pretty useless!!

%% let's try a detection algorithm

ispeak = imregionalmax(predmaxv) & predmaxv >= predthresh;
cc = bwconncomp(ispeak);
peaks = regionprops(cc,'Centroid');
peaks = cat(1,peaks.Centroid);
peaks_resize = (peaks+npadxyzt(1:2)).*rawsz(2:-1:1)./(predsz(2:-1:1)+2*npadxyzt(1:2));
peakmaxv = nan(1,cc.NumObjects);
for i = 1:cc.NumObjects,
  peakmaxv(i) = max(predmaxv(cc.PixelIdxList{i}));
end

hfig = 3;
figure(hfig);
clf;
hax = [];
hax(1) = subplot(1,2,1);
image(rawdepthim);
axis image;
hold on;
myscatter(peaks_resize(:,1),peaks_resize(:,2),[],peakmaxv,16,repmat(linspace(.25,1,16)',[1,3]),'+');
%plot(peaks_resize(:,1),peaks_resize(:,2),'w+');

hax(2) = subplot(1,2,2);
image(preddepthim_resize);
axis image;
hold on;
myscatter(peaks_resize(:,1),peaks_resize(:,2),[],peakmaxv,16,repmat(linspace(.25,1,16)',[1,3]),'+');
%plot(peaks_resize(:,1),peaks_resize(:,2),'w+');

impixelinfo;
linkaxes(hax);
colormap(cmdepth);
set(gca,'CLim',[miniall,maxiall]);
hcb = colorbar('Location','East');
set(hcb,'Color','w');

%% read in info about each time point

switch predfiletype,
  
  case 'klb',
    ntimepoints = numel(inputdatafiles);
    fprintf('Reading klb file headers\n');
    for i = 1:ntimepoints,
      fprintf('Timepoint %d / %d\n',i,ntimepoints);
      predheaderinfo = readKLBheader(preddatafiles{i});
      if i == 1,
        predsz = predheaderinfo.xyzct(1:3);
      else
        assert(all(predsz == predheaderinfo.xyzct(1:3)));
      end
    end
  case 'h5'
    ntimepoints = predheaderinfo.Dataspace.Size(4);
    predsz = predheaderinfo.Dataspace.Size(1:3);
end
nel = prod(predsz)*ntimepoints;

%% threshold at predthresh while reading into a sparse matrix

nxy = prod(predsz(1:2));
nadd = zeros(ntimepoints,predsz(3));
allidx = [];
allv = [];
t0 = 1;
for t = 1:ntimepoints-t0+1,
  fprintf('t = %d / %d\n',t,ntimepoints);
  switch predfiletype,
    case 'klb',
      [predreadframe] = get_readframe_fcn(preddatafiles{t+t0-1});
    case 'h5',
      [predreadframe] = get_readframe_fcn(preddatafile,'datasetname',preddatasetname,'slicedims',3,'fixdims',4,'fixidx',t+t0-1);
  end
  for z = 1:predsz(3),
    if mod(z,100) == 0,
      fprintf('t = %d/%d, z = %d/%d\n',t,ntimepoints,z,predsz(3));
    end
    idxcurr = nxy*(z-1+predsz(3)*(t-1)) + (1:nxy);
    imcurr = predreadframe(z);
    isthresh = imcurr >= predthresh;
    if ~any(isthresh(:)),
      continue;
    end
    ncurr = nnz(isthresh);
    allidx(end+1:end+ncurr) = idxcurr(isthresh);
    allv(end+1:end+ncurr) = imcurr(isthresh);
    nadd(t,z) = nnz(isthresh);
  end
  fprintf('Total size is %d elements\n',numel(allv));
  save tmp2.mat allidx allv;
end
imthresh = sparse(allidx,ones(size(allidx)),allv,nel,1);

% nonmaximal supporession
chunksize = [50,50,50,20];
chunkstep = chunksize/2;
[peaks_resize,scores,peaks] = selectLocalMaximaND(imthresh,[predsz,ntimepoints],...
  'thresh',predthresh,'chunksize',chunksize,'chunkstep',chunkstep,'outsz',rawsz(1:2),...
  'npadxyz',npadxyzt);

%% plot a detection

[sortedscores,order] = sort(scores,2,'descend');

for detorderi = 1:numel(scores),

  deti = order(detorderi);
  
  loc = peaks_resize(deti,:);
  %boxrad = [100,100,0,20];
  boxrad = [100,100,5,3];
  
  v0 = max(1,floor(loc)-boxrad);
  v1 = min([rawsz(1:3),ntimepoints],ceil(loc)+boxrad);
  
  rxy = rawsz(2:-1:1)./(predsz(2:-1:1)+2*npadxyzt(1:2));
  
  predv0 = v0;
  predv0(1:2) = v0(1:2) ./ rxy - npadxyzt(1:2);
  predv0(3) = predv0(3) - npadxyzt(3);
  predv0 = max(1,floor(predv0));
  
  predv1 = v1;
  predv1(1:2) = v1(1:2) ./ rxy - npadxyzt(1:2);
  predv1(3) = predv1(3) - npadxyzt(3);
  predv1 = min(ceil(predv1),[predsz(1:3),ntimepoints]);
  
  ncurr = v1 - v0 + 1;
  predncurr = predv1 - predv0 + 1;
  rawvid = zeros(ncurr,'uint16');
  for t = v0(end):v1(end),
    rawvid(:,:,:,t-v0(end)+1) = readKLBroi(inputdatafiles{t+t0-1},[v0(1:3);v1(1:3)]);
  end
  
  switch predfiletype,
    case 'klb'
      predvid = zeros(predncurr,'uint8');
      for t = v0(end):v1(end),
        predvid(:,:,:,t) = readKLBroi(preddatafiles{t+t0-1},[predv0(1:3);predv1(1:3)]);
      end
    case 'h5'
      predvid = h5read(preddatafile,preddatasetname,predv0+[0,0,0,t0-1],predv1-predv0+1);
  end
  
  hfig = 10;
  figure(hfig);
  clf;
  hax = createsubplots(2,ncurr(end),.01);
  hax = reshape(hax,[2,ncurr(end)]);
  for ti = 1:ncurr(end),
    imagesc([v0(2),v1(2)],[v0(1),v1(1)],max(rawvid(:,:,:,ti),[],3),'Parent',hax(1,ti));
    if ti == 1,
      text(v0(2),v0(1),sprintf(' Det rank %d, score = %.1f, z in [%d,%d], t = %d',...
        detorderi,scores(deti),v0(3),v1(3),ti+v0(end)-1),'Parent',hax(1,1),'Color','w',...
        'HorizontalAlignment','left','VerticalAlignment','top');
    else
      text(v0(2),v0(1),sprintf(' t = %d',ti+v0(end)-1),'Parent',hax(1,ti),'Color','w',...
        'HorizontalAlignment','left','VerticalAlignment','top');
    end
    axis(hax(1,ti),'image');
    
    imagesc([predv0(2),predv1(2)],[predv0(1),predv1(1)],max(predvid(:,:,:,ti),[],3),'Parent',hax(2,ti));
    hold(hax(1,ti),'on');
    h = nan(1,2);
    h(1) = plot(hax(1,ti),peaks_resize(deti,2),peaks_resize(deti,1),'r.');%'r+','MarkerSize',10,'LineWidth',2);
    hold(hax(2,ti),'on');
    h(2) = plot(hax(2,ti),peaks(deti,2),peaks(deti,1),'r.');%,'MarkerSize',10,'LineWidth',2);
    if abs(ti+v0(end)-1 - peaks(deti,end)) < 1,
      set(h,'Marker','+');
    end
    axis(hax(2,ti),'image');
  end
  colormap parula;
  set(hax,'XTickLabel',{},'YTickLabel',{});
  set(hax,'Box','off');
  %impixelinfo;
  drawnow;
  
  filename = fullfile('Detections20170523',sprintf('Rank%ddet%d_score%.1f_t%dx%dy%dz%d.png',...
    detorderi,deti,scores(deti),round(loc(3))+t0-1,round(loc(1)),round(loc(2)),round(loc(3))));
  set(hfig,'InvertHardCopy','off','Color','w');
  savefig_pa(filename,hfig,'png');
  
end

%% plot a label

labeli = 2;

loc = labellocs(labeli,:);
boxrad = [200,200,5,3];
v0 = max(1,floor(loc)-boxrad);
v1 = min(rawsz,ceil(loc)+boxrad);

ncurr = v1 - v0 + 1;
rawvid = zeros(ncurr,'uint16');
for t = v0(end):v1(end),  
  rawvid(:,:,:,t-v0(end)+1) = readKLBroi(inputdatafiles{t},[v0(1:3);v1(1:3)]);
end

hfig = 11;
figure(hfig);
clf;
hax = createsubplots(1,ncurr(end),.025);
for ti = 1:ncurr(end),
  imagesc([v0(2),v1(2)],[v0(1),v1(1)],max(rawvid(:,:,:,ti),[],3),'Parent',hax(ti));
  if ti == 1,
    text(v0(2),v0(1),sprintf('Label %d, z in [%d,%d], t = %d',...
      labeli,v0(3),v1(3),ti+v0(end)-1),'Parent',hax(1),'Color','r',...
      'HorizontalAlignment','left','VerticalAlignment','top');
  else
    text(v0(2),v0(1),sprintf('t = %d',ti+v0(end)-1),'Parent',hax(ti),'Color','r',...
      'HorizontalAlignment','left','VerticalAlignment','top');    
  end
  axis(hax(ti),'image');
  
  hold(hax(ti),'on');
  h = plot(hax(ti),labellocs(labeli,2),labellocs(labeli,1),'cx');%,'MarkerSize',10,'LineWidth',2);
  if abs(ti+v0(end)-1 - labellocs(labeli,end)) < 1,
    set(h,'Marker','o','Color','r');
  end
end
colormap gray;
set(hax(2:end),'YTickLabel',[]);

set(hax,'Box','off');
impixelinfo;