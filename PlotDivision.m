function [hfig,hax] = PlotDivision(allinputdatafiles,alltimestamp,allpredidx,loc,rawsz,varargin)

[scores_pert,peaks_pert,timepoints_pa,labellocs,thresh,boxrad,...
  hfig,textstr0,textstr1,filename,figpos] = ...
  myparse(varargin,...
  'scores_pert',{},'peaks_pert',{},'timepoints_pa',[],...
  'labellocs',zeros(0,4),...
  'thresh',.001,...
  'boxrad',[100,100,5,3],...
  'hfig',[],...
  'textstr0','',...
  'textstr1','',...
  'filename','',...
  'figpos',[10,10,2400,760]);

ncolors = 256;
w = linspace(0,1,ncolors)';
colors = w.*[1,0,1] + (1-w).*[1,1,1];
ntimepoints = numel(allinputdatafiles);
sz = [rawsz(1:3),max(allpredidx(:,4))];
  
v0 = double(max(1,floor(loc)-boxrad));
v1 = double(min([rawsz(1:3),ntimepoints],ceil(loc)+boxrad));
maint = loc(end);
    
ncurr = v1 - v0 + 1;

rawvid = zeros(ncurr([1,2,4]),'uint16');
for t = v0(end):v1(end),
  i = find(alltimestamp == t,1);
  rawvid(:,:,t-v0(end)+1) = max(readKLBroi(allinputdatafiles{i},[v0(1:3),1,1;v1(1:3),1,1]),[],3);
end

predvid = zeros(ncurr([1,2,4]));
for t = v0(end):v1(end),
  predvid(:,:,t-v0(end)+1) = ConstructPredMaxProj2(allpredidx,sz,t,...
    'xlim',[v0(1),v1(1)],'ylim',[v0(2),v1(2)],'zlim',[v0(3),v1(3)]);
end

if isempty(hfig),
  hfig = figure;
elseif ~ishandle(hfig),
  figure(hfig);
else
  set(0,'CurrentFigure',hfig);
end
  
clf(hfig);
set(hfig,'Units','pixels','Position',figpos);
hax = createsubplots(2,ncurr(end),.01,hfig);
hax = reshape(hax,[2,ncurr(end)]);
for ti = 1:ncurr(end),
  t = ti + v0(end)-1;
  
  % predictions in the current timepoint
  tipred = find(timepoints_pa==t);
  if isempty(tipred),
    detis = [];
  else
    scores = scores_pert{tipred};
    peaks = peaks_pert{tipred};
    detis = find(scores >= thresh & ...
      peaks(:,4) == t & ...
      peaks(:,1) >= v0(1) & peaks(:,1) <= v1(1) & ...
      peaks(:,2) >= v0(2) & peaks(:,2) <= v1(2) & ...
      peaks(:,3) >= v0(3) & peaks(:,3) <= v1(3));
  end
  
  % labels in the current timepoint
  labelis = find(labellocs(:,4) == t & ...
    labellocs(:,1) >= v0(1) & labellocs(:,1) <= v1(1) & ...
    labellocs(:,2) >= v0(2) & labellocs(:,2) <= v1(2) & ...
    labellocs(:,3) >= v0(3) & labellocs(:,3) <= v1(3));

  % draw raw image
  imagesc([v0(2),v1(2)],[v0(1),v1(1)],rawvid(:,:,ti),'Parent',hax(1,ti));
  
  % add label text
  if ti == 1,
    s = sprintf(' z in [%d,%d], t = %d',v0(3),v1(3),ti+v0(end)-1);
    if ~isempty(textstr0),
      s = [textstr0,',',s];
    end
    text(v0(2),v0(1),s,...
      'Parent',hax(1,1),'Color','w',...
      'HorizontalAlignment','left','VerticalAlignment','top');
  else
    s = sprintf(' t = %d',ti+v0(end)-1);
    if ~isempty(textstr1),
      s = [textstr1,',',s];
    end
    text(v0(2),v0(1),s,'Parent',hax(1,ti),'Color','w',...
      'HorizontalAlignment','left','VerticalAlignment','top');
  end
  axis(hax(1,ti),'image');
  hold(hax(1,ti),'on');

  % draw raw predictions
  imagesc([v0(2),v1(2)],[v0(1),v1(1)],predvid(:,:,ti),'Parent',hax(2,ti));
  axis(hax(2,ti),'image');
  hold(hax(2,ti),'on');

  % plot all predictions
  for ii = 1:numel(detis),
    i = detis(ii);
    colori = ceil(scores(i)/max(scores)*ncolors);
    for j = 1:2,
      plot(hax(j,ti),double(peaks(i,2)),double(peaks(i,1)),'.','Color',colors(colori,:));
      if j == 2,
        text(double(peaks(i,2)),double(peaks(i,1)),sprintf('%.2f',scores(i)),...
          'HorizontalAlignment','left','VerticalAlignment','middle',...
          'Color',colors(colori,:),'Parent',hax(j,ti));
      end
    end
  end

  % plot all labels
  for ii = 1:numel(labelis),
    i = labelis(ii);
    for j = 1:2,
      plot(hax(j,ti),labellocs(i,2),labellocs(i,1),'o','Color',[1,.6,0]);
    end
  end
  
end
    
colormap parula;
set(hax(1,:),'XTickLabel',{});
set(hax(:,2:end),'YTickLabel',{});
set(hax,'Box','off');
set(hax(1,:),'CLim',[0,max(rawvid(:))]);
set(hax(2,:),'CLim',[0,1]);
%impixelinfo;
drawnow;

if ~isempty(filename),
  set(hfig,'InvertHardCopy','off','Color','w');
  savefig_pa(filename,hfig,'png');
end