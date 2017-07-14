function [precision,recall,hfig,hax] = ...
  ComputePrecisionRecall(labellocs,predlocs,scores,thresholds_try,timepoints_test,varargin)

[match_rad,isoscale,hfig,hax,doplot,plotcolor] = myparse(varargin,'match_rad',[50,2],...
  'isoscale',[1,1,5],'hfig',[],'hax',[],'doplot',true,'plotcolor','k');

if doplot && (numel(hax) < 2 || ~all(ishandle(hax(1:2)))),
  if isempty(hfig),
    hfig = figure;
  elseif ~ishandle(hfig),
    figure(hfig);
  else
    set(0,'CurrentFigure',hfig);
  end
  hax = createsubplots(2,1,.05,hfig);
end

npred = numel(scores);
nlabel = size(labellocs,1);

canmatch = true(npred,nlabel);
for i = 1:nlabel,
  d = sqrt( sum((predlocs(:,1:3) - labellocs(i,1:3)).^2 .* isoscale.^2,2) );
  canmatch(:,i) = d <= match_rad(1) & ...
    (abs(predlocs(:,4) - labellocs(i,4)) <= match_rad(2));
end

istesttimepoint = ismember(round(predlocs(:,4)),timepoints_test);

% to compute precision, determine fraction of detections that *can* match,
% not one-to-one

precision = nan(1,numel(thresholds_try));
recall = nan(1,numel(thresholds_try));

if doplot,
  for i = 1:2,
    hold(hax(i),'on');
    set(hax(i),'XScale','log');
  end
  ylabel(hax(1),'Precision');
  ylabel(hax(2),'Recall');
end

for threshi = 1:numel(thresholds_try),
  
  thresh = thresholds_try(threshi);
  idxthresh = scores >= thresh;
  
  recall(threshi) = nnz(any(canmatch(idxthresh,:),1)) / nlabel;
  precision(threshi) = nnz(any(canmatch(idxthresh&istesttimepoint,:),2)) / nnz(idxthresh&istesttimepoint);

  if doplot,
    
    plot(hax(1),thresh,precision(threshi),'.','Color',plotcolor);
    plot(hax(2),thresh,recall(threshi),'.','Color',plotcolor);
    drawnow;
  end
end
