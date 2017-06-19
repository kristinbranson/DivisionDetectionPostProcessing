function [precision,recall,ntruepos,nfalsepos,nfalseneg,hfig,hax] = ...
  MatchLabelsAndPredictions(labellocs,predlocs,scores,thresholds_try,varargin)

[match_rad,hfig,hax,doplot,plotcolor] = myparse(varargin,'match_rad',[20,20,10,2],...
  'hfig',[],'hax',[],'doplot',true,'plotcolor','k');

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
  for d = 1:4,
    canmatch(:,i) = canmatch(:,i) & ...
      (predlocs(:,d) >= labellocs(i,d)-match_rad(d)) & ...
      (predlocs(:,d) <= labellocs(i,d)+match_rad(d));
  end
  
end
costmatch = inf(size(canmatch));
costmatch(canmatch) = 0;

% no matter what the threshold, these have to be false positives and false
% negatives
isfalsepos0 = sum(canmatch,2)==0;
isfalseneg0 = sum(canmatch,1)==0;

% we will then do matching only on the remaining predictions and labels
nlabel0 = nnz(~isfalseneg0);

nfalsepos = nan(1,numel(thresholds_try));

nfalseneg = nan(1,numel(thresholds_try));


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
  costmatchcurr = costmatch(idxthresh&~isfalsepos0,~isfalseneg0);
  npred0 = size(costmatchcurr,1);
  if npred0 == 0,
    
    isfalseneg1 = true(1,nlabel0);
    isfalsepos1 = [];
    
  elseif nlabel0 == 0,
    
    isfalseneg1 = [];
    isfalsepos1 = true(1,npred0);
    
  else

    A = cat(2,cat(1,costmatchcurr,ones(nlabel0)),...
      cat(1,ones(npred0),zeros(nlabel0,npred0)));
    [label2pred,~] = hungarian(A);
    isfalseneg1 = label2pred(1:nlabel0) > npred0;
    isfalsepos1 = label2pred(nlabel0+1:end) <= npred0;
    
  end
    
  nfalsepos(threshi) = nnz(isfalsepos0&idxthresh) + nnz(isfalsepos1);
  nfalseneg(threshi) = nnz(isfalseneg0) + nnz(isfalseneg1);

  if doplot,
    
    ntrueposcurr = nlabel - nfalseneg(threshi);
    precisioncurr = ntrueposcurr ./ (ntrueposcurr + nfalsepos(threshi));
    recallcurr = ntrueposcurr ./ (ntrueposcurr + nfalseneg(threshi));

    plot(hax(1),thresh,precisioncurr,'.','Color',plotcolor);
    plot(hax(2),thresh,recallcurr,'.','Color',plotcolor);
    drawnow;
  end
end

ntruepos = nlabel - nfalseneg;

precision = ntruepos ./ (ntruepos + nfalsepos);
recall = ntruepos ./ (ntruepos + nfalseneg);
