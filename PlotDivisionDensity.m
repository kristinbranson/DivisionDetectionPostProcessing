function [hfig,hax,colorim] = PlotDivisionDensity(maxprojim,maxprojdiv,varargin)

[maxrho,maxint,hax,hfig] = myparse(varargin,'maxrho',[],'maxint',[],'hax',nan,'hfig',[]);
if isempty(maxrho),
  maxrho = max(maxprojdiv(:));
end

ncolors = 256;
cm = jet(ncolors+95);
cm = cm(76:end-20,:);
%cm = autumn(ncolors);
color_div = colormap_image(maxprojdiv,cm,[0,maxrho]);

maxprojim = double(maxprojim);
if isempty(maxint),
  maxint = max(maxprojim(:));
end
maxprojim = min(1,maxprojim ./ maxint);
colorim = maxprojim .* color_div;

if ~ishandle(hax),
  if isempty(hfig),
    hfig = figure;
  elseif ~ishandle(hfig),
    figure(hfig);
  else
    set(0,'CurrentFigure',hfig);
  end
  hax = gca;
end
image(colorim,'Parent',hax);