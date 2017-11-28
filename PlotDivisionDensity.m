% [hfig,hax,colorim] = PlotDivisionDensity(maxprojim,maxprojdiv,...)
%
% Uses division density to set hue and raw image to set intensity, and
% plots the results.
% 
% Inputs:
% maxprojim: Maximum intensity projection of raw image. 
% maxprojdiv: Maximum intensity projection of division density. 
% 
% Output: 
% hfig: Handle to figure plotted on
% hax: Handle to axes plotted on
% colorim: Resulting image from combining raw and division density images. 
% 
% Optional inputs:
% maxrho: Maximum value for scaling division densities. Default =
% max(maxprojdiv(:))
% maxint: Maximum value for scaling raw image. Default = max(maxprojim(:))
% hax, hfig: Handles of axes, figure to plot on. By default, new figure and
% axes will be created. 

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