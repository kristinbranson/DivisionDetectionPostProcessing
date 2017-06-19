function [peaks_resize,peaks] = selectLocalMaxima(predim,varargin)

predsz = size(predim);
[thresh,npadxy,outsz,algorithm] = myparse(varargin,'thresh',150,...
  'npadxy',[16,16],'outsz',[2169,2048],...
  'algorithm','regionalmax');

switch lower(algorithm),
  
  case 'regionalmax',

    ispeak = imregionalmax(predim) & predim >= thresh;
    cc = bwconncomp(ispeak);
    peaks = regionprops(cc,'Centroid');
    peaks = cat(1,peaks.Centroid);
    
  otherwise,

    error('Algorithm %s not implemented',algorithm);
    
end
    
peaks_resize = (peaks+npadxy).*outsz(2:-1:1)./(predsz(2:-1:1)+2*npadxy);