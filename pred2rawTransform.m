function predim2 = pred2rawTransform(predim,varargin)

[interp,outsz,npadxy] = myparse(varargin,'interp','bilinear','outsz',[2169,2048],'npadxy',[16,16]);

[nr,nc,ncolors] = size(predim);
predclass = class(predim);

predim2 = imresize( cat(2,zeros([nr+2*npadxy(2),npadxy(1),ncolors],predclass),...
  cat(1,zeros([npadxy(2),nc,ncolors],predclass),predim,zeros([npadxy(2),nc,ncolors],predclass)),...
  zeros([nr+2*npadxy(2),npadxy(1),ncolors],predclass)), outsz, interp );