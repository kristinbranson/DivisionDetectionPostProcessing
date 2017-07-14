%% set up paths  

addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;
addpath /groups/branson/home/bransonk/codepacks/keller-lab-block-filetype/matlabWrapper/

maxprojdatadir = '/media/KellerS8/SV1/14-05-21/DivisionDetection/Kristin/Projections';
maxprojfilestr = '*yzProjection*.klb';

% where to save division density results
outdir = '/nrs/branson/MouseLineaging/ComputeDivisionDensity_yz_v2';
inmatfile = '/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/ComputeDivisionDensityInput_yz_v2.mat';
if ~exist(outdir,'dir'),
  mkdir(outdir);
end

% stuff for running on cluster

% location of compiled division density code
SCRIPT = '/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/ComputeDivisionDensity/for_redistribution_files_only/run_ComputeDivisionDensity.sh';

% location of MCR
MCR = '/groups/branson/bransonlab/share/MCR/v91';

% where to save temp data on the cluster
username = GetUserName();
TMP_ROOT_DIR = fullfile('/scratch',username);
MCR_CACHE_ROOT = fullfile(TMP_ROOT_DIR,'mcr_cache_root');

% location of division data
datafile = 'DivisionDetections20170714.mat';

% whether to save a video
dosavevideo = true;
outvideofile = 'Divisions_yz.avi';

%% load in data

load(datafile,'nonmaxpredidx_combined','sz');

%% parameters

% scales for x,y,z
isoscale = [1,1,5];

% parameters for computing division density
cdd = struct;
cdd.filsig = [12.5,12.5,2.5,4];
cdd.filrad = ceil(2.5*cdd.filsig);
cdd.filrad(4) = 4;
cdd.thresh = 0.084;
cdd.ncores = 4;
cdd.list = nonmaxpredidx_combined;
cdd.sz = sz;
cdd.stepsz = [220,220,100,9];
cdd.dim = 1;

% plotting parameters

% to compute the maximum intensity for scaling image intensity
prctile_maxint = 99.9;
intscale = .75;
% smooth the fit of maximum intensity over time
filsig_maxint = 20;
filrad_maxint = 50;

% output video
outfps = 10;

%% locations of MIPs

maxprojfiles = mydir(fullfile(maxprojdatadir,maxprojfilestr));
m = regexp(maxprojfiles,'_TM(\d+)_','once','tokens');
maxprojtimestamps = str2double([m{:}]);
[maxprojtimestamps,order] = sort(maxprojtimestamps);
maxprojfiles = maxprojfiles(order);

%% compute division density locally

params = struct2paramscell(cdd);
for t = 1:T,
  res = struct;
  [res.maxprojdiv,res.zdiv] = ComputeDivisionDensity(cdd.list,t,cdd.sz,params{:});
  outmatfile = fullfile(outdir,sprintf('%d.mat',t));
  save(outmatfile,'-struct','res');
end

%% compute division density on the cluster

save(inmatfile,'-struct','cdd');

for t = 1:T,
  
  scriptfile = fullfile(outdir,sprintf('%d.sh',t));
  logfile = fullfile(outdir,sprintf('%d.log',t));
  outmatfile = fullfile(outdir,sprintf('%d.mat',t));
  if exist(outmatfile,'file'),
    continue;
  end
  jobid = sprintf('CDD%d',t);
  fid = fopen(scriptfile,'w');
  fprintf(fid,'if [ -d %s ]\n',TMP_ROOT_DIR);
  fprintf(fid,'  then export MCR_CACHE_ROOT=%s.%s\n',MCR_CACHE_ROOT,jobid);
  fprintf(fid,'fi\n');
  fprintf(fid,'%s %s %s %d %s\n',...
    SCRIPT,MCR,inmatfile,t,outmatfile);
  fclose(fid);
      
  unix(sprintf('chmod u+x %s',scriptfile));
  cmd = sprintf('ssh login1 ''source /etc/profile; bsub -n %d -J %s -o ''%s'' ''\"%s\"''''',...
    cdd.ncores,jobid,logfile,scriptfile);
  unix(cmd);
  
end

%% choose maximum values for making color image

maxrhos = nan(1,T);
maxints = nan(1,T);
for t = 1:T,
  if isnan(maxints(t)),
    i = find(maxprojtimestamps == t);
    maxprojim = readKLBstack(maxprojfiles{i});
    maxints(t) = prctile(maxprojim(maxprojim>0),prctile_maxint);
  end
  if isnan(maxrhos(t)),
    outmatfile = fullfile(outdir,sprintf('%d.mat',t));
    if ~exist(outmatfile,'file'),
      continue;
    end
    tmp = load(outmatfile);
    maxrhos(t) = max(tmp.maxprojdiv(:));
  end
end

fil = normpdf(-filrad_maxint:filrad_maxint,0,filsig_maxint);
maxintsmooth = imfilter(maxints,fil,'same','symmetric');
maxrho = nanmedian(maxrhos);

%% save a video

if dosavevideo,
  vidobj = VideoWriter(outvideofile);
  vidobj.FrameRate = outfps;
  open(vidobj);
end

otherdims = setdiff(1:3,cdd.dim);

hfig = 1;
figure(hfig);
clf;
hax = gca;
for t = 1:T,
  
  outmatfile = fullfile(outdir,sprintf('%d.mat',t));
  if ~exist(outmatfile,'file'),
    break;
  end
  tmp = load(outmatfile);
  i = find(maxprojtimestamps == t);
  maxprojim = readKLBstack(maxprojfiles{i});

  [~,~,colorim] = PlotDivisionDensity(maxprojim,tmp.maxprojdiv,'hax',hax,'maxint',maxintsmooth(t)*intscale,'maxrho',maxrho);
  title(num2str(t));
  drawnow;
  
  if dosavevideo,
    % make sure the size is a factor of 4
    if t == 1,
      newsz = round(size(colorim).*[isoscale(otherdims),1]/4)*4;
      newsz = newsz(1:2);
    end
    colorim = max(0,min(1,imresize(colorim,newsz)));
    writeVideo(vidobj,colorim);
  end

end

if dosavevideo,
  close(vidobj);
end