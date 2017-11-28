%% Description

% Script that makes a video of the density of divisions at each
% spatiotemporal location. Note that this uses compiled versions of the
% MATLAB functions to run the code in parallel on the Janelia cluster. 

%% set up paths  

addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;
addpath /groups/branson/home/bransonk/codepacks/keller-lab-block-filetype/matlabWrapper/

% maxprojdatadir = '/media/KellerS8/SV1/14-05-21/DivisionDetection/Kristin/Projections';
maxprojdatadir = '/media/KellerS8/SV1/14-05-21/DivisionDetection/Kristin/ProjectionsWavelet';
%maxprojdatadir = '/media/KellerS8/SV1/14-05-21/DivisionDetection/Kristin/ProjectionsNew';

timestamp = '20170818';
dim = 3;
switch dim,
  case 1,
    projstr = 'yz';
  case 2,
    projstr = 'xz';
  case 3,
    projstr = 'xy';
end
    

maxprojfilestr = sprintf('*%sProjection*.klb',projstr);

% where to save division density results
% outdir = '/nrs/branson/MouseLineaging/ComputeDivisionDensity_xy_curated';
outdir = sprintf('/nrs/branson/MouseLineaging/ComputeDivisionDensity_%s_%s',projstr,timestamp);
% inmatfile = '/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/ComputeDivisionDensityInput_xy_curated.mat';
inmatfile = sprintf('/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/ComputeDivisionDensityInput_%s_%s.mat',projstr,timestamp);
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

% curateddatafile = 'Predictions Manually Curated.csv';
curateddatafile = 'Predictions Manually Curated.v2.csv';

% whether to save a video
dosavevideo = true;
%outvideofile = 'Divisions_xy_curated.avi';
outvideofile = sprintf('Divisions_%s_%s.avi',projstr,timestamp);

%% load in data

load(datafile,'nonmaxpredidx_combined','sz');

fid = fopen(curateddatafile);
s = fgetl(fid);
ss = regexp(strtrim(s),',','split');
xi = find(strcmpi(ss,'X'));
yi = find(strcmpi(ss,'Y'));
zi = find(strcmpi(ss,'Z'));
ti = find(strcmpi(ss,'T'));
scorei = find(strcmpi(ss,'Score'));
curateddata = zeros(0,5);
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  if isempty(s),
    continue;
  end
  ss = regexp(s,',','split');
  datacurr = str2double(ss([xi,yi,zi,ti,scorei]));
  curateddata(end+1,:) = datacurr;
end

fclose(fid);


%% parameters

% scales for x,y,z
isoscale = [1,1,5];

% parameters for computing division density
cdd = struct;
cdd.filsig = [12.5,12.5,2.5,4];
%cdd.filsig = [50,50,10,4];
cdd.filrad = ceil(2.5*cdd.filsig);
cdd.filrad(4) = 4;
cdd.thresh = 0.084;
cdd.ncores = 4;
cdd.list = curateddata;
%cdd.list = nonmaxpredidx_combined;
cdd.sz = sz;
cdd.stepsz = [440,440,200,9];
cdd.dim = dim;

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

T = min(max(round(cdd.list(:,4))),max(maxprojtimestamps));

%% compute division density locally

params = struct2paramscell(rmfield(cdd,{'list','sz'}));
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
  %out3dfile = fullfile(outdir,sprintf('%d.klb',t));
  out3dfile = '';
  if exist(outmatfile,'file'),
    %continue;
  end
  jobid = sprintf('CDD%d',t);
  fid = fopen(scriptfile,'w');
  fprintf(fid,'if [ -d %s ]\n',TMP_ROOT_DIR);
  fprintf(fid,'  then export MCR_CACHE_ROOT=%s.%s\n',MCR_CACHE_ROOT,jobid);
  fprintf(fid,'fi\n');
  fprintf(fid,'%s %s %s %d %s %s\n',...
    SCRIPT,MCR,inmatfile,t,outmatfile,out3dfile);
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
  vidobj = VideoWriter(outvideofile,'Uncompressed AVI');
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

%% make a file with metadata
% 
% fid = fopen('KLBFileInfo_xy_curated.csv','w');
% fprintf(fid,'Filename,MaxValue,T,X0,Y0,Z0,X1,Y1,Z1\n');
% allfiles = mydir(outdir,'name','.*_.*.klb');
% filets = regexp(allfiles,'/(\d+)_','once','tokens');
% filets = str2double([filets{:}]);
% assert(numel(allfiles)==numel(filets));
% 
% for t = 1:T,
%   
%   fprintf('t = %d\n',t);
%   
%   files = allfiles(filets==t);
% 
%   for i = 1:numel(files),
%   
%     header = readKLBheader(files{i});
%     ss = regexp(strtrim(header.metadata(1:255)),',','split');
%     
%     j = find(~cellfun(@isempty,regexp(ss,'^maxd','once')));
%     k = find(ss{j}==':',1);
%     maxd = str2double(ss{j}(k+1:end));
%     
%     j = find(~cellfun(@isempty,regexp(ss,'^v0','once')));
%     k = find(ss{j}==':',1);
%     v0 = eval(ss{j}(k+1:end));
%     
%     j = find(~cellfun(@isempty,regexp(ss,'^v1','once')));
%     k = find(ss{j}==':',1);
%     v1 = eval(ss{j}(k+1:end));
%     
%     [~,n,ext] = fileparts(files{i});
%     fprintf(fid,'%s%s,%e,%d,%d,%d,%d,%d,%d,%d\n',n,ext,maxd,t,v0(1),v0(2),v0(3),v1(1),v1(2),v1(3));
%     
%   
%   end
% end
% fclose(fid);