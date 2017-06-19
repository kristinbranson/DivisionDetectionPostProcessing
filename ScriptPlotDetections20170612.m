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
  %preddatadir = '/nrs/turaga/bergera/division_detection/prediction_outbox/mk4_large_balanced_augmented_bn/sparse';
  preddatadir = '/nrs/turaga/bergera/division_detection/prediction_outbox/mk5_large_small_batch_balanced_good_reweight/sparse';
  %preddatafile = '/nrs/turaga/bergera/division_detection/prediction_outbox/lr_large_balanced_bn.h5';
  %gtdatafiles = {'/groups/turaga/home/bergera/data/div_detect/annotations/partial/divisionAnnotations.mat'};
  gtdatafiles = {'/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/divisionAnnotations.mat'
    'AnnotatedTimePoints.csv'
    'TPs_DivAnnotations.csv'};
  % some of the z-coordinates are in the isotropic coordinate system
  gtratios = [1,1,1,1
    1,1,.2,1
    1,1,.2,1];
  savedir = 'Detections20170617';
end

tgmmfile = 'tgmmDivisions.mat';

% sudo mount //Keller-S8/Processing -t cifs -o uid=990313,gid=93319,username=SiMView /media/KellerS8/

%% parameters

% threshold raw scores here for efficiency
predthresh = .01;

predfiletype = 'sparseh5';
preddatasetname = '/coo';
predshapename = '/shape';

rawfiletype = 'klb';
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
alltimestamp = timestamp;
inputdatafiles = inputdatafiles(order);
allinputdatafiles = inputdatafiles;
assert(all(timestamp == (0:numel(timestamp)-1)));

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
    
labellocs = nan(0,4);
datafilei = nan(0,1);

for i = 1:numel(gtdatafiles),
  
  gtdatafile = gtdatafiles{i};
  [~,~,ext] = fileparts(gtdatafile);
  
  switch ext,
    case '.csv',
      fid = fopen(gtdatafile,'r');
      s = fgetl(fid);
      ss = regexp(s,',','split');
      
      xidx = find(regexpcmp(ss,'X'));
      yidx = find(regexpcmp(ss,'Y'));
      zidx = find(regexpcmp(ss,'Z'));
      tidx = find(regexpcmp(ss,'T$'));
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
        labellocs(end+1,:) = n([xidx,yidx,zidx,tidx]).*gtratios(i,:);
        datafilei(end+1) = i;
      end
      fclose(fid);
    case '.mat'
      ll = load(gtdatafile);
      ispos = ismember(ll.divisionAnnotations(:,2),[1,2,3,4,5,103]);
      labellocs(end+1:end+nnz(ispos),:) = ll.divisionAnnotations(ispos,[3,4,5,1]).*gtratios(i,:);
      datafilei(end+1:end+nnz(ispos)) = i;

  end
end

%% load in all timepoints with annotations

timepoints_annotated = unique(labellocs(:,4));
timepoints_pa = intersect(timepoints_annotated,predtimestamps);

i = 1;
[rawreadframe,rawnframes,rawfid,rawheaderinfo] = get_readframe_fcn(inputdatafiles{i});
rawsz = [rawheaderinfo.xyzct(1:3),numel(timepoints_pa)];
%allpreddata = sparse(prod(rawsz),1);
allpredidx = zeros(0,5);

for i = 1:numel(timepoints_pa),

  ipred = find(predtimestamps==timepoints_pa(i));
  
  preddata = h5read(preddatafiles{ipred},preddatasetname);
  if isempty(preddata),
    continue;
  end
  preddata = preddata([3,2,1,4],:)';
  idxthresh = preddata(:,4) >= predthresh;
  fprintf('i = %d, nthresh = %d / %d\n',i,nnz(idxthresh),size(preddata,1));
  allpredidx(end+1:end+nnz(idxthresh),:) = [preddata(idxthresh,1:3),zeros(nnz(idxthresh),1)+timepoints_pa(i),preddata(idxthresh,4)];
%   idx = sub2indv(rawsz,[preddata(idxthresh,1:3),zeros(nnz(idxthresh),1)+i]);
%   allpreddata(idx) = preddata(idxthresh,4);
  
end

%% compute max projections

for i = 1:numel(timepoints_pa),
  
  t = timepoints_pa(i);
  %idx = allpredidx(:,4)==timepoints_pa(i);
  predmaxv = ConstructPredMaxProj2(allpredidx,rawsz(1:3),timepoints_pa(i));
  figure(i);
  imagesc(predmaxv); axis image;
  truesize;
  hold on;
  
  labelis = find(labellocs(:,4) == t);
  for ii = 1:numel(labelis),
    labeli = labelis(ii);
    plot(labellocs(labeli,2),labellocs(labeli,1),'o','Color',[1,.6,0]);
  end

end

%% filtering & non-maximal suppression

% nonmaximal supporession
sz = [rawsz(1:3),max(allpredidx(:,4))];
filrad = [5,5,5,0]; % probably make t non-zero when all timepoints have results

nonmaxpredidx0_pert = {};
nonmaxpredidx_pert = {};
scores_pert = {};
peaks_pert = {};

for ti = 1:numel(timepoints_pa),
  t = timepoints_pa(ti);

  chunksize = [100,100,50,2];
  chunkstep = chunksize/2;
  
  idxcurr = allpredidx(:,4)==t;
  newallpredidxcurr = imregionalmax_list(allpredidx(idxcurr,:),[sz(1:3),t],'filrad',filrad,'chunksize',chunksize,'chunkstep',chunkstep);
  nonmaxpredidx0_pert{ti} = newallpredidxcurr;

  chunksize = [100,100,50,4];
  chunkstep = chunksize/2;

  newallpredidx2curr = bwconncomp_list(newallpredidxcurr,[sz(1:3),t],'chunksize',chunksize,'chunkstep',chunkstep);
  nonmaxpredidx_pert{ti} = newallpredidx2curr;
    
  scores_pert{ti} = newallpredidx2curr(:,end);
  peaks_pert{ti} = newallpredidx2curr(:,1:4);
  
  save detections3.mat nonmaxpredidx0_pert nonmaxpredidx_pert scores_pert peaks_pert ti
end
  
%% plot detections

ti = 2;
thresh = .001;
boxrad = [100,100,5,3];
figpos = [10,10,2400,760];
nplot = 5;

for ti = 1:numel(timepoints_pa),

  [sortedscores,order] = sort(scores_pert{ti},1,'descend');

  for detorderi = 1:min(5,numel(sortedscores)),
    
    deti = order(detorderi);
    loc = peaks_pert{ti}(deti,:);
    
    hfig = 10;%+detorderi;
    
    textstr0 = sprintf(' Det rank %d, score = %.1f',detorderi,scores(deti));
    
    filename = fullfile(savedir,sprintf('Rank%ddet%d_score%.1f_t%dx%dy%dz%d.png',...
      detorderi,deti,scores(deti),round(loc(3)),round(loc(1)),round(loc(2)),round(loc(3))));
    
    PlotDivision(allinputdatafiles,alltimestamp,allpredidx,loc,rawsz,...
      'scores_pert',scores_pert,'peaks_pert',peaks_pert,'timepoints_pa',timepoints_pa,...
      'labellocs',labellocs,...
      'thresh',thresh,...
      'boxrad',boxrad,...
      'hfig',hfig,...
      'textstr0',textstr0,...
      'figpos',figpos,'filename',filename);
    
    drawnow;
    
  end
  
end

%% plot labels

thresh = .001;
boxrad = [100,100,5,3];
figpos = [10,10,2400,760];

labelidx = find(ismember(labellocs(:,4),timepoints_pa));

for labelii = 1:numel(labelidx),
  
  labeli = labelidx(labelii);
  loc = labellocs(labeli,:);
  t = loc(4);
  
  hfig = 20;
    
  filename = fullfile(savedir,sprintf('Label%d_t%dx%dy%dz%d.png',...
    labeli,round(loc(3)),round(loc(1)),round(loc(2)),round(loc(3))));
  
  PlotDivision(allinputdatafiles,alltimestamp,allpredidx,loc,rawsz,...
    'scores_pert',scores_pert,'peaks_pert',peaks_pert,'timepoints_pa',timepoints_pa,...
    'labellocs',labellocs,...
    'thresh',thresh,...
    'boxrad',boxrad,...
    'hfig',hfig,...
    'figpos',figpos,'filename',filename);
  
  drawnow;

end

%% compute accuracy at different thresholds

test_timepoints = [120,240,360];
labeldatafiles = [2,3];

istestlabel = ismember(labellocs(:,4),test_timepoints) & ismember(datafilei,labeldatafiles);
testpredidx = find(ismember(timepoints_pa,test_timepoints));

thresholds_try = logspace(-5,log10(.5),100);
match_rad = [20,20,10,2];

scores_cat = cat(1,scores_pert{testpredidx});
peaks_cat = cat(1,peaks_pert{testpredidx});
labellocs_cat = labellocs(istestlabel,:);

hfig = 122;

[precision,recall,ntruepos,nfalsepos,nfalseneg,hfig,hax] = ...
  MatchLabelsAndPredictions(labellocs_cat,peaks_cat,scores_cat,thresholds_try,...
  'match_rad',match_rad,'hfig',hfig,'doplot',true,'plotcolor','k');

%% train timepoints

timepoints_pred = timepoints_pa(~cellfun(@isempty,scores_pert));
train_timepoints = setdiff(timepoints_pred,test_timepoints);

istrainlabel = ismember(labellocs(:,4),train_timepoints) & ismember(datafilei,labeldatafiles);
trainpredidx = find(ismember(timepoints_pa,train_timepoints));

scores_cat = cat(1,scores_pert{trainpredidx});
peaks_cat = cat(1,peaks_pert{trainpredidx});
labellocs_cat = labellocs(istrainlabel,:);

[precision_train,recall_train,ntruepos_train,nfalsepos_train,nfalseneg_train,hfig,hax] = ...
  MatchLabelsAndPredictions(labellocs_cat,peaks_cat,scores_cat,thresholds_try,...
  'match_rad',match_rad,'hfig',hfig,'hax',hax,'doplot',true,'plotcolor','r');

%% precision-recall for tgmm

load(tgmmfile);

idxtest = ismember(tgmmDivisions(:,4),test_timepoints);
npredtgmm = nnz(idxtest);
istestlabel = ismember(labellocs(:,4),test_timepoints) & ismember(datafilei,labeldatafiles);
labellocs_cat = labellocs(istestlabel,:);

[precision_tgmm,recall_tgmm,ntruepos_tgmm,nfalsepos_tgmm,nfalseneg_tgmm] = ...
  MatchLabelsAndPredictions(labellocs_cat,tgmmDivisions(idxtest,:),ones(nnz(idxtest),1),.5,...
  'match_rad',match_rad,'doplot',false);

%% results per timepoint

labeldatafiles = [2,3];
thresholds_try = logspace(-5,log10(.5),100);
match_rad = [20,20,10,2];
colors = jet(numel(timepoints_pa))*.7;

precision_pert = nan(numel(thresholds_try),numel(timepoints_pa));
recall_pert = nan(numel(thresholds_try),numel(timepoints_pa));
ntruepos_pert = nan(numel(thresholds_try),numel(timepoints_pa));
nfalsepos_pert = nan(numel(thresholds_try),numel(timepoints_pa));
nfalseneg_pert = nan(numel(thresholds_try),numel(timepoints_pa));

for i = 1:numel(timepoints_pa),

  t = timepoints_pa(i);

  islabel = (labellocs(:,4)==t) & ismember(datafilei,labeldatafiles);
  
  [precision_pert(:,i),recall_pert(:,i),ntruepos_pert(:,i),nfalsepos_pert(:,i),nfalseneg_pert(:,i),hfig,hax] = ...
    MatchLabelsAndPredictions(labellocs(islabel,:),peaks_pert{i},scores_pert{i},thresholds_try,...
    'match_rad',match_rad,'hfig',hfig,'hax',hax,'doplot',true,'plotcolor',colors(i,:));
end


%% plot results

hfig = 123;
figure(hfig);
clf;
h = [];
legs = {};
legs{end+1} = sprintf('CNN, test timepoints %s',mat2str(test_timepoints));
h(end+1) = plot(recall,precision,'ko-','MarkerFaceColor','k','LineWidth',2);
hold on;
h(end+1) = plot(recall_train,precision_train,'ro-','MarkerFaceColor','r','LineWidth',2);
legs{end+1} = sprintf('CNN, train timepoints %s',mat2str(train_timepoints));
h(end+1) = plot(recall_tgmm,precision_tgmm,'co','MarkerFaceColor','c','LineWidth',2);
legs{end+1} = sprintf('TGMM, test timepoints %s',mat2str(test_timepoints));
trainmarker = '.';
trainlinestyles = {':','--'};
testmarker = '.';
testlinestyle = '-';

for i = 1:numel(timepoints_pa),
  legs{end+1} = sprintf('CNN, timepoint %d',timepoints_pa(i));
  if ismember(timepoints_pa(i),test_timepoints),
    marker = [testmarker,testlinestyle];
  else
    j = find(timepoints_pa(i) == train_timepoints);
    marker = [trainmarker,trainlinestyles{mod(j-1,numel(trainlinestyles))+1}];
  end
  h(end+1) = plot(recall_pert(:,i),precision_pert(:,i),marker,'Color',colors(i,:));
end

xlabel('Recall');
ylabel('Precision');
legend(h,legs);
box off;
set(gca,'XLim',[-.01,1.01],'YLim',[-.01,1.01]);

%% output table of results

fid = fopen('Detections20170617/Threshold2PrecisionRecall.csv','w');
fprintf(fid,'Threshold,TestPrecision,TestRecall,TrainPrecision,TrainRecall\n');
for i = 1:numel(thresholds_try),
  fprintf(fid,'%e,%e,%e,%e,%e\n',thresholds_try(i),precision(i),recall(i),...
    precision_train(i),recall_train(i));
end
fclose(fid);

fid = fopen('Detections20170617/Predictions.csv','w');
fprintf(fid,'X,Y,Z,T,Score\n');
for ti = 1:numel(peaks_pert),
  t = timepoints_pa(ti);
  for i = 1:size(peaks_pert{ti},1),
    fprintf(fid,'%f,%f,%f,%f,%e\n',peaks_pert{ti}(i,1),...
      peaks_pert{ti}(i,2),peaks_pert{ti}(i,3),t,scores_pert{ti}(i));
  end
end
fclose(fid);



%%

save mk5_large_small_batch_balanced_good_reweight.mat 