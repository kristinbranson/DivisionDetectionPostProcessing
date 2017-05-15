addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;
addpath /groups/branson/home/bransonk/codepacks/keller-lab-block-filetype/matlabWrapper/
datadir = '/groups/branson/home/bransonk/tracking/code/Lineaging/MouseEmbryo/DivisionDetectionData_14-05-21/Predictions_0427';

%% find data files
inputdatafiles = mydir(fullfile(datadir,'*00.klb'));
preddatafiles = cell(size(inputdatafiles));
ismissing = false(1,numel(inputdatafiles));
for i = 1:numel(inputdatafiles),
  preddatafile = [inputdatafiles{i}(1:end-5),'1.klb'];
  if ~exist(preddatafile,'file'),
    ismissing(i) = true;
    [~,n] = myfileparts(inputdatafiles{i});
    fprintf('Could not find predictions for %s\n',n);
  else
    preddatafiles{i} = preddatafile;
  end
end

inputdatafiles(ismissing) = [];
preddatafiles(ismissing) = [];

%% load one

filenum = 1;
[rawreadframe,rawnframes,rawfid,rawheaderinfo] = get_readframe_fcn(inputdatafiles{filenum});
[predreadframe,prednframes,predfid,predheaderinfo] = get_readframe_fcn(preddatafiles{filenum});

