function [segs] = readSegsss(bsdsRoot,present,iid)

files = fullfile(bsdsRoot,'clear_groundTruth',[iid,'_GT.mat']); %%%groundTruth
segs_gt = load(files);
segs = cell(1,length(segs_gt));
for i = 1:length(segs_gt)
    segs{1,i} = double(segs_gt.newseg);
end