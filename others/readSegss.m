function [segs] = readSegss(bsdsRoot,present,iid)

files = fullfile(bsdsRoot,'groundTruth',present,[num2str(iid),'.mat']);
segs_gt = load(files);
segs = cell(1,length(segs_gt.groundTruth));
for i = 1:length(segs_gt.groundTruth)
    segs{1,i} = double(segs_gt.groundTruth{1,i}.Segmentation);
end