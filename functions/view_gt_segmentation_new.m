function [gt_imgs gt_cnt] = view_gt_segmentation_new(bsdsRoot,img,present,out_path,only_name,para,save)

if para.Nimgs==300
    gt_imgs = readSegs(bsdsRoot,'color',str2num(only_name));
elseif para.Nimgs==500
    gt_imgs = readSegss(bsdsRoot,present,only_name);
else
    gt_imgs = readSegsss(bsdsRoot,present,only_name);
end
    
gt_path = fullfile(out_path, 'gt'); 
if  ~exist(gt_path), mkdir(gt_path);end

gt_cnt = [];
for i=1:size(gt_imgs,2)
    if save == 1
        [imgMasks,segOutline,imgMarkup]=segoutput(img,gt_imgs{i});
        imwrite(imgMarkup,fullfile(gt_path, [only_name, '_', int2str(i), '.bmp'])); 
    end
    if iscell(gt_imgs)
        gt_cnt(i) = max(gt_imgs{i}(:));
    else
        gt_cnt(i) = max(gt_imgs(:));
    end
end
