close all;clear;
addpath 'others'
addpath 'evals'
addpath 'L0'
addpath 'functions'
addpath 'histsmooth'

%% set parameters for bipartite graph
para.alpha = 0.001; % affinity between pixels and superpixels
para.beta  =  20;   % scale factor in superpixel affinityI
para.nb = 1; % number of neighbors for superpixels
para.L = 3; % control the sparsity during solving the L0 problem 

%% read image
bsdsRoot='BSD';
load_file='bsd_300_feat';
outputpath = 'results';
fid = fopen('Nsegs.txt','r');
Nimgs = 300; % number of images in BSDS300

[BSDS_INFO] = fscanf(fid,'%d %d \n',[2,Nimgs]);
fclose(fid);
PRI_all = zeros(Nimgs,1);
VoI_all = zeros(Nimgs,1);
GCE_all = zeros(Nimgs,1);
BDE_all = zeros(Nimgs,1);

global_graph_mode={'L0'};
features_mode={'mlab'};
weight = 1;

for m=1:length(global_graph_mode)
    gmode = global_graph_mode{m};
    for idxI = 1:Nimgs
        tic; out_path= fullfile(outputpath,'BSDS300',gmode);
        if ~exist(out_path,'dir'), mkdir(out_path);  end
       %% locate image
        img_name = int2str(BSDS_INFO(1,idxI));
        img_loc = fullfile(bsdsRoot,'images','test',[img_name,'.jpg']);
        if ~exist(img_loc,'file')
            img_loc = fullfile(bsdsRoot,'images','train',[img_name,'.jpg']);
        end
        img = im2double(imread(img_loc)); [X,Y,~] = size(img);
        load_name=fullfile(load_file,[img_name '.mat']);
        load(load_name)
        
       %% construct graph
        Np = X*Y;   Nsp = 0;
        for k = 1:length(seg)
            Nsp = Nsp + size(seg{k},2);
        end      
        W_Y = sparse(Nsp,Nsp);
        edgesXY = [];
        j = 1;
        for k = 1:length(seg) 
            % for each over-segmentation
            temp = feat{k}.shape;           
            tmp1 = feat{k}.mlab;
            labels(idxI,k) = APclustering(feat{k});
            Center = clusteringcenter(tmp1',labels(idxI,k));
            index_tmp = litekmeans(tmp1,labels(idxI,k),'Start',Center');
            
            supixel_index = seg{k};
            centroid = temp(:,2:3);
            Area = temp(:,1);
            
            % superpixel division 
            local_nodes  = find(index_tmp == mode(index_tmp));
            global_nodes = find(index_tmp ~= mode(index_tmp));
            W_f = zeros(length(Area),length(Area));
            %you could change the feature descriptor here
            for n = 1:length(features_mode)
                fmode = features_mode{n};
                switch fmode
                    case 'mlab'
                        feature = feat{k}.mlab;
                    case 'ch'
                        feature = cat(2,feat{k}.chr,feat{k}.chg,feat{k}.chb);
                    case 'lbp'
                        feature = cat(2,feat{k}.lbpr(:,1:end-1),feat{k}.lbpg,feat{k}.lbpb);
                    case 'bow'
                        feature = feat{k}.siftBoW1;
                end
                feature(:,all(feature == 0, 1))=[];
                [fm,fn] = size(feature);
                feature=(feature-repmat(mean(feature),fm,1))./repmat(std(feature),fm,1); 

                % first we construct the adjacent graph over all nodes
                w = makeweights(seg_edges{k},feature,para.beta);
                W_local = adjacency(seg_edges{k},w);

                W = W_local;

                % randomly generate supperpxiels
                p = randperm(length(global_nodes));
                W_L0 = L0GRAPH(feature,para.L,centroid,Area);
                W = assignGraphValue(W,W_L0,p);
                W_f = W_f + W*weight(n);
            end
            W = sparse(W_f);
            Nk = size(seg{k},2); % number of superpixels in over-segmentation k
            W_Y(j:j+Nk-1,j:j+Nk-1) = prune_knn(W,para.nb);
            
            % affinities between pixels and superpixels
            for i = 1:Nk
                idxp = seg{k}{i}; % pixel indices in superpixel i
                Nki = length(idxp);
                idxsp = j + zeros(Nki,1);
                edgesXY = [edgesXY; [idxp, idxsp]];
                j = j + 1;
            end
        end
        W_XY = sparse(edgesXY(:,1),edgesXY(:,2),para.alpha,Np,Nsp);
        % affinity between a superpixel and itself is set to be the maximum 1.
        W_Y(1:Nsp+1:end) = 1;  B = [W_XY;W_Y];
        
       %% Transfer cut
        out_path_gt= fullfile(outputpath,'BSDS300', gmode, img_name);
        if ~exist(out_path_gt,'dir'), mkdir(out_path_gt); end
        [gt_imgs, gt_cnt] = view_gt_segmentation(bsdsRoot,img,BSDS_INFO(1,idxI),out_path_gt,img_name,0);
        E=[];   Nseg = BSDS_INFO(2,idxI);
        label_img = Tcut(B,Nseg,[X,Y]);  ti = toc; clear B;
        % display the result
        view_segmentation(img,label_img(:),out_path,img_name,0);
        
       %% Evaluation and save result
        out_vals = eval_segmentation(label_img,gt_imgs);  clear label_img gt_imgs;
        fprintf('%s %4d %6s: %2d %9.6f, %9.6f, %9.6f, %9.6f %.4fs\n',gmode,idxI,img_name,...
            Nseg, out_vals.PRI, out_vals.VoI, out_vals.GCE, out_vals.BDE, ti);
        PRI_all(idxI) = out_vals.PRI;
        VoI_all(idxI) = out_vals.VoI;
        GCE_all(idxI) = out_vals.GCE;
        BDE_all(idxI) = out_vals.BDE;
    end
    %%
    fprintf('%s Mean: %14.6f, %9.6f, %9.6f, %9.6f \n', gmode, mean(PRI_all), mean(VoI_all), mean(GCE_all), mean(BDE_all));
    fid_out = fopen(fullfile(outputpath,'BSDS300',gmode,'evaluation.txt'),'w');
    for idxI=1:Nimgs
        fprintf(fid_out,'%6d %9.6f, %9.6f, %9.6f, %9.6f \n', BSDS_INFO(1,idxI),...
            PRI_all(idxI), VoI_all(idxI), GCE_all(idxI), BDE_all(idxI));
    end
    fprintf(fid_out,'Mean: %10.6f, %9.6f, %9.6f, %9.6f \n', mean(PRI_all), mean(VoI_all), mean(GCE_all), mean(BDE_all));
    fclose(fid_out);
end
