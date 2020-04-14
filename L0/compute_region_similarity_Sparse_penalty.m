function W = compute_region_similarity_Sparse_penalty(histogram,L,Centroid,Area)
   % smooth features by the KDE method
    histogram1 = [];
    for j=1:size(histogram,2)
        histogram2 = histogram(:,j); 
        histogram2 = histogram_smoothing(histogram2, 30, 20);
        histogram1 = [histogram1,histogram2(:,2)];
    end
    histogram = histogram1';
    histogram(:,1) = mean(histogram,2);
    [m,n] = size(histogram);
    d = sqrt( sum(histogram.^2,1) );
    d(d<eps)=1;
    histogram = histogram ./ repmat( d, m,1 );
    W_Y = [];
    Dictionary = [];
    % and uses all the cores of the machine
    for i = 1:n
        y = histogram(:,i);
        Dictionary = histogram;
        Dictionary(:,i) = zeros(m,1);
%         x = OMP(Dictionary,y,L);
        x = CS_SP(Dictionary,y,L);
        ind_x=find(x~=0);
        C_neg=Centroid(ind_x,:);
        C_cur=Centroid(i,:);
        D_C=sqrt(sum((repmat(C_cur,size(C_neg,1),1)-C_neg).*(repmat(C_cur,size(C_neg,1),1)-C_neg),2));

        ind_d=find(D_C);
        if numel(ind_d)>1
            idxD_C = find(D_C==max(D_C));
            for t=1:length(idxD_C)
                bb= ind_d == idxD_C(t);
                ind_d(bb)=[];
            end
            %ind_x(ind_x(find(D_C==max(D_C))))=[];
            %     [~,idx]=sort(D_C,'descend');
            ind_new=ind_x(ind_d);
            Dictionary_new=histogram(:,ind_new);
    %         Dictionary_new=Dictionary;
    %         for kk=1:numel(ind_new)
    %             Dictionary_new(:,ind_new(kk))=zeros(m,1);
    %         end
            %  Dictionary(:,i) = zeros(m,1);
             AreaWeights=ones(numel(ind_new),1);
             ind_area=find(Area(ind_new)<200);
             if numel(ind_area)>1
                %     ind_area=find(Area(ind_new));
                   AreaWeights(ind_area)=2;
             end    
    %         AreaWeights=1-AreaWeights./sum(AreaWeights);
             x1 = pinv(Dictionary_new)*y.*AreaWeights;
    %         x1 = OMP(Dictionary_new,y,2);
    %         x0=x1;
            %     x=sparse(x);
            x0=sparse(n,1);
            if x1~=0, x0(ind_new)=x1(x1~=0);end
        else
            x0=x;
        end
        R_Lable = [1:n];
        [Region_Coeffi] = ProjectCoefficient(x0,R_Lable);
        Region = Dictionary*Region_Coeffi;
        Error_Region = repmat(y,1,n) - Region;
        similarity = [];
        similarity = sqrt( sum(Error_Region.^2,1) );
        similarity = (1 - similarity);
        similarity(i) = 0;
        index = similarity<0.00001;
        similarity(index) = 0;
        W_Y =[W_Y; similarity];
    end
    W = W_Y + W_Y';
    index = find(W>1);
    W(index) =  W(index)/2;


    % WWY = symetricy_WY(WY);
    % W_Y = (W_Y + W_Y')/2;

    end