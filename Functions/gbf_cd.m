function change_map = gbf_cd(Xl, locations)

modalities = 2;

[rows,cols,~] = size(Xl{1});

Xl_AA{1} = Xl{1}(locations)'; %Ir
Xl_AA{2} = Xl{2}(locations)'; %Ig

n = length(locations);

% local graph
wl = cell(modalities,1);

complement = setdiff(1:(rows*cols), locations);

metric = @(x,kernelstda) exp(-(x.^2) ./ ((kernelstda ^ 2)));


for i = 1 : modalities
    
    %% L2 norm
    
    distlAA = pdist2(Xl_AA{i}',Xl_AA{i}','euclidean');
    distlAB = pdist2(Xl{i}(complement)',Xl_AA{i}','euclidean').^3;
    
    
    %% L1 norm
    %                 x=repmat(Xl_AA{i}',[size(Xl_AA{i}',1),1]);
    %                 y=repmat(Xl_AA{i},[size(Xl_AA{i},2),1]);
    %                 y = y(:);
    %                 distlAA = reshape(sum(abs(x-y),2).^2,n,n);
    %                 x=repmat(Xl_AA{i},[size(Xl{i}(complement),2),1]);
    %                 x = x(:);
    %                 y=repmat(Xl{i}(complement)',[size(Xl_AA{i}',1),1]);
    %                 distlAB = reshape(sum(abs(x-y),2).^2,length(complement),n);
    %
    %% Normalization and Degree
    
    D1 = distlAA*ones(n,1) + distlAB'*ones(length(complement),1);%diag(sum(distlAA, 2));
    distlAA = distlAA./repmat((D1),1,n);
    D2 = distlAB*ones(n,1) + (distlAB*pinv(distlAA))*(distlAB'*ones(length(complement),1));%sum(distlAB, 2);
    distlAB = distlAB./repmat((D2),1,n);
    
    
    %clear x y D1 D2
    
    %% Metric computation
    
    %auxstd = [distlAA ; distlAB];
    %mean data
    kernelstd = mean(distlAB(:));
    %silverman
    %kernelstd = 2.07*std(distlAB(:))/(n^(1/5));
    %clear auxstd
    
    distlAA = metric(distlAA,kernelstd);
    distlAB = metric(distlAB,kernelstd);
    
    
    wl{i} = [distlAA ;distlAB];
    
    %clear distlAA distlAB distlAAsmooth
    
    
end
clear Xl_AA x y D1 D2 distlAA distlAB distlAAsmooth
%% Multimodal Weights

n = length(locations);
WL = min(cat(3,wl{1} , wl{2}),[],3);


% One shot method of Nystrom
W_AA = WL(1:n,1:end);
W_BA = WL(n+1:end,1:end);


%% Nystrom aproximation

W_AA_sqrtinv = pinv(sqrtm(W_AA));%W_AA^-0.5;
%
S = W_AA + (W_AA_sqrtinv*(W_BA'*W_BA)*W_AA_sqrtinv);

[U_s,D_s] = eig(S);

Uhat_W = [ W_AA ; (W_AA_sqrtinv*W_BA')']*(U_s*pinv(sqrtm(D_s)));%%%[U_s; W_BA*U_s*(pinv(D_s))];%[W_AA;W_BA*W_AA_sqrtinv]*(U_s*pinv(sqrtm(D_s)));%%%

clear S W_AA_sqrtinv W_AA W_BA WL

%E = zeros(1,n);
MI = zeros(1,n);
prior1  = (Xl{1} - Xl{2})./((Xl{1} + Xl{2}));
prior2  = (Xl{2} - Xl{1})./((Xl{1} + Xl{2}));
prior = imbinarize(prior1) + imbinarize(prior2);


for j = 1 : n
    Iaux = Uhat_W(:,j)*sqrt(D_s(j,j));
    
    A = Iaux(1:n);
    AB = Iaux(n+1:end);
    Iaux(locations) = A;
    Iaux(complement) = AB;
    
    Iaux = imbinarize(abs(Iaux));
    
    %     figure, imshow(Iaux,[]),colorbar%, colormap hot
    %     title(['Eigenvector sample ' , num2str(i) ])
    %     set(gca,'FontSize',12)
    
    MI(j) =  mutInfo(prior(:),Iaux);
    clear Iaux A AB;
end

%     figure, plot(MI,'linewidth',2)
%     title('MI of eigenvectors')
%     xlabel('$u_{i} \sqrt{d_{i}}$','FontSize',13,'interpreter','latex')
%     ylabel('$MI(I_{u_{i}},I_{Prior})$','FontSize',13,'interpreter','latex')
%     set(gca,'FontSize',12)
%
j = find(MI == max(MI),1,'first');
%
%     hold on
%
%     plot(i,MI(i),'ro','MarkerSize',12)

Iaux = Uhat_W(:,j)*sqrt(D_s(j,j));

A = Iaux(1:n);
AB = Iaux(n+1:end);
Iaux(locations) = A;
Iaux(complement) = AB;
Iaux = ((reshape(Iaux,rows,cols)));

change_map = imbinarize(abs(Iaux));
