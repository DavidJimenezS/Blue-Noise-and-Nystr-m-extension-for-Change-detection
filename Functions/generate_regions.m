function [centroids, region1, region2] = generate_regions(before, after, stop, Xl)

[rows,cols,~] = size(before);

prior1  = (before(:,:,1) - after)./((before(:,:,1) + after));
prior2  = (after - before(:,:,1))./((before(:,:,1) + after));
prior = imbinarize(prior1) + imbinarize(prior2);

[~, B ,~, biaoji] = imgsegment(prior,100,stop);     %Segmentaci?n

figure
imshow(B,[]);

biaoji = reshape(biaoji,rows,cols);

a = min(min(biaoji));
b = max(max(biaoji));

if a == 0
    region1 = zeros(b+1,1);
    region2 = zeros(b+1,1);
    centroids = zeros(b+1,2);
    %kappas = zeros(b + 1,1);
else
    region1 = zeros(b,1);
    region2 = zeros(b,1);
    centroids = zeros(b,2);
    %kappas = zeros(1,b);
end


for j = a:b               %????????????????????????
    
    tmp = find(biaoji==j);
    [I , J] = ind2sub([rows , cols],tmp);
    
    if a == 0
        centroids(j + 1,1) = sum(I)/length(I);
        centroids(j + 1,2) = sum(J)/length(J);
    else
        centroids(j,1) = sum(I)/length(I);
        centroids(j,2) = sum(J)/length(J);
    end
    
    b1 = Xl{1}(tmp);
    b2 = Xl{2}(tmp);
    
    if a == 0
        region1(j + 1)= mean(b1(:));
        region2(j + 1)= mean(b2(:));
    else
        region1(j)= mean(b1(:));
        region2(j)= mean(b2(:));
    end
end


centroids = floor(centroids);