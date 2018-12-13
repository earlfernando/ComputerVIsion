clear all
close all
A=imread('a.jpg');
B=imread('b.jpg');
[ fA dA ] = vl_sift ( single ( rgb2gray ( A )) );
[ fB dB ] = vl_sift ( single ( rgb2gray ( B )) );
matches = vl_ubcmatch ( dA , dB );
xA = fA (1:2 , matches (1 ,:));
xB = fB (1:2 , matches (2 ,:));
total_points = size(matches,2);
points_subset =4;

best_inlier_indicies = [];
best_sum = 0;
best_homography=[];
dimension_of_points = 3;
d=0;

while(best_sum)<=200
  randind = randi(total_points,[points_subset,1]);  
  ground = [xA(:,randind);ones(1,points_subset)];
  truth = [xB(:,randind);ones(1,points_subset)];
  M = zeros(3*points_subset,3*dimension_of_points+points_subset);
  for i= 1: points_subset
    M(3*(i-1)+1,1:dimension_of_points)=ground(:,i)';
    M(3*(i-1)+2,dimension_of_points+1:2*dimension_of_points)=ground(:,i)';
    M(3*i,dimension_of_points*2+1:dimension_of_points*3)=ground(:,i)';
    M(3*(i-1)+1:3*i,3*dimension_of_points+i)=-truth(:,i);
  end
  [U S V] = svd(M);
  sol = V(:,end);
  H = reshape(sol(1:9),[3 3])';
  inliers =0;
  inlier_indicies=[];
  for v =1:total_points
      if norm(pflat(H*[xA(:,v);1])-[xB(:,v);1])<=5
          inliers = inliers+1;
          inlier_indicies = [inlier_indicies v];
      end
  end
  if inliers > best_sum
      best_sum = inliers;
      best_homography = H;
      best_inlier_indicies = inlier_indicies;
      
  end
  d=d+1;
end

tform = maketform('projective', best_homography');
transfbounds = findbounds(tform ,[1 1; size(A ,2) size(A ,1)]);
xdata =[ min([transfbounds(: ,1); 1]) max([ transfbounds(: ,1); size(B ,2)])];
ydata =[ min([transfbounds(: ,2); 1]) max([ transfbounds(: ,2); size(B ,1)])];

[newA] = imtransform(A ,tform , 'xdata', xdata , 'ydata ', ydata );
tform2 = maketform('projective', eye(3));
[ newB ] = imtransform(B , tform2 , 'xdata', xdata ,'ydata', ydata ,'size', size( newA ));

newAB = newB ;
newAB ( newB < newA ) = newA ( newB < newA );
%display stiching
figure;
imshow(newAB)
