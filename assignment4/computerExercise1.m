
clear all
close all
load('compEx1data.mat')
%variables 
iterations = 10;
points_subset_planes = 3;


X= pflat(X);
total_points = length(X);


%Least squares calculations
plane = leastSquare(X,"overall");



%Ransac
best_inlier_indicies = [];
best_sum = 0;
best_plane=[];

for i =1 : iterations
    randind = randi(total_points,[points_subset_planes,1]);
    plane = null (X (:, randind )');
    plane = plane ./ norm ( plane (1:3));
    inliers = abs ( plane'* X ) <= 0.1;
    total_inlier = sum(inliers);
    inidicies = find(inliers==1);  
    if (best_sum<= total_inlier)
       best_sum = total_inlier;
       best_plane = plane;
       best_inlier_indicies = inliers;        
    end  
end
X_inliers_ransac= X(:,best_inlier_indicies);
ransac_rms = sqrt ( sum (( plane'* X_inliers_ransac).^2)/ size (X_inliers_ransac,2))
%plot
figure;
histogram(abs(plane'*X_inliers_ransac),100);


%Least squares on the ransac inliers
plane = leastSquare(X_inliers_ransac,"Ransac");
%plot
figure;
histogram(abs(plane'*X_inliers_ransac),100);

%inlier projection
P1 = cell2mat(P(1));
P2 = cell2mat(P(2));
x1 =pflat( P1* X_inliers_ransac);
x2 = pflat(P2 *X_inliers_ransac);
im1=imread('house1.jpg');
im2=imread('house2.jpg');
%plot
figure;
imshow(im1);
hold on
plot(x1(1,:),x1(2,:),'b.','markersize', 20)
hold off
figure(4);
imshow(im2);
hold on
plot(x2(1,:),x2(2,:),'b.','markersize', 20)

%Compute Homographies
norm_P1 = inv(K)*P1;
norm_P2 = inv(K)*P2;
R = norm_P2(1:3,1:3);
t = norm_P2(:,4);
pi = pflat(plane);
H = (R-t*pi(1:3)');
points1_2d = pflat(norm_P1*X);
points2_2d = pflat(norm_P2*X);
new_points2_2d = K*pflat(H*points1_2d);
points2_2d = K*points2_2d;
points1_2d = K*points1_2d;
%plot
figure;
imshow(im2);
hold on
plot(points2_2d(1,:),points2_2d(2,:),'bo');
plot(new_points2_2d(1,:),new_points2_2d(2,:),'g+');
figure;
hold on
imshow(im1);
hold on
plot(points1_2d(1,:),points1_2d(2,:),'bo');


%functions
function plane = leastSquare(X,term)
meanX = mean (X ,2);
Xtilde = (X-repmat(meanX,[1 size(X,2)]));
M = Xtilde (1:3 ,:)* Xtilde (1:3 ,:)';
[V , D ] = eig ( M );
coefficients = V(:,1); %smallest eigenvector of the M 
local_calculation =meanX(1:3);
d= -(coefficients'*local_calculation); %given in the optional exercise
plane = [coefficients; d];
plane = plane ./ norm ( plane (1:3));
disp(term);
RMS = sqrt ( sum (( plane'* X ).^2)/ size (X ,2));
disp(RMS);
end





