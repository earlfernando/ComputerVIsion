clear all 
close all

load('compEx3data.mat')
img_1 = imread('cube1.JPG');
img_2 = imread('cube2.JPG');

points_1 = cell2mat(x(1));
points_2 = cell2mat(x(2));

i = 1;
%Get the mean, std and compute s:
m_ = mean(x {i }(1:2 ,:) ,2); % Computes the mean value of the 1 st and 2 nd rows ow x{ i}
std_ = std(x {i }(1:2 ,:) ,0 ,2); % Standard deviation of the 1 st and 2 nd rows ow x{ i}
s = 1./std_;

%Construct transformation matrix N
N = [s(1),0,-s(1)*m_(1);0,s(2),-s(2)*m_(2);0,0,1];
%Get the normalized points
points_1_normalized = N*points_1;
%Plot those points
plot(points_1_normalized(1,:),points_1_normalized(2,:),'*')
%Check the mean distance from the center:
mean(sqrt(points_1_normalized(1,:).^2 + points_1_normalized(2,:).^2) )

%Construct M:
M = zeros()
for i = 1:37
    
    
    
end


% [U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M