%clear all 
%close all

load('compEx3data.mat')
img_1 = imread('cube1.JPG');
img_2 = imread('cube2.JPG');

[f1,d1]= vl_sift( single( rgb2gray( img_1 )) , 'PeakThresh', 1);
[f2,d2]= vl_sift( single( rgb2gray( img_2 )) , 'PeakThresh', 1);

imshow(img_1)
vl_plotframe( f1 );

[ matches , scores ] = vl_ubcmatch (d1 , d2 );


x1 = [ f1(1, matches(1 ,:)); f1(2, matches(1 ,:))];
x2 = [ f2(1, matches(2 ,:)); f2(2, matches(2 ,:))];


perm = randperm( size( matches ,2));
figure ;
imagesc([ img_1 img_2 ]);
hold on ;
plot([ x1(1 , perm (1:10)); x2(1 , perm (1:10))+ size( img_1 ,2)] , ...
[ x1(2 , perm (1:10)); x2(2 , perm (1:10))] , '-' );
hold off ;