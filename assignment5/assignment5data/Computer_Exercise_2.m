clear all
close all
run('Computer_Exercise_1.m');
im1=imread('kronan1.JPG');
im2=imread('kronan2.JPG');
load('compex2data.mat')
d = linspace (5 ,11 ,200);
sc = 0.25;
[ ncc , outside_image ] = compute_ncc (d , im2 , P {2} , im1 , segm_kronan1 , P {1} ,3 , sc );
[ maxval , maxpos ] = max ( ncc ,[] ,3);
disp_result ( im2 , P {2} , segm_kronan2 , d ( maxpos ) ,0.25 , sc );
