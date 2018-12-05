clear all
close all
clc

load('compEx1data.mat')
load('compEx3data.mat')

x1 = cell2mat(x(1));
x2 = cell2mat(x(2));

x1_n = inv(K)*x1;
x2_n = inv(K)*x2;


M = zeros(size(x1,2),9);
for i =1:size(x1,2)
    xc1 = x1_n(1:3,i);
    xc2 = x2_n(1:3,i);
    M(i,1:end) = [xc2(1)*xc1(1),xc2(1)*xc1(2),xc2(1)*xc1(3),xc2(2)*xc1(1),xc2(2)*xc1(2),xc2(2)*xc1(3),xc2(3)*xc1(1),xc2(3)*xc1(2),xc2(3)*xc1(3)];
end
[U ,S ,V] = svd (M); % Computes the singular value decomposition of M
sol = V(1:end,end);
fprintf("\nthe ||v|| is :%f",norm(sol));
fprintf("\nthe eigenvalue is :%f",S(end,end));
fprintf("\nthe ||Mv|| is :%f\n",norm(M*sol));

Eapprox= reshape( sol ,[3 3])' ;
[U ,S ,V] = svd ( Eapprox );
if det (U *V') >0
    E = U* diag ([1 1 0])* V';
else
    V = -V;
    E = U* diag ([1 1 0])* V';
end

fprintf("Epipolar constrains fullified? %d \n", (x2_n(1:3,9))'*E* x1_n(1:3,9) )

E = E./E(3,3);
F= inv(K)'*E*inv(K);


l = F*x1;


random_indexes = floor(rand(1,20)*2008) +1 ;

img = imread('kronan2.JPG');
imshow(img)
hold on
for i=random_indexes
    point = x2(1:3,i);
    l_point = l(1:3,i);
    plot(point(1),point(2),'xr','LineWidth',2.5)
    rital(l_point)
end

figure
hist ( abs ( sum ( l .* x {2})) ,100);


