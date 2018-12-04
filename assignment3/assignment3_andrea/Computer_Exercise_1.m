clear all
close all
clc

load('compEx1data.mat')
img_1 = imread('kronan1.JPG');

x1 = cell2mat(x(1));
x2 = cell2mat(x(2));

for i=1:2
    %Get the mean, std and compute s:
    m_ = mean(x {i }(1:2 ,:) ,2); %Computes the mean value of the 1 st and 2 nd rows ow x{ i}
    std_ = std(x {i }(1:2 ,:) ,0 ,2);% Standard deviation of the 1 st and 2 nd rows ow x{ i}
    s = 1./std_;
    N = [s(1),0,-s(1)*m_(1);0,s(2),-s(2)*m_(2);0,0,1];
    if i==1
        N1 = N;
        x1_norm = N*x1;
    else
        N2 = N;
        x2_norm = N*x2;
    end
end

M = zeros(size(x1,2),9);
for i =1:size(x1,2)
    
    xc1 = x1_norm(1:3,i);
    xc2 = x2_norm(1:3,i);
    M(i,1:end) = [xc2(1)*xc1(1),xc2(1)*xc1(2),xc2(1)*xc1(3),xc2(2)*xc1(1),xc2(2)*xc1(2),xc2(2)*xc1(3),xc2(3)*xc1(1),xc2(3)*xc1(2),xc2(3)*xc1(3)];
end


[U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M
sol = V(1:end,end);
fprintf("\nthe ||v|| is :%f",norm(sol));
fprintf("\nthe eigenvalue is :%f",S(end,end));
fprintf("\nthe ||Mv|| is :%f",norm(M*sol));

F = reshape( sol ,[3 3])';
det(F)                                                 %Maybe do it better
indice = 4;
x2_norm(1:3,indice)'*F*x1_norm(1:3,indice)

FF = N2'*F*N1 ; 
F = FF;
l = F*x {1}; % Computes the epipolar lines
l = l ./ sqrt ( repmat ( l (1 ,:).^2 + l (2 ,:).^2 ,[3 1]));

hist ( abs ( sum ( l .* x {2})) ,100);

