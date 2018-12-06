clear all
close all
clc

load('compEx1data.mat')
img_1 = imread('kronan1.JPG');

normalization = true;

x1 = cell2mat(x(1));
x2 = cell2mat(x(2));

if normalization
    for i=1:2
    %Get the mean, std and compute s:
    m_ = mean(x {i }(1:2 ,:) ,2); %Computes the mean value of the 1 st and 2 nd rows ow x{ i}
    std_ = std(x {i }(1:2 ,:) ,0 ,2);% Standard deviation of the 1 st and 2 nd rows ow x{ i}
    s = 1./std_;
    N = [s(1),0,-s(1)*m_(1);0,s(2),-s(2)*m_(2);0,0,1];
        if i==1
            N1 = N;
            x1_norm = N1*x1;
        else
            N2 = N;
            x2_norm = N2*x2;
        end
    end
else
        N = eye(3);
        N1= N; N2 = N;
        x1_norm = N*x1;
        x2_norm = N*x2;
end

M = zeros(size(x1,2),9);
for i =1:size(x1,2)
    xc1 = x1_norm(1:3,i);
    xc2 = x2_norm(1:3,i);
    M(i,1:end) = [xc2(1)*xc1(1),xc2(1)*xc1(2),xc2(1)*xc1(3),xc2(2)*xc1(1),xc2(2)*xc1(2),xc2(2)*xc1(3),xc2(3)*xc1(1),xc2(3)*xc1(2),xc2(3)*xc1(3)];
end


[U ,S ,V] = svd (M); % Computes the singular value decomposition of M
sol = V(1:end,end);
fprintf("\nthe ||v|| is :%f",norm(sol));
fprintf("\nthe eigenvalue is :%f",S(end,end));
fprintf("\nthe ||Mv|| is :%f",norm(M*sol));

F_old = reshape( sol ,[3 3])'
[U ,S ,V] = svd (F_old);
S(3,3) = 0;
F_old = U*S*V' ;
fprintf("determinant of F Normalized is %d\n",det(F_old)  )
indice = 2;
fprintf("Checking for the %d point the epipolar constrain we get: %d\n",...
        indice,x2_norm(1:3,indice)'*F_old*x1_norm(1:3,indice) )

%Getting the original Matrix:    
F = N2'*F_old*N1 ;
F = F./F(3,3);

l = F*x {1}; % Computes the epipolar lines
l = l ./ sqrt ( repmat ( l (1 ,:).^2 + l (2 ,:).^2 ,[3 1]));

hist ( abs ( sum ( l .* x {2})) ,100);

