clear all 
close all

%Run previous scripts:
run('Computer_Exercise_3.m')
run('Computer_Exercise_4.m')
close all

%Compute X:

for i =1:size(x1,2)

    point_1 = [x1(1:2,i);1];
    point_2 = [x2(1:2,i);1];
    M = [P_1;P_2];
    M(1:3,5)= -point_1;
    M(4:end,6)= -point_2;

    [U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M
    sol = V(1:end,end);

    point_3D = sol(1:4);
    X(1:4, i) = pflat(point_3D); 
    
end

xproj1 = pflat( P_1 * X );
xproj2 = pflat( P_2 * X );

figure
imshow(img_1)
hold on
for i=1:size(X,2)
    plot([xproj1(1,i);x1(1,i)],[xproj1(2,i);x1(2,i)],'o-')
end
    

if 1==0
    % Computes the projections
    good_points = ( sqrt( sum(( x1 - xproj1(1:2 ,:)).^2)) < 3 & ...
    sqrt( sum(( x2 - xproj2(1:2 ,:)).^2)) < 3);
    % Finds the points with reprojection error less than 3 pixels in both images
    X = X(: , good_points );
    % Removes points that are not good enough .

xproj1 = pflat( P_1 * X );
xproj2 = pflat( P_2 * X );

end