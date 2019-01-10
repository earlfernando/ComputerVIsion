clear all 
close all

%Run previous scripts:
run('Computer_Exercise_3.m')
run('Computer_Exercise_4.m')
close all
one_or_two = 1;
not_good = false;
actual_1 =P_1;
actual_2 =P_2;

%Compute X:
for i =1:size(x1,2)

    point_1 = [x1(1:2,i);1];
    point_2 = [x2(1:2,i);1];

    M = [P_1;P_2];
    M(1:3,5)= -[point_1];
    M(4:end,6)= -[point_2];

    [U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M
    sol = V(1:end,end);

    point_3D = sol(1:4);
    X(1:4, i) = pflat(point_3D); 
    
end

xproj1 = pflat( P_1 * X );
xproj2 = pflat( P_2 * X );
if not_good
    if one_or_two==1
        figure(1);
        imshow(img_1);
        hold on;
        plot(xproj1(1,:),xproj1(2,:),'+g','Markersize',3);
        plot(x1(1,:),x1(2,:),'ro')
        hold off;
    else


        figure(2)
        imshow(img_2);
        hold on;
        plot(xproj2(1,:),xproj2(2,:),'+g','Markersize',3);
        plot(x2(1,:),x2(2,:),'ro')

        hold off;
    end
        % Computes the projections
    good_points = ( sqrt( sum(( x1 - xproj1(1:2 ,:)).^2)) < 3 & ...
    sqrt( sum(( x2 - xproj2(1:2 ,:)).^2)) < 3);
    % Finds the points with reprojection error less than 3 pixels in both images
    X = X(: , good_points );
%     x1 = x1(: , good_points );
%     x2 = x2(: , good_points );

    xproj1 = pflat( P_1 * X );
    xproj2 = pflat( P_2 * X );
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Remove this to all matcjes

    % Computes the projections
    good_points = ( sqrt( sum(( x1 - xproj1(1:2 ,:)).^2)) < 3 & ...
    sqrt( sum(( x2 - xproj2(1:2 ,:)).^2)) < 3);
    % Finds the points with reprojection error less than 3 pixels in both images
    X = X(: , good_points );

%     x1 = x1(: , good_points );
%     x2 = x2(: , good_points );

    xproj1 =pflat( P_1 * X );
    xproj2 =pflat( P_2 * X );
    
    X=X(1:3,:);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot the two images
    if one_or_two==1
        figure(1);
        imshow(img_1);
        hold on;
        plot(xproj1(1,:),xproj1(2,:),'+g','Markersize',3);
        plot(x1(1,:),x1(2,:),'ro')
        hold off;
    else


        figure(2)
        imshow(img_2);
        hold on;
        plot(xproj2(1,:),xproj2(2,:),'+g','Markersize',3);
        plot(x2(1,:),x2(2,:),'ro')

        hold off;
    end
end



figure(3);
plot3(X(1,:),X(2,:),X(3,:),'.r','Markersize',2);
hold on;

plot3([Xmodel(1,startind);Xmodel(1,endind)],[Xmodel(2,startind);Xmodel(2,endind)],[Xmodel(3,startind);Xmodel(3,endind)],'b-');
grid on;
plotcams({actual_1,actual_2})
axis equal

hold off;



function e = compute_error(xm, xp)
 
n = size(xm,2);
e = norm(xm-xp)^2;
e = sqrt(e/n);

end