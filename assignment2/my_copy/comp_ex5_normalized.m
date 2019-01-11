clear all 
close all

%Run previous scripts:
run('Computer_Exercise_3.m')
run('Computer_Exercise_4.m')
close all
one_or_two = 2;
not_good = true;
actual_1 =P_1;
actual_2 =P_2;
P_1 = inv(K_1)*P_1;
P_2= inv(K_2)*P_2;
%Compute X:
for i =1:size(x1,2)

    point_1 = [x1(1:2,i);1];
    point_2 = [x2(1:2,i);1];
    point_1=pflat( inv(K_1)*point_1);
    point_2 = pflat(inv(K_2)*point_2);
    x1(1:2,i)=  point_1(1:2,:);
    point_1 = point_1(1:2,:);
    x2(1:2,i) = point_2(1:2,:);
    point_2 = point_2(1:2,:);
    M = [P_1;P_2];
    M(1:3,5)= -[point_1;1];
    M(4:end,6)= -[point_2;1];

    [U ,S ,V] = svd ( M ); % Computes the singular value decomposition of M
    sol = V(1:end,end);

    point_3D = sol(1:4);
    X(1:4, i) = pflat(point_3D); 
    
end

xproj1 = pflat( P_1 * X );
xproj2 = pflat( P_2 * X );
    
if not_good
    good_points = ( sqrt( sum(( x1 - xproj1(1:2 ,:)).^2)) < 3/3000 & ...
    sqrt( sum(( x2 - xproj2(1:2 ,:)).^2)) < 3/3000);
    % Finds the points with reprojection error less than 3 pixels in both images
    
    xproj1 = K_1*xproj1;
    xproj2 =K_2*xproj2;
    x1 = [x1;ones(1,size(x1,2))];
    x2 = [x2;ones(1,size(x2,2))];
    x1= K_1*x1(1:3,:);
    x2 = K_2*x2(1:3,:);
    % error calculation 
x1_error = x1-xproj1;
x2_error = x2-xproj2;
x1_true = ~isnan(x1_error);
x2_true = ~isnan(x2_error);
n_1 =size(x1_true,2);
n_2= size(x2_true,2);
e1_n = sqrt(norm(x1_error(x1_true),2));
e2_n = sqrt(norm(x2_error(x2_true),2));
error_1 = compute_error(x1,xproj1);
error_2 = compute_error(x2,xproj2);
    
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

%     x1 = x1(: , good_points );
%     x2 = x2(: , good_points );
X = X(: , good_points );
    xproj1 =K_1* pflat( P_1 * X );
    xproj2 =K_2* pflat( P_2 * X );
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Remove this to all matcjes

    % Computes the projections
    good_points = ( sqrt( sum(( x1 - xproj1(1:2 ,:)).^2)) < 3/3000 & ...
    sqrt( sum(( x2 - xproj2(1:2 ,:)).^2)) < 3/3000);
    % Finds the points with reprojection error less than 3 pixels in both images
    X = X(: , good_points );

%     x1 = x1(: , good_points );
%     x2 = x2(: , good_points );

    xproj1 =K_1* pflat( P_1 * X );
    xproj2 =K_2* pflat( P_2 * X );
    x1 = [x1;ones(1,size(x1,2))];
    x2 = [x2;ones(1,size(x2,2))];
    x1= K_1*x1(1:3,:);
    x2 = K_2*x2(1:3,:);
    
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
limit = 30;
xlim([-limit limit])
ylim([-limit limit])
zlim([-limit limit])

hold off;



function e = compute_error(xm, xp)
truth =xm-xp;
ind=~isnan(truth);
diff = truth(ind);
n = size(diff,2);
e = norm(diff)^2;
e = sqrt(e/n);

end