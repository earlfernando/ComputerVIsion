clear all;
loaded_image = imread("compEx2.JPG");
imagesc(loaded_image);
colormap(gray);
points = load('compEx2.mat');
point1 = points.p1;
point2 = points.p2;
point3 = points.p3;
points_arraynames= fieldnames(points);
all_points = [];
for i =1:numel(points_arraynames)
    local_points_pairs = points.(points_arraynames{i});
    local_size = size(local_points_pairs);
    local_size_columns = local_size(2);
    local_array = [];
    for j = 1:local_size_columns
        local_point = local_points_pairs(:,j);
%         local_point(end) = [];
        all_points = [all_points local_point];
        local_array = [local_array local_point];
        hold on ;
        plot(local_point(1),local_point(2),'r*');
    end
    coeff = polyfit(local_array(:,1)',local_array(:,2)',1);
    line =[ coeff(1) , 1 , coeff(2)]';
    hold on ;
    disp('rital');
    rital(line,'-o');
    
    
end
