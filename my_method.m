clear all
close all
clc

data_dir = strcat(pwd,'/data/');
files = dir(strcat(data_dir,'data*.mat'));
image = dir(strcat(data_dir,'img*.jpg'));

for j =1 : numel(files)
    images{j}= imread(strcat(data_dir,image(j).name));
    data{j} = load(strcat(data_dir,files(j).name));
end
number_of_objects =numel( data{1}.U);
number_of_images = numel(files);
ld_iterations=30;
lambda = 0.0000000000000000001;
maxDistance = 0.01;
iterations= 10;
P_est ={};
poses ={};
bounding_boxes ={};
for i =1 : number_of_objects   
    for j = 1: number_of_images
    best_inlier = 0;
     % normalize the 3D points
    x = data{j}.U{i};
    meanU = mean(x,2);
    normalized3Dpoint=(x-repmat(meanU,[1 size(x,2)]));
    normalized2Dpoint =data{j}.u{i};
    ptCloudIn = pointCloud(normalized3Dpoint');
    [model,inlierIndices,outlierIndices] = pcfitplane(ptCloudIn,maxDistance);
    inliers_2d=normalized2Dpoint(:,inlierIndices);
    inliers_3d=normalized3Dpoint(:,inlierIndices);
       for k=1:iterations
            ind = randsample(size(inliers_2d,2),3);
            Ps = minimalCameraPose(pextend(normalized2Dpoint(:,ind)),normalized3Dpoint(:,ind));
            if numel(Ps) >0
%                 temp_ones = ones(1,size(normalized2Dpoint,2));
%                 for s =1 :numel(Ps)
%                     local_u{s} = [normalized2Dpoint(1:2,:);temp_ones];
%                 end
%                                   %ld implementation
%                   P =Ps;
%                   error = zeros(1,ld_iterations);
%                   intial_camera_matrixes = {P};
%                   local_inliers_3d = [normalized3Dpoint;temp_ones];
%                   intial_3d_points = {normalized3Dpoint};
                
                 for r =1:numel(Ps) 
                   for it =1 : ld_iterations   
                    temp_ones = ones(1,size(normalized2Dpoint,2));
               
                    local_u = {[normalized2Dpoint(1:2,:);temp_ones]};
                
                                  %ld implementation
                  P =Ps(r);
                  error = zeros(1,ld_iterations);
                  intial_camera_matrixes = {P};
                  local_inliers_3d = [normalized3Dpoint;temp_ones];
                  intial_3d_points = {inliers_3d};
                    [ err , res ] = ComputeReprojectionError (P ,local_inliers_3d , local_u );
                    [r_local , J ] = LinearizeReprojErr (P ,local_inliers_3d , local_u);
                    C = J'* J + lambda * speye ( size (J ,2));
                    c = J'* r_local ;
                    deltav = -C \ c ;
                    [ P , local_inliers_3d ] = update_solution ( deltav ,P , local_inliers_3d );   
                    intial_camera_matrixes{it+1} = P;
                    intial_3d_points{it+1} = local_inliers_3d;
                    error(it) = err;
                  end
  
                  %%% end of implementation
                     local_camera =cell2mat(Ps(r));
                     num_inlier = 0;
                    for h =1:size(normalized2Dpoint,2)
                     % local_3d_point = [inliers_3d(:,h);1];
                        local_3d_point = local_inliers_3d(:,h);
                       local_2d_point = normalized2Dpoint(:,h);
                       computed_2d_point = pflat(local_camera * local_3d_point);
                       computed_2d_point = computed_2d_point(1:2,:);
                       resids = computed_2d_point - local_2d_point;
                       errs = sqrt(resids(1,:).^2 + resids(2,:).^2);
                       if errs >0.01
                       num_inlier = 1+num_inlier; 
                       end
                    end
                    if num_inlier > best_inlier
                       best_inlier = num_inlier;
                       best_camera = local_camera;                        
                    end            
                 end
            end 
       end
       P_est{i,j}= best_camera;
       bounding_boxes{i,j}=data{j}.bounding_boxes{i};
       poses{i,j}= data{j}.poses{i};      
    end
    
end
%evaluation
for i =1:number_of_images
    draw_bounding_boxes (images{i}, poses(:,i) , P_est(:,i), bounding_boxes(:,i) );
    scores = eval_pose_estimates (poses(:,i),P_est(:,i),bounding_boxes(:,i));    
end
