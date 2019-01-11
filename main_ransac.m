clear
close all
clc

%% Load and prepare data
curdir = strcat(pwd,'/data/');
image_info = dir(fullfile(curdir,'img*.jpg'));
data_info = dir(fullfile(curdir,'data*.mat'));
images = cell(1,numel(image_info));
data = cell(1,numel(data_info));
for i=1:numel(image_info)
 images{i} = imread(strcat(curdir,image_info(i).name));
 data{i} = load(strcat(curdir,data_info(i).name));
end

for i=1:7
 for j=1:numel(image_info)
     U{i,j}=(data{j}.U{i}-repmat(mean(data{j}.U{i},2),[1 size(data{j}.U{i},2)])); %rows = objects, cols = cameras
     u{i,j}=data{j}.u{i};
     bboxes{i,j}=data{j}.bounding_boxes{i};
     trues{i,j}= data{j}.poses{i};
 end
end

%%
% obj1 = buda
% obj2 = pot
% obj3 = cat
% obj4 = duck
% obj5 = eggbox
% obj6 = bottle
% obj7 = stapler
% l = 5;
%  for j=1:7
%      temp = U{j,l};
%      tempu = u{j,l};
%      tempPose = trues{j,l};
%      tempTranslated = tempPose*[temp;ones(1,length(temp))];
%      reprojected = (([eye(3),zeros(3,1)]*[tempTranslated;ones(1,length(tempTranslated))]));
%      figure
%      subplot(1,3,1)
%      plot3(temp(1,:),temp(2,:),temp(3,:),'.')
%      axis ij equal
%      subplot(1,3,2)
%      plot(tempu(1,:),tempu(2,:),'.')
%      axis ij equal
%      subplot(1,3,3)
%      plot(reprojected(1,:),reprojected(2,:),'.')
%      axis ij equal
%      check{j}=tempu-reprojected(1:2,:);
%  end

%% Mean algorithm try
iterations = 10;

for i=1:7
 for j=1:numel(image_info) 
    scores = zeros(1,iterations);
    candidateCams = cell(1,iterations);
%     goodPoints=good_points(U{i,j});
%      ptCloud = pointCloud(U{i,j}(:,goodPoints)');
    meanX=mean(U{i,j},2);%Computes the mean 3D point 
    Xtilde=(U{i,j}-repmat(meanX,[1 size(U{i,j},2)]))

     ptCloud = pointCloud(Xtilde');
    [~,goodPoints,~] = pcfitplane(ptCloud,0.008);
    imagePoints=u{i,j}(:,goodPoints);
    worldPoints=Xtilde(:,goodPoints);
   for k=1:iterations
            ind = randsample(size(imagePoints,2),3);
            Ps = minimalCameraPose(pextend(u{i,j}(:,ind)),U{i,j}(:,ind));
            if ~isempty(Ps)
                [candidateCams{k},scores(k)] = ransac_solver(imagePoints,worldPoints,Ps);
            end
   end
    [~,ind]=max(scores);
     poses{i,j}=candidateCams{ind};
  end
end

for i=1:numel(image_info)
    draw_bounding_boxes(images{i}, trues(:,i) , poses(:,i), bboxes(:,i) );
    scores = eval_pose_estimates (trues(:,i),poses(:,i),bboxes(:,i)  );
end

clearvars data_info image_info


