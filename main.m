close all;
clear all;
addpath('vlfeat-0.9.21/toolbox')
vl_setup
addpath('Givenfunctions');
%% Attention TA! Please set these variables first. Thank you 
dataSet =1; % Select 1 for your provided dataset, select 2 for my own dataset.
highRes = 0; % Select 1 for high Resolution (20 to 30 minutes), 0 for low resolution (3 to 4 minutes)
optim = 0; % set to 1 to perform optimization step for two view as demonstration.
%% Define constants and parameters
if dataSet ==1

K = [ 1698.873755 0.000000     971.7497705;
      0.000000    1698.8796645 647.7488275;
      0.000000    0.000000     1.000000 ];
 
else

K = [ 2291.56576 0.000000     1019.93238;
      0.000000    2296.52628 1022.24806;
      0.000000    0.000000     1.000000 ];
end
  
%% Feature extraction and matching
% Load images and extract features and find correspondences.
disp('Loading Dataset')
if dataSet ==1
 type = '*.JPG';
  path = 'Data/';
else
    type = '*.jpg';
    path = 'Data/self/data/';
end

% Load Images 
images= loadImages(path,type);

% Finding Features and Matches for all the images
disp('Loading Dataset Complete')
disp('Finding All Feature and Matches ,Please wait :)')
if dataSet == 1  
    if highRes ==1
    [MF,M,F,D] = matchFeat(images,1,10,-1);% Args:  Input Images, Peak Threshold, Edge Threshold, First Octave
    else
    [MF,M,F,D] = matchFeat(images,3,10,0);% Args:  Input Images, Peak Threshold, Edge Threshold, First Octave
    end
else
    if highRes ==1
        [MF,M,F,D] = matchFeat(images,1,10,-1);
    else
        [MF,M,F,D] = matchFeat(images,1,10,0);
    end
end

disp('Plotting Some Sample Feature Maps')
% Plot some random images and overlap features to observe  
k=2;
n1 =5;
imagesc(rgb2gray(images{k}));colormap(gray(255))

im = single(rgb2gray(images{k}));
im2 = single(rgb2gray(images{k+n1}));
hold on
f = F{k};
%f(3,:) = 1.0;
perm = randperm(size(f,2)) ;
sel = perm;
h1 = vl_plotframe(f(1:2,sel)) ;
set(h1,'color','r','linewidth',3) ;
hold off
title('Feature Extraction Step')
figure;imagesc(rgb2gray(images{k+n1}));colormap(gray(255))
f2 = F{k+n1}; 
f2(3,:) = 1.0;
perm = randperm(size(f2,2)) ;
sel = perm;%(1:50) ;
h2 = vl_plotframe(f2(1:3,sel)) ;
set(h2,'color','r','linewidth',3) ;
title('Feature Extraction Step')

%% Initialization step (for 2 views)
disp('Two View Step in Progress....')

[ref,pair] = refIm(M);
disp('Reference Pair Image :');disp(pair)


max_try = 1000; % Maximum Steps for RANSAC
t= 0.01;% Setting Threshold for RANSAC
[E_best,inliers_best,Q_feat,colornew,reffeat] = init(MF,F,ref,pair,K,t*1e-5,max_try,images);
disp('Best Essential Matrix:');disp(E_best)
    
% Camera Matrix
[P,X] = cmat(E_best,Q_feat);

disp('Correct Projection Matrix:');disp(P)

% Storing the Projection matrix for each image in a cell
Pcell{1,length(MF)} = [];
Pcell{1,pair(1)} = [eye(3,3) zeros(3,1)]; 
Pcell{1,pair(2)} = P;


% Making new cell which has all images and their 2d and 3d correspondences 
X = X(1:3,:);
X_exist = X(:,find(X(3,:)>=0));
Q_exist = Q_feat(:,find(X(3,:)>=0));
ref_exist = reffeat(:,find(X(3,:)>=0));

x2dX3d{1,length(MF)} = [];
x2dX3d{1,pair(1)} =[Q_exist(1:3,:);X_exist;ref_exist(1,:)]; 
x2dX3d{1,pair(2)} =[Q_exist(4:6,:);X_exist;ref_exist(2,:)];


% Plotting the final two view with colors
X_exist = X(:,find(X(3,:)>=0));
color_exist = colornew(:,find(X(3,:)>=0));
X_exist = [X_exist; color_exist/255];
figure;pcshow(X_exist(1:3,:)',X_exist(4:6,:)','MarkerSize',20);shg;xlabel('x');ylabel('y');zlabel('z');
if dataSet ==1
axis([-4 4 -4 4 0 20]);shg
else
axis([-4 4 -4 4 0 10]);shg
    
end
title('Final Two View SFM with Color')

% saving the two view file
disp('Saving Ply file for Two View but please use Matlab visualizer')

filename = sprintf('%02dviews.ply', 2);

SavePLY(filename, X_exist);
%% Growing Step 
disp('Starting the Growing Step')

[Pcell,test2,Xp,colorp] = grow(pair,M,F,K,MF,x2dX3d,Pcell,images,10*max_try,t*100,color_exist,2);

Xp1 = Xp(1:3,:);
rot = [-1  0  0;  0 -1  0; 0  0  1]; % Flipping the point cloud for better view

Xp1 = rot * Xp1;
X_existp = Xp1(:,find(Xp1(3,:)>=0));
color_existp = colorp(:,find(Xp1(3,:)>=0));
X_existp = [X_existp; color_existp/255];

if dataSet ==1 % Plotting accordingly depending on dataset
    if highRes ==1
        figure;pcshow(X_existp(1:3,:)',X_existp(4:6,:)','MarkerSize',20);shg;axis([-5 5 -5 5 0 30]);xlabel('x');ylabel('y');zlabel('z');shg
    else
        figure;pcshow(X_existp(1:3,:)',X_existp(4:6,:)','MarkerSize',20);shg;axis([-5 5 -5 5 0 30]);xlabel('x');ylabel('y');zlabel('z');shg
    end
else
    if highRes ==1
    figure;pcshow(X_existp(1:3,:)',X_existp(4:6,:)','MarkerSize',10);shg;axis([-5 5 -5 5 0 30]);xlabel('x');ylabel('y');zlabel('z');shg
    else
    figure;pcshow(X_existp(1:3,:)',X_existp(4:6,:)','MarkerSize',40);shg;axis([-5 5 -5 5 0 30]);xlabel('x');ylabel('y');zlabel('z');shg    
    end
end

title('Final Grown Image');
hold on
for kdx = 1:length(Pcell)
if ~isempty(Pcell{1,kdx})
pcshow(Pcell{1,kdx}(:,4)','MarkerSize',1000)
dx = 0.5; dy = 0.5; dz = 0.5; % displacement so the text does not overlay the data points
c = num2str(kdx);
text(Pcell{1,kdx}(1,4)'+dx, Pcell{1,kdx}(2,4)'+dy, Pcell{1,kdx}(3,4)'+dz ,c);
end
end
if dataSet ==1
axis([-20 20 -20 20 0 30])
else
axis([-4 4 -7 7 0 10])
end
hold off

% saving ply file for the Growing Step
disp('Saving Ply file but please use Matlab visualizer')
filename = 'growviews.ply';
SavePLY(filename, X_existp);

%% Two view optimization step (Using Online repository)

%Please note that this step has been mostly borrowed from the online repository 
% Author: Riccardo Giubilato, Padova (Italy) 2016
% mail:   riccardo.giubilato@gmail.com
% https://www.researchgate.net/profile/Riccardo_Giubilato
%Performing bundle adjustment for two view 
if optim==1
%ind = randsample(size(Q_feat,2), 405);
ob1 = (K*Q_feat(1:3,:))';
%ob1 = ob1(ind,:);
ob1(:,3) = 1:length(ob1); 
ob2 = (K* Q_feat(4:6,:))';
%ob2 = ob2(ind,:);
ob2(:,3) = 1:length(ob2); 
obs1 = {ob1,ob2};
%app = [0 0 0 1];
[Points, Camera] = BA(X(1:3,:)', cat(3,eye(3,4),[P]), obs1, K);%cat(3,eye(4,4),[P;app])

Q1 = Q_feat(1:3,:);
Q2 = Q_feat(4:6,:);
CP1 = Camera(:,:,2);
CP2 = Camera(:,:,1);

for i=1:size(Q_feat,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
end

Xtemp = Xtemp./Xtemp(4,:);
X44 =  Xtemp;
X44 = X44(1:3,:);
X_exist44 = X44(:,find(X44(3,:)>=0));

figure;pcshow(X_exist44(1:3,:)','MarkerSize',20);shg;xlabel('x');ylabel('y');zlabel('z');axis([-4 4 -4 4 0 20]);shg

title('Two View Step After Bundle Adjustment')
end