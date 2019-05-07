function [bestE,bestInliers,Q_feat,colornew,ref_feat] = init(MF,F,ref,pair,K,t,max_try,images)
% All the matching features in the pair
ref_feat = MF{pair(1),pair(2)}; % All features of pair 1 and pair 2

% x and y co-ord of all the normalized features in image 1 and imag2 
f_norm1 = inv(K)*[F{pair(1)}(1:2,:); ones(1, size(F{pair(1)}(1:2,:),2))];
f_norm2 = inv(K)*[F{pair(2)}(1:2,:); ones(1, size(F{pair(2)}(1:2,:),2))];


% Extracting feature one more time for color extraction
f_unnorm1 = [F{pair(1)}(1:2,:); ones(1, size(F{pair(1)}(1:2,:),2))];
Q_unnorm = round(f_unnorm1(:,ref_feat(1,:)));
colornew = zeros(3,length(Q_unnorm));  
for i = 1 : length(Q_unnorm)
    
color = images{pair(1)};%(Q_unnorm(2,:),Q_unnorm(1,:),:);
colornew(:,i) = [color(Q_unnorm(2,i),Q_unnorm(1,i),1), color(Q_unnorm(2,i),Q_unnorm(1,i),2) ,color(Q_unnorm(2,i),Q_unnorm(1,i),3) ];
end

% x and y cord of matching features  
Q1_feat = f_norm1(:,ref_feat(1,:));
Q2_feat = f_norm2(:,ref_feat(2,:));

% Plotting Matched features on The Pair of Images for visualization
figure;imagesc(rgb2gray(images{pair(1)}));colormap(gray(255))
im = single(rgb2gray(images{pair(1)}));
im2 = single(rgb2gray(images{pair(2)}));
hold on
f =  K*Q1_feat;%K*f_norm1;%
%f(3,:) = 1.0;
perm = randperm(size(f,2)) ;
sel = perm;%(1:50) ;
h1 = vl_plotframe(f(1:2,sel)) ;
set(h1,'color','r','linewidth',3) ;
hold off
title('Feature Matching of Best Image Pair')

figure;imagesc(rgb2gray(images{pair(2)}));colormap(gray(255))
f2 = K*Q2_feat; %K*f_norm2;%
f2(3,:) = 1.0;
perm = randperm(size(f2,2)) ;
sel = perm;%(1:50) ;
h2 = vl_plotframe(f2(1:3,sel)) ;
set(h2,'color','r','linewidth',3) ;
title('Feature Matching of Best Image Pair')


Q_feat = [Q1_feat;Q2_feat];
ninliers = 0;     % Number of inliers

% Starting code for RANSAC and Calibrated five point
    for i = 1: max_try

        % Sampling to get 5 points repeatedly
        ind = randsample(size(Q_feat,2), 5);

        Q1 = Q_feat(1:3,ind);
        Q2 = Q_feat(4:6,ind);

        Evec = calibrated_fivepoint(Q1,Q2);

        nsol = size(Evec,2);
        Emat = permute(reshape(Evec, 3, 3, nsol), [2,1,3]); % Reshaping to get each solution
        E = mat2cell(Emat, 3, 3, ones(1, nsol));

        x1 = Q1_feat;    % Extract x1 and x2 from x
        x2 = Q2_feat;

        nE = length(E);   % Number of solutions to test
        
        for k = 1:nE
          Ex1 = E{k}*x1;
          Etx2 = E{k}'*x2;     

          x2tEx1 = sum(x2.*Ex1);

          d =  x2tEx1.^2 ./ ...
              (Ex1(1,:).^2 + Ex1(2,:).^2 + Etx2(1,:).^2 + Etx2(2,:).^2);

          inliers = find(abs(d) < t);     % Checking all inliers

          if length(inliers) > ninliers   % Storing the best solution till now
            ninliers = length(inliers);
            bestE = E{k};
            bestInliers = inliers;
          end
        end
    end
end