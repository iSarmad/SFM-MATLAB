function [d2,d3,d2_check,ref_feat_match1] = corr2d3d(pair,MF,F,K,x2dX3d,Pcell) %[d2,d3,d2_check]

% All the matching features in the pair
temp1= pair(1);
temp2 = pair(2);
pair(1) =temp2; pair(2) =temp1;
s=1;
ref_feat = MF{pair(1),pair(2)}; % All features of pair 1 and pair 2

% x and y co-ord of all the normalized features in image 1 and imag2 
f_norm1 = inv(K)*[F{pair(1)}(1:2,:); ones(1, size(F{pair(1)}(1:2,:),2))];
f_norm2 = inv(K)*[F{pair(2)}(1:2,:); ones(1, size(F{pair(2)}(1:2,:),2))];

old_feat = x2dX3d{1,pair(s)}(7,:); 

% Comparing old image pair features and new image pair feature to find a subset of the features
% which are common in old pair and new pair
[tf1, idx1] = ismembertol(old_feat,ref_feat(s,:));
ind2 = find(tf1(1,:));
ref_feat_match = old_feat(s,ind2);

% Finding all the subset image point in old image
d2_check = f_norm1(:,ref_feat_match(1,:)); 


% Finding all the subset image point in new image
ref_feat_match1 = ref_feat(2,idx1(idx1>0));
d2 = f_norm2(:,ref_feat_match1(1,:));

% Finding all the subset World point in new image
d3 = x2dX3d{1,pair(s)}(4:6,ind2);

end