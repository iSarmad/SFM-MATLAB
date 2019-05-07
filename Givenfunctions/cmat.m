function [Cmat,X] = cmat(E_best,Q_feat)

% Perform Singular value decomposition
[u,s,v]=svd(E_best);
diag_110 = [1 0 0; 0 1 0; 0 0 0];
newE = u*diag_110*v';
[u,s,v] = svd(newE);
w = [0 -1 0; 
     1  0 0;
     0  0 1];
 
% Find Four Possible Solutions  
R1 = u*w*v';
R2 = u*w'*v';
t1 = u(:,3);%./max(abs(u(:,3)));%./max(abs(u(:,3)));
t2 = -u(:,3);%./max(abs(u(:,3)));%./max(abs(u(:,3)));

P1 =   [R1 t1];
P2 =   [R1 t2];

P3 =   [R2 t1];
P4 =   [R2 t2];

P = {P1,P2,P3,P4};

Q1 = Q_feat(1:3,:);
Q2 = Q_feat(4:6,:);
CP1 = eye(3,4);
CP2 = P1;

% Perform the traingulation step to show that not all results are correct
for i=1:size(Q_feat,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
end

Xtemp = Xtemp./Xtemp(4,:);
X41 =  Xtemp;
X41 = X41(1:3,:);


CP2 = P2;
for i=1:size(Q_feat,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
end

Xtemp = Xtemp./Xtemp(4,:);
X42 =  Xtemp;
X42 = X42(1:3,:);


CP2 = P3;
for i=1:size(Q_feat,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
end

Xtemp = Xtemp./Xtemp(4,:);
X43 =  Xtemp;
X43 = X43(1:3,:);

CP2 = P4;

for i=1:size(Q_feat,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
end

Xtemp = Xtemp./Xtemp(4,:);
X44 =  Xtemp;
X44 = X44(1:3,:);

X_exist = {X41,X42,X43,X44};

P ={P1,P2,P3,P4};
CP1 = [eye(3,3),zeros(3,1)];

Depth = [];

% Now Checking which P is correct by finding camera depth for both cameras

for j = 1:4
X = X_exist{1,j}; 
figure;pcshow(X');shg;xlabel('x');ylabel('y');zlabel('z');axis([-4 4 -4 4 0 20]);shg
title('Triangulation Using All possible R and T (Only one will be correct)')

CP2 = P{1,j};

x =  det([ CP1(:,2), CP1(:,3), CP1(:,4) ]);
y = -det([ CP1(:,1), CP1(:,3), CP1(:,4) ]);
z =  det([ CP1(:,1), CP1(:,2), CP1(:,4) ]);
t = -det([ CP1(:,1), CP1(:,2), CP1(:,3) ]);
c1 = [ x/t; y/t; z/t ]; % finding camera center

CP2 = P{1,j};

x =  det([ CP2(:,2), CP2(:,3), CP2(:,4) ]);
y = -det([ CP2(:,1), CP2(:,3), CP2(:,4) ]);
z =  det([ CP2(:,1), CP2(:,2), CP2(:,4) ]);
t = -det([ CP2(:,1), CP2(:,2), CP2(:,3) ]);
C2{1,j} = [ x/t; y/t; z/t ];  % finding camera center

% Now checking Which cameras have positive depth 
rot2 = CP2(:,1:3); %taking only the rotation matrix
rot1 = CP1(:,1:3);
c2 = C2{1,j};

    for i=1:size(Q_feat,2)
   % Camera 1
    wc1 = rot1(3,:) * (X(1:3,i) - c1(1:3,:));
    depthc1 = (sign(det(rot1)) * wc1) / 1 * norm(rot1(3,:));
    DC1(i) = depthc1(1,1);
    
   %Camera 2
    wc2 = rot2(3,:) * (X(1:3,i) - c2(1:3,:));
    depthc2 = (sign(det(rot2)) * wc2) / 1 * norm(rot2(3,:)); 
    DC2(i) = depthc2(1,1);
    
    end
    D1 = sum(sign(DC1));
    D2 = sum(sign(DC2));
    
    Depth = [Depth D1+D2];


end
[arg,val] = max(Depth);

X = cell2mat(X_exist(val));
X = [X;ones(1,length(X))];
Cmat = cell2mat(P(val));
end