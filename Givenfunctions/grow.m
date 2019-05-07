function [Pcell,x2dX3d,Xp,colorp] = grow(pair,M,F,K,MF,x2dX3d,Pcell,images,max_try,t,colorold,grownum)

% assigning old pair of images
old = pair;
Xp = [x2dX3d{1,pair(1)}(4:6,:)];

% color information for old images
colorp = [colorold];

    for num = 1: length(MF)-grownum
        % Finding a new image and its pair for the growing step
        [new,pairIm] =newIm(M,old);
        [xd2,Xd3,xd2_check,ref_feat]=corr2d3d([new,pairIm],MF,F,K,x2dX3d,Pcell);
        Xd3 = [Xd3;ones(1,size(Xd3,2))];
        if(size(Xd3,2)<3)
            break;
        end
        disp('Growing Image Number: ');disp(num)
        xd2Xd3 = [xd2;Xd3]';
        
        x2dX3d{1,new} = [xd2Xd3(:,1:6)' ;ref_feat];
        ninliers = 0;     % Number of inliers
        
        % Starting P3P and RANSAC
        for i = 1: max_try

            ind = randsample(size(xd2Xd3,1), 3);

            data_p3p = xd2Xd3(ind,:);

            RT = PerspectiveThreePoint(data_p3p(:,1:6));

            nP = size(RT,1)/4;
            for j = 1:nP
                Pcheck = RT(1+4*(j-1):4+4*(j-1),:);
                Pcheck = Pcheck(1:3,:);
                % Now implementing the distance matrix as indicated in
                % literature 
                KP1X = K*Pcheck*Xd3;
                KP1X = KP1X./KP1X(3,:);
                
                KP2X = K*(Pcell{1,pairIm})*Xd3;
                KP2X = KP2X./KP2X(3,:);
                     
                temp1 = (K*xd2 - KP1X).^2;
                temp2 = (K*xd2_check - KP2X).^2;
                
                temp1 = temp1(1,:)+temp1(2,:);
                temp2 = temp2(1,:)+temp2(2,:);
                
                dsq = temp1 +temp2;
                d = sqrt(dsq);
                
                inliers = find(abs(d) < t);     % Indices of inlying points

                if length(inliers) > ninliers   % Record best solution
                    ninliers = length(inliers);
                    bestP = Pcheck;
                    bestInliers = inliers;
                end
            end    
        end
        
        % Storing the best P yet in cell structure
        Pcell{1,new} = bestP;

        % Updating the list of old images as the current new image is not
        % used for growing
        old = [old new];
        
        % Now starting thetriangulation procedure to to append to existing
        % point thus growing the view effectively
        
        CP1 = cell2mat(Pcell(1,pairIm)); % camera matrix

        CP2 = bestP; % Best camera matrix from P3P


        % All the matching features in the pair
        ref_feat = MF{pairIm,new}; % All features of pair 1 and pair 2

        % x and y co-ord of all the normalized features in image 1 and imag2 
        f_norm1 = inv(K)*[F{pairIm}(1:2,:); ones(1, size(F{pairIm}(1:2,:),2))];
        f_norm2 = inv(K)*[F{new}(1:2,:); ones(1, size(F{new}(1:2,:),2))];

        f_unnorm1 = [F{pairIm}(1:2,:); ones(1, size(F{pairIm}(1:2,:),2))];
        Q_unnorm = round(f_unnorm1(:,ref_feat(1,:)));
        colornew = zeros(3,length(Q_unnorm));  

        for i = 1 : length(Q_unnorm)

        color = images{pairIm};%(Q_unnorm(2,:),Q_unnorm(1,:),:);
        colornew(:,i) = [color(Q_unnorm(2,i),Q_unnorm(1,i),1), color(Q_unnorm(2,i),Q_unnorm(1,i),2) ,color(Q_unnorm(2,i),Q_unnorm(1,i),3) ];
        end

        % x and y cord of matching features  
        Q1_feat = f_norm1(:,ref_feat(1,:));
        Q2_feat = f_norm2(:,ref_feat(2,:));


        Q_feat = [Q1_feat;Q2_feat];

        Q1 = Q_feat(1:3,:);
        Q2 = Q_feat(4:6,:);
        X1 = []; 
        % Traingulation step
        for i=1:size(Q_feat,2)
            A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
            [U,S,V] = svd(A);
            X1(:,i) = V(:,4);
        end
        X1 = X1./X1(4,:);
        
        X1 = X1(1:3,:);
        X_exist = X1(:,find(X1(3,:)>=0));
        Q_exist = Q_feat(:,find(X1(3,:)>=0));
        ref_exist = ref_feat(:,find(X1(3,:)>=0));
        colornew = colornew(:,find(X1(3,:)>=0));
        
        % Storing a new set of features so that growing step doesnot starve
        % for points
        x2dX3d{1,new} =[Q_exist(4:6,:);X_exist;ref_exist(2,:)];
       
        
        % This the duplicate point removal step 
        test1 = round(X_exist(1:3,:)*5);
        [a,indx,c]= unique(test1','rows');
        X_exist = X_exist(1:3,indx);
        colornew = colornew(:,indx);
        Xp = [Xp,X_exist(1:3,:)];
        colorp = [colorp colornew];
    end
end


