function [maxIm,pairIm] = newIm(M,old)
% This function helps find the new pair of images for the growing step
feats = zeros(1,length(old));
maxfeats = 0;
for j = 1 : length(M)
    if ~ismember(j,old)
        for i = 1: length(old)
            feats(i) = M(j,old(i));
        end
        featsS = sum(feats);
        if maxfeats<featsS
            maxIm = j;
            [arg,ind]= max(feats);
            pairIm = old(ind);
        end
        maxfeats = max(maxfeats,featsS);
    end
end
end