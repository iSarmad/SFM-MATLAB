function [MF,M,F,D] = matchFeat(images,peak_thres,edge_thres,foctave)
N = length(images); % Number of Picture 
D{N,1} = []; % Descriptor Cell 
F{N,1} = []; % Feature Cell
MF{N,N} = []; % Best Match Cell
M = zeros(N,N);


% Extract All features and store in cell
Scores{N,N} = [];
    for i = 1 : N
        im = single(rgb2gray(images{i}));
        [f,d] = vl_sift(im,'PeakThresh',peak_thres,'EdgeThresh',edge_thres,'FirstOctave',foctave);
        F{i}=f;        
        D{i}=d;
    end
% Match All features and store in cell
    for i = 1 : N
        for j = 1 : N
            if i~= j
            [MF{i,j}, Scores{i,j}] = vl_ubcmatch(D{i}, D{j});
            end
         end
    end
% Store number of matched features in a matrix for ease
    for i = 1 : N
        for j = 1 : N
            M(i,j) =length(MF{i,j}) ;
            if i==j
                M(i,j) = 0;
            end
         end
    end
    
    
end