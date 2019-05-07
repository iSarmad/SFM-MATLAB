function [argmax,pair] = refIm(M)
sumc = sum(M');
 [argvalue, argmax] = max(sumc(:));
%[I_row, I_col] = ind2sub(size(M),argmax);
 [a , b] =  max(M(argmax,:));
 pair = [argmax b];
end