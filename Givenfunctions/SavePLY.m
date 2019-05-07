% This code is for saving the reconstructed 3D points into PLY file.
% You can open PLY files using 3D viewer like 'MeshLab'.

function SavePLY(filename, X)
    out=fopen(filename,'w');
    fprintf(out,'ply\n');
    fprintf(out,'format ascii 1.0\n');
    fprintf(out,'element vertex %d\n',size(X,2));
    fprintf(out,'property float x\n');
    fprintf(out,'property float y\n');
    fprintf(out,'property float z\n');
    fprintf(out,'property uchar diffuse_red\n');
    fprintf(out,'property uchar diffuse_green\n');
    fprintf(out,'property uchar diffuse_blue\n');
    fprintf(out,'end_header\n');
    for i=1:size(X,2)
        fprintf(out,'%f %f %f %d %d %d\n',[X(1,i),X(2,i),X(3,i),min(round(X(4,i)*255),255),min(round(X(5,i)*255),255),min(round(X(6,i)*255),255)]);
    end
    fclose(out);
end