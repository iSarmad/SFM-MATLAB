function ProjectionMatrix=PerspectiveThreePoint(Data)

% Input : 3 x 6 matrix
% [ nx1, ny1, 1, X1, Y1, Z1 ;
%    nx2, ny2, 1, X2, Y2, Z2 ;
%    nx3, ny3, 1, X3, Y3, Z3 ]
% 
% nx, ny : normalized coordinate
% X, Y, Z : 3D space
% 
% Output : 4n * 4 matrix
% [         R1          t1 ;
%     0,     0,    0,    1 ;
%           R2          t2 ;
%     0,     0,    0,    1 ;
%               ...
%           Rn          tn ;
%     0,     0,     0,   1 ];
% that contains " n   RT matrices "
% 
% 
% Example : 
% data=zeros(3, 6);
% data(1,:)=[nx1' X1'];
% data(2,:)=[nx2' X2'];
% data(3,:)=[nx3' X3'];
% 
% RT=PerspectiveThreePoint(data);

a=norm(Data(2,4:6)-Data(3,4:6)); a2=a*a;
b=norm(Data(1,4:6)-Data(3,4:6)); b2=b*b;
c=norm(Data(1,4:6)-Data(2,4:6)); c2=c*c;
t1=1/norm(Data(1,1:3));
t2=1/norm(Data(2,1:3));
t3=1/norm(Data(3,1:3));
cosa=Data(2,1:3)*Data(3,1:3)'*t2*t3;
cosb=Data(1,1:3)*Data(3,1:3)'*t1*t3;
cosc=Data(1,1:3)*Data(2,1:3)'*t1*t2;
cos2a=cosa*cosa;
cos2b=cosb*cosb;
cos2c=cosc*cosc;
ab=a2/b2;
cb=c2/b2;
acb1=ab-cb;
acb2=ab+cb;

Coeff4=1/((acb1-1)*(acb1-1)-4*cb*cos2a);
Coeff3=4*(acb1*(1-acb1)*cosb-(1-acb2)*cosa*cosc+2*cb*cos2a*cosb);
Coeff2=2*(acb1*acb1*(1+2*cos2b)-1+2*(1-cb)*cos2a-4*acb2*cosa*cosb*cosc+2*(1-ab)*cos2c);
Coeff1=4*(-acb1*(1+acb1)*cosb-(1-acb2)*cosa*cosc+2*ab*cos2c*cosb);
Coeff0=(1+acb1)*(1+acb1)-4*ab*cos2c;

A=zeros(4,4);
A(1,1)=-Coeff3*Coeff4;
A(1,2)=-Coeff2*Coeff4;
A(1,3)=-Coeff1*Coeff4;
A(1,4)=-Coeff0*Coeff4;
A(2,1)=1;
A(3,2)=1;
A(4,3)=1;
%-----------------modify---------------------
for i=1:numel(A)
    if isnan(A(i)) || isinf(A(i))
        ProjectionMatrix = -1;
        return;
    end
end
A(isnan(A)) = 0;

%---------------------------------------------
[V,D]=eig(A);
ProjectionMatrix=[];

for n=1:4
    v=D(n,n);
    if imag(v)==0
        u=2*(cosc-v*cosa);
        if abs(u)>1.0e-12
            u=((-1+acb1)*v*v-2*acb1*cosb*v+1+acb1)/u;
            if u>0 & v>0
                s1=a/sqrt(u*u+v*v-2*u*v*cosa);
                s2=u*s1;
                s3=v*s1;
                P=Data(:,4:6);
                Q=Data(:,1:3);
                Q(1,:)=Q(1,:)*s1*t1;
                Q(2,:)=Q(2,:)*s2*t2;
                Q(3,:)=Q(3,:)*s3*t3;
                v1=P(1,:)-P(2,:); v1=v1/norm(v1);
                v2=cross(v1,P(1,:)-P(3,:)); v2=v2/norm(v2);
                v3=cross(v1,v2);
                v4=Q(1,:)-Q(2,:); v4=v4/norm(v4);
                v5=cross(v4,Q(1,:)-Q(3,:)); v5=v5/norm(v5);
                v6=cross(v4,v5);
                R=[v4;v5;v6]'*[v1;v2;v3];
                T=Q(1,:)'-R*P(1,:)';
                ProjectionMatrix=[ProjectionMatrix;R,T;0,0,0,1];
            end
        end
    end
end
