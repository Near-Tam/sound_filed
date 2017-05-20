function [ P ] = surf_cube( la,lb,lc,x,y,z )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
a=la*[0 1];
b=lb*[0 1];
c=lc*[0 1];
[A1,B1]=meshgrid(a,b);
[A2,C2]=meshgrid(a,c);
[B3,C3]=meshgrid(b,c);
P(1)=surf(A1+x,B1+y,0*A1+z);
hold on
P(2)=surf(A1+x,B1+y,0*A1+lc+z);
P(3)=surf(A2+x,0*B1+y,C2+z);
P(4)=surf(A2+x,0*B1+lb+y,C2+z);
P(5)=surf(0*A2+x,B3+y,C3+z);
P(6)=surf(0*A2+la+x,B3+y,C3+z);
for i=1:6
    set(P(i),'facecolor','k');
    set(P(i),'facealpha',0.3);
end

end

