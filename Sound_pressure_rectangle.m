%矩形活塞的声压分布
function [ P] = Sound_pressure_rectangle( p0,Fs,lambda,r,k,a,b,sita,phi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
P=(p0*Fs./(lambda.*r)).*(sin(k*a*sin(sita).*cos(phi))./(k*a*sin(sita).*cos(phi))).*(sin(k*b*sin(phi))./(k*b*sin(phi)));

end

