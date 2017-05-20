function [ P ] = P_rectangle_acoustic_axis( Fs,a,k,lambda,r ,sita,phi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
P=(Fs./(lambda.*r)).*(sin(k*a*sin(sita).*cos(phi))./(k*a*sin(sita).*cos(phi)));
P=abs(P);
end

