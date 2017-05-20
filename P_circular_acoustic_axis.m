%Ô²ĞÎ»îÈûÖáÏßÉùÑ¹
function [ P ] = P_circular_acoustic_axis(lambda,Rs,r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
P=2*sin(pi*((Rs^2+r.^2).^0.5-r)/lambda);
end

