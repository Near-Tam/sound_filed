%Բ�λ�������ѹ�ֲ�
function [ P ] = Sound_pressure_circular( k,Rs,Sita, w,p0,c,r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    X=k*Rs*sin(Sita);
    J1=besselj(1,X);
    %�ɽ��ۣ�ֱ�Ӵ������
    P=(w*p0*Rs^2./(2*c*r)).*(2*J1./X);
end

