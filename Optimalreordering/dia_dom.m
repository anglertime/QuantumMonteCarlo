function [ cc ] = dia_dom( B )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = size(B, 1);
for i = 1:N
    d = abs(B(i, i));
    su = sum(abs(B(i, :))) - d;
    cc(i) = d/su;
end
% plot(cc, '--*b');

