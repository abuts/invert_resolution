function [t_seq,dt,Np] = adjust_step(t_min,t_max,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t_seq = t_min:dt:t_max;
Np = numel(t_seq);
dt = (t_max-t_min)/(Np-1);
t_seq = t_min:dt:t_max;
end

