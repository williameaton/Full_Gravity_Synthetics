clc
clear all 
close all 

% Import MT file produced in Python script plot_source.py
% This has columns [lat, lon, M0, delay_time]
Data = importdata('Tohoku_CMT_for_MATLAB.txt');

norm = 2e24;


geoscatter(Data(1,:),Data(2,:), Data(3,:)/norm,Data(4,:), 'o', 'filled')
title('Tohoku 2011 Moment Tensor progression Ji Chen ')

geobasemap colorterrain
