clc;
clear all;
close all;
warning off;

files = dir('*.tif');

hole_area_sum= {};
for ii=1:length(files)
% for ii = [5]
filename = files(ii).name;
[hole_area,hole_area_hist]=hole_count_2(filename);
hole_area_sum = {hole_area_sum; hole_area};
end