clear all
close all
clc
%%

addpath(genpath(strcat(pwd,'\Functions\'))); %Please add to the path or in this folder 

%--------------------Available datasets-----------------------
%   sardania_dataset
%   omodeo_dataset
%   alaska_dataset
%   Madeirinha_dataset
%   katios_dataset
%   dique_dataset
%   SF_dataset
%   Wenchuan_dataset
%   toulouse_dataset
%   canada_dataset
%   california_flood
%   contest_dataset
%   Bastrop_dataset
%   gloucester_dataset
%-------------------------------------------------------------

dataset = "sardania_dataset";

%if you want to try your own dataset it must be like:
%   dataset{1} = before; image with one chanel
%   dataset{2} = after;  image with one chanel
%% Parameters

n_samples = 96;
Ksmooth = 10;
stop = 0.001;

%% Detection.

Change_map = bs_smooth_cd(dataset, n_samples, Ksmooth, stop);

%% Metrics and error map

%all available datasets has the gt
load(dataset,'gt');

[MA, FA, Precision, recall, kappa, OE] = cohensKappa(gt,Change_map);