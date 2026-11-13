%% METHODS Pipeline
clear all; clc;

%% 0 load variables (structural connectome + functional timecourses)
% compute & decompose Graph Laplacian
mypath=findpath; %set code folder path
% Load preprocessed and z-scored fMRI time series
ts = load('fmri_data.mat', 'subj_timeseries');
% Load volume-normalized individual SC matrices
sc = load('SC_data.mat', 'SC_norm');
% Individual connectome harmonic decomposition
GSP_Laplacian
%% 1 analyze graph signal
GSanalysis

%% 2 Energy
% 2.1 
Event_energy
% roiResults
% netResults
% brainRsults
%% 2.2 
% 3.2.1 roi
Energy_ROI_Statistics
% 3.2.2 network 
Energy_Network_Statistics

%%  NDR
clc;clear;
Event_NRD
LME_NRD
%% PLS
PLS_NRD