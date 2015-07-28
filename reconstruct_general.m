% Post processing on averaged images

% clc
% close all
clear

addpath 'scripts/'
addpath 'misc/'
addpath 'utils/'
addpath 'utils/export_fig/';
addpath 'solvers/lsqrSOL' 
addpath(genpath('solvers/TVAL3_beta2.4'))
addpath(genpath('solvers/unlocbox-mat'));
addpath 'solvers/spgl1-master/';
addpath(genpath('solvers/NESTA_v1.1/'));
addpath(genpath('AliTestingReg'));

% addpath 'misc\export_fig\';

saveFigs = false;

% calibration sequence parameters
calibrationSequence = 'hadamardSeparable'; % hadamardSeparable
screenPatchSize = 128; 
maskType = '255M030C';
switch maskType
    case {'16Mseq','08Mseq','031M100','127M060','255M030','255M030C','127M060C','255M010C'}
        op_makeSeparable = @(z) subtract_mean(z);
    case {'16Sseq'}
        op_makeSeparable = @(z) z;
end 
%
downSamplingFactor = 4;
downSamplingMode = 'average';
% calibMode = 'pm';

solver = 'lsqr'; 
testMode = 'testImages';


msg1 = sprintf(['reconstruct_general separable (color): calibration-%s, maskType=%s, \n',...
        'screenPatchSize=%d, downSampling=%02d (%s), solver-%s, testMode-%s'],...
    calibrationSequence, maskType, screenPatchSize, downSamplingFactor,downSamplingMode,solver, testMode);
fprintf([msg1,'\n']);

%%
calibDir = sprintf('Data/mask%s_%s%d_14in/',maskType,calibrationSequence,screenPatchSize);
calibCompute_mat = sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode);
if ~exist(calibCompute_mat)
    script_computeCalibMatrices    
else
    fprintf('using calibration data recorded earlier...\n');
end


%% setup inverse problem
script_inverseProblemSetup;

%% 

switch testMode
    case 'simulatedImages';
        script_reconstructSimulatedImages;
    case 'testImages'
        testDir = sprintf('Data/mask%s_%dx%d_testImages_14in/',maskType,screenPatchSize,screenPatchSize);
        script_reconstructTestImages;        
    case 'realImages'
        testDir = sprintf('Data/mask%s_realImages/',maskType);
        script_reconstructRealImages;        
end

%% BM3D or denoising... 

fprintf('done! \n');
