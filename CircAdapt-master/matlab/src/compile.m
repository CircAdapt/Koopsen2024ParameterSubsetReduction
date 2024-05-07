%% Example1 matlab
% Nick van Osta - April 2020
% Compile CircAdapt for Matlab

%Set path to C++ files
CircAdaptCpp = '../../core/';

% Set path to CircAdapt Matlab Course
CircAdaptMatlab = './';

% compile cpp, should be done only once or when you change the code
CircAdaptCompile(CircAdaptCpp,CircAdaptMatlab);
