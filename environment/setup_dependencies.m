%% CVX Install
%
addpath(genpath('/CVX'));
cvx_setup();
disp('Installed CVX');

%% CORA Install
%
addpath(genpath('/CORA'));
disp('Installed CORA');

%% MPT3 Install
%
addpath(genpath('/MPT3'));
mpt_init();
disp('Installed MPT3');

%% YALMIP Install
%
addpath(genpath('/YALMIP'));
disp('Installed YALMIP');

%% SReachTools Install
%
addpath(genpath('/SReachTools'));
disp('Installed SReachTools');


savepath; % For compatibility with CodeOcean.
