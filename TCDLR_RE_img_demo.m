addpath(genpath(cd))
clear all;
clc;

method={'TCDLR_RE'};
short_name={'TCDLR-RE'};
load small_img
I=data;

option.lp=0.8;
option.rate=0.5;

p=0.3; % sampling rate
j=1; % methods
PSNR_=zeros(1,50);
TIME_=zeros(1,50);
for i=1:length(I)% data
        X=double(I{i});
        normalize=max(X(:));
        X=X/normalize;
        funcc=str2func (method{j});
        [Xhat,time,psnr_result,ssim_result,fsim_result] = funcc(X,p,i,option);
        PSNR_(i)=psnr_result;
        TIME_(i)=time;
end


