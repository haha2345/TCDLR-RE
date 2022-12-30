%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluate background initialization results as compared to ground thruth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Changes:
% 2015-04-14 - SSIM substituted by MS-SSIM
% 2016-01-25 - added 7 more image sequences 
%            - eliminated (redundant) EPs and PEPs metrics
%
clear all
close all
%
% Adapt here your parameters:
%
%   NameDirGT   = directory containing reference background images
%   NameDirBC   = directory containing computed background images
%   Format      = type of input images
%   Method      = background initialization method to evaluate. E.g., 'Mean'
%   OutputStyle = format for output table of results. Available: 
%                 'latex' or anything else
%   Show        = to display images (GT, BC, |GT-BC|, EPs, CEPs)

NameDirGT   = '../GT/';
Method      = 'Mean';
NameDirBC   = strcat('../Results/',Method,'/');
Format      = '.png';
OutputStyle ='excel';
Show        = 0;

%
% Sequence names (do not need to be changed, if used with provided sequences)
%
Sequence = struct('Name', {});
Sequence(1).Name  = 'Board';
Sequence(2).Name  = 'Candela_m1.10';
Sequence(3).Name  = 'CAVIAR1';
Sequence(4).Name  = 'CAVIAR2';
Sequence(5).Name  = 'CaVignal';
Sequence(6).Name  = 'Foliage';          
Sequence(7).Name  = 'HallAndMonitor'; 
Sequence(8).Name  = 'HighwayI';
Sequence(9).Name  = 'HighwayII';
Sequence(10).Name = 'HumanBody2';
Sequence(11).Name = 'IBMtest2';
Sequence(12).Name = 'PeopleAndFoliage'; 
Sequence(13).Name = 'Snellen';          
Sequence(14).Name = 'Toscana';          
numSeq=length(Sequence);

AGE=zeros(numSeq);
pEPs=zeros(numSeq);
pCEPs=zeros(numSeq);
MSSSIM=zeros(numSeq);
PSNR=zeros(numSeq);
CQM=zeros(numSeq);
%
%Evaluate performance measures
%
for i=1:numSeq
    NameGT=strcat(NameDirGT,'GT_',Sequence(i).Name,Format);
    NameBC=strcat(NameDirBC,'BC_',Sequence(i).Name,Format);

    [AGE(i),pEPs(i),pCEPs(i),MSSSIM(i),PSNR(i),CQM(i)]=Evaluate3(NameGT, NameBC, Show);
end

%
%Print performance results
%

if (strcmp(OutputStyle,'latex'))
	OutputFormat='%16s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\n';
else
	OutputFormat='%16s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n';
end

fprintf('**************** Results of Method: %s ******************\n',Method);

disp('Sequence            AGE       pEPs%     pCEPS%    MSSSIM    PSNR       CQM');
for i=1:numSeq
    fprintf(OutputFormat,Sequence(i).Name,AGE(i),pEPs(i)*100,pCEPs(i)*100,MSSSIM(i),PSNR(i),CQM(i));
end
