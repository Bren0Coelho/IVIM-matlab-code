%% SNR_EST
% =======================================================================
% =======================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   11/15/2020: ver 1.0
% ========================================================================
% ========================================================================
% This function estimates the SNR of the image
% ========================================================================
% ========================================================================

clear;
clc;
close;
load('workspace_simulation_RETICULADO_v3_NOISE.mat');
pMode = 'db';
snr = 20;

% Catch information from image
iSNR = find(SNR == snr);
MR_b0 = Voxel_noise(:,:,1,iSNR);                                                         % image for b = 0 s/mm² 
text1 = ['Acquisition for b = ', num2str(b(1)),...
    ' s/mm^{2}, SNR = ', num2str(SNR(iSNR))];     % title of image           
MR_gray = mat2gray(MR_b0);
image(MR_gray, 'CDataMapping','scaled');                                      % image
colormap(gray);
title(text1);
set(gcf, 'Position', get(0, 'Screensize'));
set(gca, 'XTick', [], 'YTick', []);
h1 = imrect; BW1 = createMask(h1);
h2 = imrect; BW2 = createMask(h2);
[s1,q1] = find(BW1);
[s2,q2] = find(BW2);
noise = MR_b0(s1,q1);
sig = MR_b0(s2,q2);

% Signal power and noise power 
sigPower = sum(abs(sig(:)).^2)/length(sig(:));
noisePower = sum(abs(noise(:)).^2)/length(noise(:));

% SNR
if(strcmp(pMode,'db'))
    % Decibel
    sigPower = 10*log10(sigPower);
    noisePower = 10*log10(noisePower);      
    reqSNR = sigPower-noisePower;
else
    reqSNR = sigPower/noisePower;
end

disp(['SNR: ', num2str(reqSNR)]);