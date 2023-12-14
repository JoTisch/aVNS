%TODO: Revisit the Code, Simplify things,
%TODO: Add description to functions, improve var names
%TODO: Improve inspiration/expiration detection 
%TODO: Add title to plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     choose patient and protocoll                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select Patient_ID between 1-5

clear all; close all; clc; 

patient_id=3;  
addpath './Data'; 
addpath './src'; 
addpath './src/analyse_ecg'; 

[raw, resp, ecg, bp, stim, rr] = load_data(patient_id); 
fs=1000; %Sample rate signal in Hz
l=length(ecg); 
t_signal=((0:l-1)/fs)';%time vector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        DETECT R PEAKS & STIM                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rpeaks,~,~]=analyze_ecg_offline_r(ecg, fs);  
win_stim=window_function(stim,0.1,fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   Categorize stimulation protocoll                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectrogram to visualize different stimulation protocol (frequency-time
% domain) 
nfields=10;
nffts=fs;
[s,~,~]=spectrogram(stim,nffts,[],nfields,fs,'yaxis');
%further progress with the lowest frequencies (s(1,:)) 
win_ptcl=window_function(s(1,:),2,fs);
win_ptcl.beg=(win_ptcl.beg/2) * fs; 
edges=[1,win_ptcl.beg,length(stim)]; 
  
% Label the data 
stimstr=["Non-sync aVNS", "Systole-sync aVNS", "Diastole-sync aVNS"]; 
sflag(edges(1):max(edges),1)=1; 
sflag(edges(2):edges(3),1)=2;
sflag(edges(4):edges(5),1)=3; 
% the first patient has a longer acqusition time, with an additional
% systolic-syn VNS  
if patient_id==1 
sflag(edges(6):edges(7),1)=2; 
end 
stim_flag=sflag(win_stim.beg);

[stimRR,stim_alpha,stim_tR,stim_flag,win_stim.beg]=initRR(win_stim.beg, ...
                                                   rpeaks, ...
                                                   t_signal, ...
                                                   stim_flag); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       Normalize ECG and histogram                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
nflag=zeros(length(sflag),1); 
nflag(win_stim.beg)=1; 
x1=linspace(0,1,length(rpeaks)); 
%
for j=1:length(rpeaks)-1
% split ECG signal and flag from R to R peak 
segment= rpeaks(j):rpeaks(j+1); 
split_ecg{j,:}=[ecg(segment),nflag(segment),sflag(segment)];
%find stimulation point in splitted ECG signal 
split_ecg_sloc{j,1}=find(split_ecg{j,1}(:,2) == 1,1,'last'); 
label{j,1} = split_ecg{j,1}(split_ecg_sloc{j,1},3); 
% calulate timepoint of stimulation in degrees
beta{j,1}=split_ecg_sloc{j,1}  / length(split_ecg{j,:}(:,1)); 
% Interpolate = Normalization of ECG signals
% temporal normalization 
split_ecg_res(j,:)=interp1(linspace(0,1,length(split_ecg{j,1}(:,1))), ...
    split_ecg{j,:}(:,1), x1);
% amplitude normalization 
split_ecg_res(j,:)=rescale(split_ecg_res(j,:),0,1); 
end 

colMeans=mean(split_ecg_res,1);
prc=prctile(split_ecg_res,[25 75],1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                Classification Inspiration / Expiration                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Moving Average Filter
windowSize=1.5*fs; 
b=(1/windowSize)*ones(1,windowSize);
a=1; 
maf=filtfilt(b,a,resp);
dmaf=gradient(maf); 
% window function inspiration/ expiration 
for i=1:length(maf)
    if dmaf(i)<=0 
        wresp(i)=-1;
    elseif dmaf(i)> 0
        wresp(i)=1; 
    else 
        wresp(i)=wresp(i-1); 
    end 
end 
[~,seg.begIn]=findpeaks(wresp); 
[~,seg.begEx]=findpeaks(-wresp); 

% Ensure same number of Insp/ Exp. 
if length(seg.begEx)>length(seg.begIn)
    seg.begEx(end)=[]; 
elseif length(seg.begIn)>length(seg.begEx)
    seg.begIn(end)=[]; 
end 
% Assign label 
resp_flag(1:length(seg.begEx),1)="Expiration"; 
resp_flag(length(seg.begEx)+1:2*length(seg.begEx),1)="Inspiration"; 
[start_InEx,I]=sort([seg.begEx, seg.begIn]); 
resp_flag=resp_flag(I); 

[respRR,resp_alpha,resp_tR,resp_flag,start_InEx]=initRR(start_InEx,...
                                                rpeaks,...
                                                t_signal,...
                                                resp_flag); 
% classify stimulation
if seg.begEx(1)<seg.begIn(1)
stack =[seg.begEx;seg.begIn]';
else 
stack=[seg.begIn;seg.begEx]';
end
 
for kk = 1:length(win_stim.beg)
    if find(win_stim.beg(kk) > stack(:,1) & win_stim.beg(kk) < stack(:,2))
        class(kk) = "Expiration"; 
    else 
        class(kk) = "Inspiration"; 
    end 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   Filter signal with mean resp frequency                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanF=meanfreq(resp,fs);
q=3;%Hz
[mf_respRR,~]=RR_filt(resp_tR,respRR,meanF,q); 
[mf_stimRR,~]=RR_filt(stim_tR,stimRR,meanF,q); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           Adaptive Filter                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
af_respRR= adapt_filt(seg.begIn,resp_tR,t_signal,respRR,meanF,q); 
af_stimRR= adapt_filt(seg.begIn,stim_tR,t_signal,stimRR,meanF,q); 

% SAVE for further investigations 
save(['./src/out/Patient',num2str(patient_id),'.mat'],"stim_alpha","stimRR","mf_stimRR","af_stimRR","stim_flag")

