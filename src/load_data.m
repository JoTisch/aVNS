%% Load data 
function [y, Resp, ECG, BP, Stim, RR] = load_data(Patient_ID)

id1 = ['./Data/ID1';'./Data/ID2';'./Data/ID3';'./Data/ID4';'./Data/ID5'];
cd(id1(Patient_ID,:))
y = load("protocol_3_2.mat"); 
Resp=y.data(:,1);
% ECG data of patient one is inverted 
if Patient_ID == 1 
ECG =  -y.data(:,2); 
else 
ECG = y.data(:,2); 
end 
%
BP= y.data(:,3); 
Stim= y.data(:,4);
RR= y.data(:,5); 
