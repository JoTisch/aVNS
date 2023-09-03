function result = window_function(signal,threshold,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Description-
%Get window function of alternating signal, incl. start & end 
%
%-Input-
%signal: raw signal 
%threshold: threshold for window
%fs: sample frequency 
%
%-Output-
%Result: struct
%       - .all: window function of stimulation 
%       - .beg: beg of stimulation pulse 
%       - .end: end of stimulation pulse 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initialize window function 
avg_stim = mean(signal);
std_stim = std(signal);
norm_stim = signal - avg_stim;
ratio = round(norm_stim./std_stim);
win_stim = zeros(1,length(signal));
win_stim(ratio> threshold | ratio< -threshold ) = 1;
% when signal changes between pos and neg a timestep is not detected as
% stimulus. So overcome this problem it is manually asigned to signal (=1). 
for i=1:length(win_stim)-10
    if win_stim(i+1) == 0 && win_stim(i+2) ==1 && win_stim(i) == 1 || ...
       win_stim(i+1) == 0 && win_stim(i+5) ==1 && win_stim(i) == 1 || ...
       win_stim(i+1) == 0 && win_stim(i+10) ==1 && win_stim(i) == 1 
       win_stim(i+1)=1;  
    end
end 
%now findpeaks will work properly  
[~,start_loc] = findpeaks(win_stim); 
[~,end_loc] = findpeaks(-win_stim);
start_stim = (start_loc-1)'/fs;

result.all = win_stim; 
result.beg = start_loc; 
result.end = end_loc; 