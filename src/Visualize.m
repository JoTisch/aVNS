%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                OUTPUT                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data


% --- Data Overview ---
output_plot(raw,rpeaks, win_stim,maf,seg,edges,fs,x1,prc,colMeans,beta,label,stimstr); 

%- Raw Signal Stimulation -
fit_IE(stim_alpha, stimRR, class,stimstr(stim_flag),60)
RdR_fit(stim_alpha, stimRR, class,1,60,0) 
RdR_fit(stim_alpha, stimRR, stimstr(stim_flag),1,60,1)

%--- Mean Filter Stimulation ----
RdR_fit(stim_alpha, mf_stimRR, class,1,60,0) 
RdR_fit(stim_alpha, mf_stimRR, stimstr(stim_flag),1,60,1) 

%--- Adaptive Filter Stimulation ---
RdR_fit(stim_alpha, af_stimRR, stimstr(stim_flag),1,60,1)
