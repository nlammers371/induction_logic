clear
close all

% set basic paths
addpath('utilities')
DataPath = 'C:\Users\nlamm\projects\induction_project\dat\20200807\';

% load WT data
wt_table = readtable([DataPath 'WT.xlsx']);
wt_array = wt_table{:,:};
dT = 0.5;
time_vec = 0:dT:dT*size(wt_array,1)-dT;
% set parameters for autocorr analysis
n_lags = 20;
n_boots = 100;

% % %% calculate autocorrelation
% % [wt_autocorr, a_boot_errors, wt_dd, dd_boot_errors] = ...
% %     weighted_autocorrelation(wt_array, n_lags, 1,n_boots,ones(size(wt_array)));

% employ mild smoothing filter to reduce noisy fluctuations
% wt_traces_filtered = wt_array;
% for i = 1:size(wt_array,2)
%   wt_traces_filtered(:,i) = imgaussfilt(wt_traces_filtered(:,i),.1);
% end

%% calculate autocorrelation
min_dp = 50;
autocorr_mat = NaN(n_lags+1,1);
weight_mat = NaN(n_lags+1,1);
fragment_array = NaN(size(wt_array,1)+1,1);
% contig_vec = [];
iter = 1;
for i = 1:size(wt_array,2)
  trace_raw = wt_array(:,i)';
  trace_trunc = [0 trace_raw(find(trace_raw~=0,1):find(trace_raw~=0,1,'last')) 0];
  
  % subdivide trace into contiguous active perods
  off_ids = find(trace_trunc==0);
  contig_fragments = diff(off_ids);  
  long_fragments = find(contig_fragments > n_lags);
%   contig_vec = [contig_vec contig_fragments];
  
  for j = 1:length(long_fragments)    
    fragment = trace_trunc(off_ids(long_fragments(j)):off_ids(long_fragments(j)+1));
    
    autocorr_mat(:,iter) = autocorr(fragment,n_lags);
    weight_mat(:,iter) = length(fragment):-1:length(fragment)-n_lags;
    
    fragment_array(1:length(fragment),iter) = fragment;
    
    iter = iter + 1;
  end

end

% calculate weighted average
autocorr_boot_array = NaN(n_lags+1,n_boots);
autocorr_dd_boot_array = NaN(n_lags-1,n_boots);
sample_index = 1:size(fragment_array,2);
for n = 1:n_boots
  boot_indices = randsample(sample_index,length(sample_index),true);
  autocorr_boot_array(:,n) = nansum(autocorr_mat(:,boot_indices).*weight_mat(:,boot_indices),2) ./ nansum(weight_mat(:,boot_indices),2);
  autocorr_dd_boot_array(:,n) = diff(autocorr_boot_array(:,n),2);
end

autocorr_mean = nanmean(autocorr_boot_array,2);
autocorr_ste = nanstd(autocorr_boot_array,[],2);
autocorr_dd_mean = nanmean(autocorr_dd_boot_array,2);
autocorr_dd_ste = nanstd(autocorr_dd_boot_array,[],2);
%% Looks like there's a ton of heterogeneity
slow_ids = find(autocorr_mat(6,:)>0.5);
fast_ids = find(autocorr_mat(6,:)<=0.5);

figure(1);
yyaxis left
plot(nanmean(autocorr_mat(:,fast_ids),2))

yyaxis right
plot(diff(nanmean(autocorr_mat(:,fast_ids),2),2))

