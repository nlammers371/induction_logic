% script to write csv data to data structures
clear
close all

% set basic paths
addpath('utilities')
DataPath = 'C:\Users\nlamm\projects\induction_project\dat\20200807\';

% get list of data sets
FileList = dir([DataPath '*.xlsx']);
expStrings = {FileList.name};
expStrings = expStrings(~contains(expStrings,'~$'));
expIDs = 1:length(expStrings);
% specify time res
dT = 30;
% define key model architecture parameters 
alpha_frac = 1.302/(1.302 +  1.085/2); % Memory (approximated for now)

% initialize structure
trace_structure = struct;
i_pass = 1;
for e = expIDs
  raw_table = readtable([DataPath expStrings{e}]);
  raw_array = raw_table{:,:};
  time_vec = 0:dT:dT*size(raw_array,1)-1;
  divFactor = 10^ceil(log10(size(raw_array,2)));
  for i = 1:size(raw_array,2)
    trace_structure(i_pass).fluo = raw_array(:,i)';
    trace_structure(i_pass).time = time_vec;
    trace_structure(i_pass).setID = e;
    trace_structure(i_pass).setName = expStrings{e};
    trace_structure(i_pass).ParticleID = e + i/divFactor;
    trace_structure(i_pass).alpha_frac = alpha_frac;
    trace_structure(i_pass).Tres = dT;
    i_pass = i_pass + 1;
  end
end

% save
save([DataPath 'trace_structure.mat'],'trace_structure')