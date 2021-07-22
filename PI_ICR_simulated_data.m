%% PI-ICR_simulated_data.m
% 
% Simulating PI-ICR data using expected cyclotron frequencies/beam makeups
%
% Language: MATLAB
%
% Sam Porter
% University of British Columbia, TRIUMF
%
% last updated on 4/22/2021
%
%% User Inputs %%

% From Self-Generated Tacc List %%

clearvars

tacc_list = [.6000002,.6000141,.600030,.6000479,.6000658,.6001155,.6001513,.6001702,.600203]; % List of Accumulation Times

new_filename = '/Users/wsporter/Documents/Physics_Research/TITAN/PIICR_Analysis/Simulated_Data/56V_small_forplot/56V_sim';

spot_freq_list = [1005878.022,1005860.729,1005856.879,1005842.204,1005750.655]; % Expected w_c of spots in the trap
R = 4; % Average radius of spots
pos_deviation = 0.5; % Sigma of Gaussian distribution for X/Y positions
TOF = 0.000006; % Mean of Gaussian distribution for time-of-flight 
tof_deviation = 0.000003 % Sigma of Guassian distribution for time-of-flight

counts_per_spot = 100;
w_minus = 6150; 
amp = 0.01; % Amplitude of Sine Dependency
phase_const = 10; % Phase Constant of Sine Dependency 
ref_phase = 150; % Phase of Reference Spot

% NOTE: TRAP CENTER IS TREATED AS (0,0)
%% Code Body %% 

% Determine Number of Actual vs. Ref Files %%

ctr = 0;
for object = tacc_list
    if object > 0
        ctr = ctr + 1;
    end
end

tacc_list = tacc_list(1:ctr);

data_filename_list = {};
ref_phase = ref_phase*pi/180;

ii = 1;

% From Tacc & w_c Determine X/Y %%

for i = tacc_list
    x_data = [];
    y_data = [];
    TOF_data = [];
    trigger_data = [];
    
    for j = spot_freq_list
        w_c_sin_shifted = j + amp*sin((2*pi*w_minus)*i + phase_const); % Sine dependence of spots on w_minus/Tacc
        phase_i = 2*pi*i*w_c_sin_shifted + ref_phase;
        
        if phase_i > pi
            phase_i = phase_i - 2*pi;
        end
        
        x = R*cos(phase_i); % NOTE: we're shifting the phase but not the radius, which means simulated radii will not follow a sinusoid.
        y = R*sin(phase_i);
        
        for k = 1:counts_per_spot % For the number of points per spot, generate data within the specific distribution
            x_data = cat(1,x_data,random('Normal',x,pos_deviation));
            y_data = cat(1,y_data,random('Normal',y,pos_deviation));
            TOF_data = cat(1,TOF_data,random('Normal',TOF,tof_deviation));
            trigger_data = cat(1,trigger_data,k+random('Uniform',0,1)); 
        end
        
    end
    
    % Write New Data File %%
    
    file_data = cat(2,x_data,y_data,TOF_data,trigger_data);
    data_filename = strcat(new_filename,num2str(i),'.csv');
    dlmwrite(data_filename,file_data)
    
    data_filename_list{ii,1} = data_filename;
    
    ii = ii + 1;
end

% Write New Reference Data File %%

if ref_phase > pi
    ref_phase = ref_phase - 2*pi;
end

x = R*cos(ref_phase);
y = R*sin(ref_phase);

x_data = [];
y_data = [];
TOF_data = [];
trigger_data = [];

for k = 1:counts_per_spot
    x_data = cat(1,x_data,random('Normal',x,pos_deviation));
    y_data = cat(1,y_data,random('Normal',y,pos_deviation));
    TOF_data = cat(1,TOF_data,random('Normal',TOF,tof_deviation));
    trigger_data = cat(1,trigger_data,k+random('Uniform',0,1)); 
end

file_data = cat(2,x_data,y_data,TOF_data,trigger_data);
data_filename = strcat(new_filename,'ref.csv');
dlmwrite(data_filename,file_data)

% Construct List of Files .csv %%

final_filename_list = {};
ii = 1;

while ii <= ctr
    final_filename_list{ii,1} = data_filename_list{ii,:};
    
    ii = ii + 1;
end

final_filename_list{ctr+1,1} = data_filename;

tacc_list = transpose(tacc_list)
old_tacclist = cat(1,tacc_list,0);
ref_assign_list(1:ctr+1,1) = 1;
ref_assign_list(ctr+1,1) = NaN;

T = table(final_filename_list,old_tacclist,ref_assign_list);

data_filename = strcat(new_filename,'.csv');
writetable(T,data_filename,'WriteVariableNames',false)


