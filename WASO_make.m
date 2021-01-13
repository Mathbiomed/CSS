%% WASO_make (written by Jaehyoung Hong)
%% Relieving daytime sleepiness with personalized sleep-wake patterns aligned with the circadian rhythm
clc;
clear;

time_interval = 2;                                                          %Change according to duration of the epoch of users (min)

% Read input file
WASO_sleep = readtable('Input2_WASO.csv');                                  %Read input1: Sleep(0)/Wake(1) (1st column) and sleep-wake pattern (2nd column)

WASO = table2array(WASO_sleep(:,1));                                        %Get WASO 
patt = categorical(WASO_sleep.Var2);                                        %Convert sleep-wake to categorical array

sleep = find(patt=='Sleep');                                                %To calculate WASO, we only need sleep status
sleep_change = find(diff(sleep)~=1);                                        %Find sleep end
sleep_start = [1;sleep_change+1];                                           %Find sleep start
sleep_end = [sleep_change;length(sleep)];                                   %Find sleep end

maked_WASO = zeros(length(sleep_start),1);                                  %Base of WASO for all sleep episodes

for j = 1 : length(sleep_start)
    maked_WASO(j) = time_interval * length(find(WASO(sleep(sleep_start(j)):sleep(sleep_end(j)),:)==1)); 
end                                                                         %1 indicate wake during sleep (WASO) * duration of the epoch

writematrix(maked_WASO,'Maked_WASO.csv')                                    %Make csv file






