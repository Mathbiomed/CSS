%% Sleep_make (written by Jaehyoung Hong)
%% Relieving daytime sleepiness with personalized sleep-wake patterns aligned with the circadian rhythm

clc;
clear;

% User dependant parameter
time_interval = 2;                                                          %Change according to duration of the epoch of users (min)

% Make schedule & pattern
data = readtable('Input1_sleep.csv');
sleep_start_time =  mod(table2array(data(:,2)),24);
sleep_end_time =  mod(table2array(data(:,4)),24);

sleep_start_day = table2array(data(:,1));
sleep_start_day = hours(sleep_start_day-sleep_start_day(1))/24;
sleep_end_day = table2array(data(:,3));
sleep_end_day = hours(sleep_end_day-sleep_end_day(1))/24;
if sleep_end_time(1) < sleep_start_time(1)
    sleep_end_day = sleep_end_day+1;
end

time_start = floor(sleep_start_time(1))-1;
if time_start == 0 
    time_start = 24;
end

start_day = sleep_start_day(1);
real_sl_time = [];
real_wk_time = [];

for j = 1:length(sleep_start_day)
    real_sl_time = [real_sl_time;sleep_start_day(j)*24*60/time_interval+round(sleep_start_time(j)*60/time_interval)-time_start*60/time_interval+1]; 
    real_wk_time = [real_wk_time;sleep_end_day(j)*24*60/time_interval+round(sleep_end_time(j)*60/time_interval)-time_start*60/time_interval+1];
end

tspan = zeros(real_wk_time(end)+1,1);
te =  time_interval;
tspan(1) = time_start;
for j = 2 : length(tspan)
    tspan(j) = tspan(j-1) + time_interval/60;
end

patt = string({'Wake'});
for j = 2:length(tspan)
    patt(j,1) = 'Wake';
end
for j = 1:length(real_wk_time)
    for k = real_sl_time(j):real_wk_time(j)-1
        patt(k) = 'Sleep';
    end
end
patt = categorical(patt);

writematrix(patt,'Maked_sleep.csv')                                          %Make csv file

time_start