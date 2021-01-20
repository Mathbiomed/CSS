%% CSS_package (written by Jaehyoung Hong and Hyukpyo Hong)
%% Relieving daytime sleepiness with personalized sleep-wake patterns aligned with the circadian rhythm
clc;
clear;

% User dependant parameter
time_interval = 2;                                                          %Change according to duration of the epoch of users (min)
time_start = 24;                                                            %Change according to the starting time of sleep-wake state (h)
tau_c = 24.09;                                                              %Change according to chronotype of users
                                                                            %Morning type = 23.85, evening = 24.35 (h)
    
% read two input files
disp("---CSS_calculation---")
sleep_light = readtable('Input1_sleep_light.csv');                          %Read input1: sleep-wake pattern (1st column) and light exposure (2nd column)
waso_main = readtable('Input2_WASO_main.csv');                              %Read input2: WASO (1st column) and whether main sleep or not (2nd column)

patt = categorical(sleep_light.Var1);                                       %Convert sleep-wake to categorical array
light = table2array(sleep_light(:,2));                                      %Get light exposure

if patt(end) == 'Sleep'
    patt = [patt;'Wake'];
    light = [light;250];
end                                                                         %If sleep-wake pattern ended with 'Sleep', add 'Wake' to the last

waso = table2array(waso_main(:,1));                                         %Get WASO
main_sleep_temp = categorical(waso_main.Var2);                              %Convert main sleep to categorical array
main_sleep = find(main_sleep_temp=='M');

% Parameters
Q_max = 100; theta = 10; sigma = 3; 
mu = 4.2;coef_y = 0.8;const = 1;v_vh = 1.01;coef_x = -0.16;                 %The values of parameters in the mathematical model adopted 
gate = 1;                                                                   %in the computational package

% Baseline simulation before provided sleep-wake pattern to get stable
% sleep-wake pattern entrained to light-dark cycle
start_num = 120;                                                             %Choose 120 days which is usually enough to get stable sleep-wake patterns

if time_start < 12
    time_start = time_start + 24;
end                                                                         %When time start is less than 12 increase it as much as 24 to make natural baseline sleep

% Initial light time space 
it_first= linspace(0, 24*start_num, (24*start_num)*30+1);                   %Baseline time_interval is fixed as 2-min
it2 = it_first;
for k = 1 : length(it2)
    it2(k) =  it2(k) - 24 * (floor(it2(k)/24));
end                                                                         %Mod 24 to get exact time

i_first = zeros(1,length(it_first));
for j = 1 : length(i_first)
    if it2(j) >= 6 && it2(j) < 22-0.3
        i_first(j) = 250;
    end
end                                                                         %Light on from 6 a.m. to 21:42 a.m. as 250 lux
                                                                            %Reasonable light on timing according to            
                                                                            %baseline sleep-wake schedule (sleep occurs between 22:00-6:00) 

tspan_first = it_first;                                                     %Time 

V_00 =  [-0.5; -0.25; 5; 0; -11; 14.25];                                    %Initial for simulation (if 6th column is not too big or small, any initial is ok) 
[t0, y0] = ode15s(@(t,V) PCR(t,V,it_first,i_first,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan_first, V_00);
                                                                            %Baseline 120-day simulation
 
                                                                            
y_maximum = max(y0(end-24*30:end,5));                                       
max_index = find(y0(end-24*30:end,5)==y_maximum);                           %Save MA when MA has its highest value for forced wakefulness

y_minimum = min(y0(end-24*30:end,5));
min_index = find(y0(end-24*30:end,5)==y_minimum);                           %Save MA when MA has its lowest value for forced sleep

wake_effort_Vm = y_maximum ;
wake_effort_Vv = min(y0(end-24*30+max_index-1,4)) ;                         %Save VLPO when MA has its highest value for forced wakefulness

sleep_effort_Vm = y_minimum ;
sleep_effort_Vv = max(y0(end-24*30+min_index-1,4)) ;                        %Save VLPO when MA has its lowest value for forced sleep

% Track homeostatic sleep pressure and circadian rhythm according to       
% sleep-wake pattern and light exposure

tspan = zeros(size(patt,1),1);                                              %Base for time of users

tspan(1) = 24*start_num+time_start;                                         %Starting time of users

for j = 2 : length(tspan)
    tspan(j) = tspan(j-1) + time_interval/60;
end                                                                         %Time of users

it = tspan;                                                                 %Light exposure has same timescale with sleep-wake pattern

% Get sleep time / wake time from data 
sleep = find(patt == 'Sleep');                                              %Find sleep_finish time & sleep_num
real_time = find(diff(sleep)~=1);
real_wk_time = sleep(real_time)+1;
real_sl_time = sleep(real_time+1);                                          %Find sleep-starting, finishing time, and #(sleep)


real_wk_time = [real_wk_time;sleep(end)+1];
real_sl_time = [sleep(1);real_sl_time];                                     %Add first sleep-starting, finishing time

real_wk_time_origin = real_wk_time;
real_sl_time_origin = real_sl_time;                                         

for j = 1 : length(real_sl_time)
    if ceil(waso(j,1)/time_interval) ~= 0
    patt(real_wk_time(j,1)-ceil(waso(j,1)/time_interval):real_wk_time(j,1)-1,1)='Wake';
    real_wk_time(j,1) = real_wk_time(j,1) - ceil(waso(j,1)/time_interval);
    end   
end                                                                         %Decrease wake time as much as waso

for j = 1 : length(real_sl_time)
    real_wk_time(j) = 1 + 8 * (real_wk_time(j) - 1);
    real_sl_time(j) = 1 + 8 * (real_sl_time(j) - 1);
    real_wk_time_origin(j) = 1 + 8 * (real_wk_time_origin(j) - 1);
    real_sl_time_origin(j) = 1 + 8 * (real_sl_time_origin(j) - 1);
end                                                                         %Use more short time interval for model simulation than that of 
                                                                            %actigraphy (Use 1/8)

tspan_temp = zeros(8*(length(tspan)-1) + 1, 1);
tspan_temp(1) = tspan(1);
for j = 1 : length(tspan_temp)
    tspan_temp(j) = tspan(1) + (j-1) * (time_interval/60)/8;
end                                                                         %Use more short time interval for model simulation than that of 
                                                                            %actigraphy (Use 1/8))
tspan = tspan_temp;
st_fi = [1;real_wk_time];

% Between Initial and Day1
ts = 24 * start_num;
it2 = linspace(ts, 24*start_num+time_start, round(8*60/time_interval*(24*start_num+time_start-ts))+1);
it3 = it2;                                                                  %Time between starting time and finshing time of baseline sleep-wake schedule                                                  

 for j = 1 : length(it2)
    if it2(j) >= 24*(start_num)
        it3(j) = it3(j) - 24*start_num;   
    end
end

i2 = zeros(1,length(it3));
for j = 1 : length(i2)
    if it3(j) >= 6
        i2(j) = 250;
    end
end                                                                         %Light on after 6:00

tspan2 = it2;

tspan_total = [tspan_first' ; tspan2(2:end)' ; tspan(2:end)];               %Make total time

i_total = [i_first' ; i2(2:end)' ; light(2:end)];                           %Make total light
it_total = [it_first' ; it2(2:end)' ; it(2:end)];                           %Make time of total light

% Simulation at the gap between the end time of intial 120 days and the
% starting time of actigraphy 
V_0 = y0(end,:)';
[t1_1, y1_1] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan2, V_0);

% When users wear actiwatch very lately we assume no sleep after baseline
% sleep
D_v = -10.2 - (3.37 * 0.5) * ( const + coef_y * y1_1(:,2) + coef_x * y1_1(:,1) ) + v_vh * y1_1(:,6);
wakeup = find(D_v <= 2.46);
if isempty(wakeup) ~= 1
    if isempty(find(D_v(wakeup(1):end) > 2.46, 1)) ~= 1
        sleep_re = find(D_v(wakeup(1):end) > 2.46);
        V_0 =  [y1_1(wakeup(1)+sleep_re(1)-1,1); y1_1(wakeup(1)+sleep_re(1)-1,2); y1_1(wakeup(1)+sleep_re(1)-1,3); wake_effort_Vv; wake_effort_Vm; y1_1(wakeup(1)+sleep_re(1)-1,6)];
        [t1_1(wakeup(1)+sleep_re(1)-1:end,1), y1_1(wakeup(1)+sleep_re(1)-1:end,:)] = ode15s(@(t,V) PCR_shift(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan2(wakeup(1)+sleep_re(1)-1:end), V_0);       
    end
end


flag = 0;                                                                   %Dummy index for checking match between real and predicted sleep-
                                                                            %wake patterns
                                                                            
for day = 1 : length(st_fi)-1
    % Simulation only with light data
    V_0 = y1_1(end,:)';   
    [t1_2,y1_2] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day):st_fi(day+1),1), V_0);
    a= 0;
while(1)
    a = a + 1;
    matching_temp = Data_matching( day, tau_c, wake_effort_Vm, wake_effort_Vv, sleep_effort_Vm, sleep_effort_Vv,  y1_2, st_fi, patt, t1_2, tspan, real_sl_time, it_total,i_total );
    t1_2 = matching_temp(:,1);
    y1_2 = matching_temp(:,2:7);
    Qm = ( Q_max ./ (1+exp(-(y1_2(:,5)-theta)/sigma)) );
    % Check the difference between data and simulation
    y_temp = zeros(1+(length(y1_2(:,6))-1)/8,1);
    patt_simul = categorical(y_temp(1:end,1),1);
    for j = 1 : length(patt_simul)
        if Qm(1+8*(j-1),1) > 1
            patt_simul(j,1) = 'Wake';                                       %Qm-firing rate of MA population-is closely related with arousal (Qm > 1 : awake)
        elseif Qm(1+8*(j-1),1) <= 1
            patt_simul(j,1) = 'Sleep';                                      %Qm <= 1 : asleep
        end
    end
    diff_patt = find(patt(1+(st_fi(day)-1)/8:1+(st_fi(day+1)-1)/8,1) ~= patt_simul);                  % Compare Wake/Sleep state simulated by model with Wake/Sleep state from actigraphy
    if isempty(diff_patt) == 1 
        break;
    end
    if diff_patt(1) == length(patt_simul)
        break;
    end
    if a == 10
        disp('Too much time, Error occurs');
        flag = 1;
        break
    end
end
    if flag == 1
        break;
    end
    y1_1 = [y1_1(1:end-1,:);y1_2(1:end,:)];
    t1_1 = [t1_1(1:end-1,:);t1_2(1:end,:)];
end

t_total = t1_1(length(tspan2):end,1);
y_total = y1_1(length(tspan2):end,:);

% Use first sleep time in 12:00-12:00 rather than 0:00-24:00
sleep_time = mod(tspan(real_sl_time(main_sleep(1))),24);
if sleep_time<12
    sleep_time = sleep_time + 24;
end

sleep_time = sleep_time - 12;

day_temp = floor((tspan(real_sl_time(main_sleep(end)))-tspan(real_sl_time(main_sleep(1)))-(24-sleep_time))/24);
day_remain = mod(tspan(real_sl_time(main_sleep(end)))-tspan(real_sl_time(main_sleep(1)))-(24-sleep_time),24);

day_new = 1 + day_temp;                                                     %First day + day_num

if day_remain > 0
   day_new = day_new + 1; 
end

base_of_sleep = floor((tspan(real_sl_time(main_sleep(1)))-(24*start_num+12))/24); %Find the day which does not have main sleep before first main sleep

no_sleep = [];                                                              %We need to find the day which does not have main sleep
if length(main_sleep) < day_new
    day_set = zeros(length(main_sleep),1);
   for j =  1: length(main_sleep)
       day_set(j) = floor((tspan(real_sl_time(main_sleep(j)))-tspan(real_sl_time(main_sleep(1)))-(24-sleep_time))/24)+2;
   end
   
   for  j = 1 : day_new
       trick_day = find(day_set == j);
       if length(trick_day) == 2
           if (isempty(find(day_set==j-1,1)) == 1 && isempty(find(day_set==j+1,1)) == 0)
               day_set(trick_day(1)) = day_set(trick_day(1)) - 1;
           elseif (isempty(find(day_set==j-1,1)) == 0 && isempty(find(day_set==j+1,1)) == 1)
               day_set(trick_day(2)) = day_set(trick_day(2)) + 1;    
           else
               trick1 = mod(tspan(real_sl_time(main_sleep(trick_day(1))))-tspan(real_sl_time(main_sleep(1)))-(24-sleep_time),24);
               trick2 = mod(tspan(real_sl_time(main_sleep(trick_day(2))))-tspan(real_sl_time(main_sleep(1)))-(24-sleep_time),24);
               if trick1 <= 6
                   day_set(trick_day(1)) = day_set(trick_day(1)) - 1;
               end
               if trick2 >= 6+12
                   day_set(trick_day(2)) = day_set(trick_day(2)) + 1;
               end
           end
       end
   end
   
   for j = 1 : day_new
      if (isempty(find(day_set==j,1)) == 1)
        no_sleep = [no_sleep; j];
      end
   end
end

[CSS, CSS_temp,normal_sleep_amount] = Calcul_CSS(length(no_sleep)+base_of_sleep, coef_x, coef_y, v_vh, tau_c, gate, main_sleep, y1_1, t1_1, tspan2, real_sl_time, real_wk_time, real_sl_time_origin, real_wk_time_origin);
                                                                            %CSS calculation

%% H_C.png

mkdir Output

Font_size = 15;
Font_weight = 'bold'; 
axis_thick = 1.5;                                                           %Font, axis                
                                                                            
t_total = t1_1(length(tspan2):end,1);
y_total = y1_1(length(tspan2):end,:);

Q_m = ( Q_max ./ (1+exp(-(y_total(:,5)-theta)/sigma)) );
sleep_model = find(Q_m <= 1);
sleep_change = find(diff(sleep_model)~=1);                                  %Find sleep end
sleep_start = [1;sleep_change+1];                                           %Find sleep start
sleep_end = [sleep_change;length(sleep_model)];                             %Find sleep end

H = y_total(:,6);
C = (3.37 * 0.5) * ( const + coef_y * y_total(:,2) + coef_x * y_total(:,1) );
 
D_v = -10.2 - (3.37 * 0.5) * ( const + coef_y * y_total(:,2) + coef_x * y_total(:,1) ) + v_vh * y_total(:,6);
D_up = (2.46 + 10.2 + C)/v_vh;

my_xlim = [0,t_total(end)-t_total(1)];
my_ylim = [floor(min(H))-1, floor(max(H))+2];
temp_tick = find(mod(t_total,6)==0);
if temp_tick(1) ~= 1
    my_xtick = linspace(t_total(temp_tick(1))-t_total(1),t_total(temp_tick(end))-t_total(1),length(temp_tick));
else
    my_xtick = linspace(t_total(temp_tick(2))-t_total(temp_tick(1)),t_total(temp_tick(end))-t_total(temp_tick(1)),length(temp_tick)-1);
end

my_ytick = linspace(floor(min(H)),floor(max(H))+1,3);
my_xlabel = zeros(1,length(my_xtick));
for j = 1 : length(my_xlabel) 
    my_xlabel(j) = mod(my_xtick(j)+t_total(1),24);
end

figure(1)
lw = 5;
h1 = plot(t_total-t_total(1),D_up,'Color',[255,194,71]/255, 'LineWidth',lw);
hold on
h2 = plot(t_total-t_total(1),H,'Color','k', 'LineWidth',lw);
hold on
for j = 1 : length(real_sl_time)
    v = [t_total(sleep_model(sleep_start(j)))-t_total(1) floor(min(H))-1; t_total(sleep_model(sleep_end(j)))-t_total(1) floor(min(H))-1 ; t_total(sleep_model(sleep_end(j)))-t_total(1)  floor(max(H))+2 ; t_total(sleep_model(sleep_start(j)))-t_total(1) floor(max(H))+2];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor',[112/255 219/255 217/255],'Edgealpha',0)
    alpha(0.4)
    hold on
end

set(gca,'XLim',my_xlim)
set(gca,'YLim',my_ylim)
set(gca,'xtick',my_xtick)
set(gca,'ytick',my_ytick)
xticklabels(my_xlabel)
xlabel('Time (h)')
ylabel('Homeostatic sleep pressure')
set(gca, 'FontSize', Font_size)
set(gca, 'FontWeight', Font_weight)
set(gca,'linewidth',axis_thick)

saveas(gcf,[pwd '/Output/H_C.png'])

%% CSS.csv

actual_sleep = zeros(length(main_sleep),1);
nece_sleep = normal_sleep_amount;
temp = CSS_temp;
suff = string(temp);
CSS_change = zeros(length(main_sleep),1);
 
for j = 1 : length(suff)
    if temp(j) == 1
        suff(j) = 'S';
    else
        suff(j) = 'I';
    end
end
 
for j = 1 : length(main_sleep)
     actual_sleep(j) = (real_wk_time(main_sleep(j))-real_sl_time(main_sleep(j)))*(time_interval/60)/8;
end

for j = 1 : length(CSS_change)
    no_sleep_day = length(find(no_sleep <= j));                             %Find number of no_sleep_day that must be insufficiet sleep 
    CSS_change(j) = 100*sum(CSS_temp(1:j))/(j+no_sleep_day+base_of_sleep);                  
end
 
writematrix([nece_sleep, actual_sleep,suff,CSS_change],[pwd '/Output/CSS.csv'])
 
normal_time = normal_sleep_amount;

%% CSS.png

real_color = [0,0,0]/255;
model_color = [156,156,156]/255;
edge_color = [0,0,0]/255;
line_width = 1.5;

% space for all sleep (including nap)
sleep_interval = zeros(length(real_sl_time),1);
duration_real = zeros(length(real_sl_time),1);

% space for min sufficient sleep
duration_model = zeros(length(main_sleep),1);

% Get sleep time and duration for all sleep
for j = 1 : length(real_sl_time)
    sleep_interval(j) = tspan(real_sl_time(j))-tspan(real_sl_time(1));
    duration_real(j) = tspan(real_wk_time(j))-tspan(real_sl_time(j));
end

% Min sufficient sleep for main sleep
for j = 1 : length(main_sleep)
    duration_model(j) = normal_time(j);
end

% Use first sleep time in 17:00-17:00 rather than 0:00-24:00
sleep_time = mod(tspan(real_sl_time(1)),24);
if sleep_time<17
    sleep_time = sleep_time + 24;
end

sleep_time = sleep_time - 17;

% Calculate number of day to display

if (main_sleep(end) == length(real_sl_time) && duration_model(end)>duration_real(end))  
    day_temp = floor((duration_model(end)-duration_real(end)+tspan(real_wk_time(end))-tspan(real_sl_time(1))-(24-sleep_time(1)))/24);
    day_remain = mod(duration_model(end)-duration_real(end)+tspan(real_wk_time(end))-tspan(real_sl_time(1))-(24-sleep_time(1)),24);
else
    day_temp = floor((tspan(real_wk_time(end))-tspan(real_sl_time(1))-(24-sleep_time))/24);
    day_remain = mod(tspan(real_wk_time(end))-tspan(real_sl_time(1))-(24-sleep_time),24);
end

day_new = 1 + day_temp;                                                     %First day + day_num

if day_remain > 0
   day_new = day_new + 1; 
end

day_new = day_new + base_of_sleep;                                          %Add the day before first main sleep

% Position & scale
day_diff_bot = 1.5;
bar_thick = 2;

day_thick = day_diff_bot + bar_thick*2;

main_index = 1;

my_ytick = zeros(day_new,1);
for j = 1 : day_new
   my_ytick(j) = day_thick*(j-1) + day_diff_bot + bar_thick;
end

ytick_label = string({''});
for j = 1 : day_new
    tickname = sprintf('Day%i', day_new-j+1);
    ytick_label(j,1) = tickname;
end

% Draw Actual sleep and necessary sleep
figure(3)
for j = 1 : length(real_sl_time)   
    find_day = day_new-2-floor((sleep_interval(j)-(24-sleep_time))/24)-base_of_sleep;     %Find y-axis
    find_time = mod(sleep_time + sleep_interval(j),24);                     %Find x-axis
    if find_time + duration_real(j) <= 24
        rectangle('Position',[find_time day_thick*find_day+bar_thick+day_diff_bot duration_real(j) bar_thick],'Facecolor',real_color,'Edgecolor',edge_color,'Linewidth',line_width)
        hold on                                                             %Draw actual sleep
    else
        rectangle('Position',[find_time day_thick*find_day+bar_thick+day_diff_bot 24-find_time bar_thick],'Facecolor',real_color,'Edgecolor',edge_color,'Linewidth',line_width)
        hold on                                                             %If wake after 17:00 draw remaining sleep on next day
        rectangle('Position',[0 day_thick*(find_day-1)+bar_thick+day_diff_bot duration_real(j)-(24-find_time) bar_thick],'Facecolor',real_color,'Edgecolor',edge_color,'Linewidth',line_width)
        hold on
    end
    if main_index <= length(main_sleep)                                     %If we have more calculated necessary sleep, draw it
        if main_sleep(main_index) == j                                      %Draw necessary sleep if main sleep episode
            if find_time + duration_model(main_index) <= 24
                rectangle('Position',[find_time day_thick*find_day+day_diff_bot duration_model(main_index) bar_thick],'Facecolor',model_color,'Edgecolor',edge_color,'Linewidth',line_width)
                hold on                                                     %Draw actual sleep
                main_index = main_index + 1;
            else
                rectangle('Position',[find_time day_thick*find_day+day_diff_bot 24-find_time bar_thick],'Facecolor',model_color,'Edgecolor',edge_color,'Linewidth',line_width)
                hold on                                                     %If wake after 17:00 draw remaining sleep on next day
                rectangle('Position',[0 day_thick*(find_day-1)+day_diff_bot duration_model(main_index)-(24-find_time) bar_thick],'Facecolor',model_color,'Edgecolor',edge_color,'Linewidth',line_width)
                hold on
                main_index = main_index + 1;
            end
        end
    end
end

set(gca, 'FontSize', Font_size)
set(gca, 'FontWeight', Font_weight)
set(gca,'linewidth',axis_thick)
axis([0 24 0 day_thick*day_new+day_diff_bot])
yticks(my_ytick)
xticks([0 6 12 18 24])
xticklabels({'17:00','','5:00','','17:00'})
yticklabels([ytick_label])
xlabel('Time (h)')

txt = sprintf('CSS = %.1f (%%)', floor(CSS)+round(10*(CSS-floor(CSS)))/10);
text(16,day_thick*day_new+day_diff_bot/2,txt,'FontSize',Font_size,'FontWeight', Font_weight)

saveas(gcf,[pwd '/Output/CSS.png'])
