function [SS,SS_temp,necc_sleep_amount] = Calcul_SS(no_sleep,coef_x, coef_y, v_vh, tau_c, gate, main_sleep, y1_1, t1_1, tspan2, real_sl_time, real_wk_time, real_sl_time_origin, real_wk_time_origin)

Q_max = 100;
theta = 10;
sigma = 3;
mu = 4.2;
const =1;

% Calculate Dv to determine sufficient or not
D_v = -10.2 - (3.37 * 0.5) * ( const + coef_y * y1_1(length(tspan2):end,2) + coef_x * y1_1(length(tspan2):end,1) ) + v_vh * y1_1(length(tspan2):end,6);

SS_Dv = zeros(length(real_wk_time),1);                                      %Base for Dv at each sleep offset of 'all' sleep
SS_temp = zeros(length(main_sleep),1);                                      %Base for sufficient or not for each main sleep (1='sufficient', 0='insufficient)  
SS_temp2 = zeros(length(main_sleep),1);                                     %Base for Dv at each sleep offset of 'main' sleep

for j = 1 : length(real_wk_time)
    SS_Dv(j) = D_v(real_wk_time(j));
end                                                                         %Calculate Dv at each sleep offset of 'all' sleep

for j = 1 : length(main_sleep)
    SS_temp2(j) = SS_Dv(main_sleep(j));
end                                                                         %Aggregate Dv at each sleep offset of 'main' sleep

early_confirm = zeros(length(main_sleep),1);
for j = 1 : length(main_sleep)
    confirm_early = find(D_v(real_sl_time(main_sleep(j)):real_wk_time(main_sleep(j)),:) > 2.46);
    if isempty(confirm_early) == 1
        early_confirm(j) = 1;
    end
end                                                                         %If Dv <= 2.46 during sleep episode, we define it as early sleep

y_temp = y1_1(length(tspan2):end,:);
time_temp = t1_1(length(tspan2):end,:);

necc_sleep_amount = zeros(length(main_sleep),1);                          %Calculate necessary sleep for each main_sleep_episode
for j = 1 : length(main_sleep)
    if early_confirm(j) == 1
        V_0 = y_temp(real_sl_time(main_sleep(j))-8,:);                      %Start from wake before sleep
        time_early = zeros(8*30*24,1);                                      %Assume reach sleep threshold in 24h-wake state
        for t_index = 1 : length(time_early)
            time_early(t_index) = time_temp(real_sl_time(main_sleep(j))-8,:) + 2/(60*8) * (t_index-1);
        end
        
        i_early = 250 * ones(8*30*24,1);                                    %Assume light is 250lux which is baseline lux for wake
        
        [t_early,y_early] = ode15s(@(t,V) PCR(t,V,time_early,i_early,tau_c,mu, v_vh, coef_x, coef_y, const, gate), time_early, V_0);
        Q_temp = ( Q_max ./ (1+exp(-(y_early(:,5)-theta)/sigma)) );
        sleep_start = find(Q_temp <= 1);                                    %Find the time when it normally triggered sleep state
        
        V_0 = y_early(sleep_start(1),:);                                    %Start from sleep onset
        
        time_temp2 = t_early(sleep_start(1));
        
        time_early2 = zeros(8*30*24,1);                                     %Assume reach sleep threshold in 24h-sleep state
        for t_index = 1 : length(time_early2)
            time_early2(t_index) = time_temp2 + 2/(60*8) * (t_index-1);
        end
        
        i_early2 = zeros(8*30*24,1);                                        %Assume light is 0lux which is baseline lux for sleep
        
        [t_early2,y_early2] = ode15s(@(t,V) PCR(t,V,time_early2,i_early2,tau_c,mu, v_vh, coef_x, coef_y, const, gate), time_early2, V_0);
        
        Dv_temp = -10.2 - (3.37 * 0.5) * ( const + coef_y * y_early2(:,2) + coef_x * y_early2(:,1) ) + v_vh * y_early2(:,6);
        
        wake_start = find(Dv_temp <= 2.46);                                 %Find the time when it normally reaching wake threshold
        
        necc_sleep_amount(j) =  (wake_start(1)-1)*2/(60*8);                 %Calculate necessary sleep
    else                                                                    %If not early sleep
        if SS_temp2(j) <= 2.46                                              %When sufficient sleep
            Dv_temp = D_v(real_sl_time(main_sleep(j)):real_wk_time(main_sleep(j)));
            over_sleep_thres = find(Dv_temp<=2.46);
            early_normal = find(diff(over_sleep_thres)~=1);                 %When sleep start in bistable region and enter sleep region during sleep
                                                                            %This kind of difference can be happened
            if isempty(early_normal) == 1
                necc_sleep_amount(j) = (over_sleep_thres(1)-1)*2/(60*8);
            else
                necc_sleep_amount(j) = (over_sleep_thres(early_normal(1)+1)-1)*2/(60*8);
            end
        else                                                                %If insufficient sleep
            V_0 = y_temp(real_wk_time(main_sleep(j))-8,:);                  %Start from sleep
            time_more = zeros(8*30*24,1);                                   %Assume reach sleep threshold in 24h-sleep state
            for t_index = 1 : length(time_more)
                time_more(t_index) = time_temp(real_wk_time(main_sleep(j))-8,:) + 2/(60*8) * (t_index-1);
            end
            i_more = zeros(8*30*24,1);                                      %Assume light is 0lux which is baseline lux for sleep       
            
            [t_more,y_more] = ode15s(@(t,V) PCR(t,V,time_more,i_more,tau_c,mu, v_vh, coef_x, coef_y, const, gate), time_more, V_0);
            Dv_temp = -10.2 - (3.37 * 0.5) * ( const + coef_y * y_more(:,2) + coef_x * y_more(:,1) ) + v_vh * y_more(:,6);
            
            wake_start = find(Dv_temp <= 2.46);                             %Find the time when it normally reaching wake threshold
            necc_sleep_amount(j) =  (wake_start(1)-9)*2/(60*8) + (real_wk_time(main_sleep(j))-real_sl_time(main_sleep(j)))*2/(60*8);
                                                                            %Calculate necessary sleep = Actual sleep duration + more needed time
        end
    end
end

% For early sleep, sleep duration must be larger than normal sleep amount
% For the others, Dv<=2.46 is enough for being sufficient sleep
for j = 1 : length(main_sleep)
    if early_confirm(j) == 1
        if (real_wk_time(main_sleep(j))-real_sl_time(main_sleep(j)))*2/(60*8) > necc_sleep_amount(j)
            SS_temp(j) = 1;
        end
    else
        if SS_Dv(main_sleep(j)) <= 2.46
            SS_temp(j) = 1;
        end
    end
end

for j = 1 : length(main_sleep)
    if j ~= length(main_sleep)
        if (SS_temp(j) == 0) && (main_sleep(j) + 1 < main_sleep(j+1))
            near_range = real_sl_time_origin(main_sleep(j):main_sleep(j+1)-1)-real_wk_time_origin(main_sleep(j));
            near_main = find(near_range*2/(60*8)<=3);
            sl_du = (real_wk_time(main_sleep(j)+near_main(1)-1:main_sleep(j)+near_main(end)-1) - real_sl_time(main_sleep(j)+near_main(1)-1:main_sleep(j)+near_main(end)-1))*2/(60*8);
            if necc_sleep_amount(j) <= sum(sum(sl_du))
                SS_temp(j) = 1;
            end
        end
    else
        if (SS_temp(j) == 0) && (main_sleep(j) < length(real_sl_time))
            near_range = real_sl_time_origin(main_sleep(j):length(real_sl_time))-real_wk_time_origin(main_sleep(j));
            near_main = find(near_range*2/(60*8)<=3);
            sl_du = (real_wk_time(main_sleep(j)+near_main(1)-1:main_sleep(j)+near_main(end)-1) - real_sl_time(main_sleep(j)+near_main(1)-1:main_sleep(j)+near_main(end)-1))*2/(60*8);
            if necc_sleep_amount(j) <= sum(sum(sl_du))
                SS_temp(j) = 1;
            end
        end        
    end
end


SS = 100 * sum(sum(SS_temp)) / (length(main_sleep)+no_sleep);   
end