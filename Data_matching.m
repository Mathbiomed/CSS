function matching_data = Data_matching( day, tau_c, wake_effort_Vm, wake_effort_Vv, sleep_effort_Vm, sleep_effort_Vv,  y1_2, st_fi, patt, t1_2, tspan, real_sl_time, it_total,i_total )

Q_max = 100;
theta = 10;
sigma = 3;
mu = 4.2;
const = 1;
coef_x = -0.16;
coef_y = 0.8;
v_vh = 1.01;
gate = 1;

% Calculate MA firing rate which decide whether the model state is wake or
% not
Qm = ( Q_max ./ (1+exp(-(y1_2(:,5)-theta)/sigma)) );

% Check the difference between data and simulation
y_temp = zeros(1+(length(y1_2(:,6))-1)/8,1);
patt_simul = categorical(y_temp(1:end,1),1);

for j = 1 : length(patt_simul)
    % Qm-firing rate of MA population-is closely related with arousal (Qm > 1 : awake)
    if Qm(1+8*(j-1),1) > 1                             
        patt_simul(j,1) = 'Wake';     
    % Qm <= 1 : asleep
    elseif Qm(1+8*(j-1),1) <= 1                                                                   
        patt_simul(j,1) = 'Sleep';                                         
    end
end

% Compare Wake/Sleep state simulated by model with Wake/Sleep state from actigraphy
diff_patt = find(patt(1+(st_fi(day)-1)/8:1+(st_fi(day+1)-1)/8,1) ~= patt_simul);     

if isempty(diff_patt) == 1
    
% Last time point will be fixed at next wake-sleep block
elseif diff_patt(1) ~= length(patt_simul)
    % When first dismatch occur, model is at the wake state
    if patt_simul(diff_patt(1)) == 'Wake'
        % If one-time point before dismatching has different sleep/wake
        % state with that of dismatching find the exact time trasition of
        % wake state occur
        if patt_simul(diff_patt(1)-1) == 'Sleep'
            Qm_temp = ( Q_max ./ (1+exp(-(y1_2(1+(diff_patt(1)-1)*8-7:1+(diff_patt(1)-1)*8,5)-theta)/sigma)) );
            bifur = find(Qm_temp > 1);
            diff_modi = 1 + 8 * (diff_patt(1)-1) - 8 + bifur(1);
        else
            diff_modi = 1 + 8 * (diff_patt(1)-1);
        end
        % Calculate Dv which decide whether model is in the bistable region
        % or not
        Dv = -10.2 - (3.37 * 0.5) * ( const + coef_y * y1_2(diff_modi,2) + coef_x * y1_2(diff_modi,1) ) + v_vh * y1_2(diff_modi,6);
        % When model is in the bistable region
        if Dv >= 1.45
            % Change Vv,Vm as the value whenn the bifurcation occur
            V_0 = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);0.00;-4.78;y1_2(diff_modi,6)];
            [t1_2(diff_modi:end,1),y1_2(diff_modi:end,:)] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day)+diff_modi-1:st_fi(day+1),1), V_0);
        else
            if patt_simul(diff_patt(1)-1) == 'Wake'
                % Take the value of sleep effort as min of the Vm during
                % 120th day of First 120 days-initialization
                if t1_2(diff_modi,1) < t1_2(end-9,1)
                    V_0 = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);sleep_effort_Vv;sleep_effort_Vm;y1_2(diff_modi,6)];
                    % Give sleep effort untill the time 'REST-S' of actigraphy is finished
                    [t1_2(diff_modi:end-8,1),y1_2(diff_modi:end-8,:)] = ode15s(@(t,V) PCR_shift(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day)+diff_modi-1:st_fi(day+1)-8,1), V_0);
                elseif t1_2(diff_modi,1) == t1_2(end-9,1)
                    V_0 = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);sleep_effort_Vv;sleep_effort_Vm;y1_2(diff_modi,6)];
                    % Give sleep effort untill the time 'REST-S' of actigraphy is finished
                    t_temp_error = linspace(tspan(st_fi(day)+diff_modi-1),tspan(st_fi(day+1)-8),20);
                    [t_temp_err,y_temp_err] = ode15s(@(t,V) PCR_shift(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), t_temp_error, V_0);
                    t1_2(diff_modi,1) = t_temp_err(1,1);
                    t1_2(end-8,1) = t_temp_err(end,1);
                    y1_2(diff_modi,:) = y_temp_err(1,:);
                    y1_2(end-8,:) = y_temp_err(end,:);
                elseif t1_2(diff_modi,1) == t1_2(end-8,1)
                    y1_2(diff_modi,:) = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);sleep_effort_Vv;sleep_effort_Vm;y1_2(diff_modi,6)];
                end
                % Stop using Sleep effort when the time 'REST-S' of actigraphy is finished
                V_0 = [y1_2(end-8,1);y1_2(end-8,2);y1_2(end-8,3);y1_2(end-8,4);y1_2(end-8,5);y1_2(end-8,6)];
                [t1_2(end-8:end,1),y1_2(end-8:end,:)] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day+1)-8:st_fi(day+1),1), V_0);
            elseif patt_simul(diff_patt(1)-1) == 'Sleep'
                % Take the value of sleep effort as min of the Vm during
                % such interval
                if t1_2(diff_modi,1) < t1_2(end-9,1)
                    V_0 = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);sleep_effort_Vv;sleep_effort_Vm;y1_2(diff_modi,6)];
                    % Give sleep effort untill the time 'REST-S' of actigraphy is finished
                    [t1_2(diff_modi:end-8,1),y1_2(diff_modi:end-8,:)] = ode15s(@(t,V) PCR_shift(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day)+diff_modi-1:st_fi(day+1)-8,1), V_0);
                elseif (t1_2(diff_modi,1) == t1_2(end-9,1)) 
                    V_0 = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);sleep_effort_Vv;sleep_effort_Vm;y1_2(diff_modi,6)];
                    % Give sleep effort untill the time 'REST-S' of actigraphy is finished
                    t_temp_error = linspace(tspan(st_fi(day)+diff_modi-1),tspan(st_fi(day+1)-8),20);
                    [t_temp_err,y_temp_err] = ode15s(@(t,V) PCR_shift(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), t_temp_error, V_0);
                    t1_2(diff_modi,1) = t_temp_err(1,1);
                    t1_2(end-8,1) = t_temp_err(end,1);
                    y1_2(diff_modi,:) = y_temp_err(1,:);
                    y1_2(end-8,:) = y_temp_err(end,:);
                elseif t1_2(diff_modi,1) == t1_2(end-8,1)
                    y1_2(diff_modi,:) = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);sleep_effort_Vv;sleep_effort_Vm;y1_2(diff_modi,6)];
                end
                % Stop using Sleep effort when the time 'REST-S' of actigraphy is finished
                % When stop using sleep effort, change Vv, Vm as the value
                % at the transition of sleep to wake occur
                V_0 = [y1_2(end-8,1);y1_2(end-8,2);y1_2(end-8,3);y1_2(end-8,4);y1_2(end-8,5);y1_2(end-8,6)];
                [t1_2(end-8:end,1),y1_2(end-8:end,:)] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day+1)-8:st_fi(day+1),1), V_0);
            end
        end
    % When first dismatch occur, model is at the sleep state
    elseif patt_simul(diff_patt(1)) == 'Sleep'
        % If one-time point before dismatching has different sleep/wake
        % state with that of dismatching find the exact time trasition of
        % sleep state occur
        if diff_patt(1) ~= 1
            if patt_simul(diff_patt(1)-1) == 'Wake'
                Qm_temp = ( Q_max ./ (1+exp(-(y1_2(1+(diff_patt(1)-1)*8-7:1+(diff_patt(1)-1)*8,5)-theta)/sigma)) );
                bifur = find(Qm_temp <= 1);
                diff_modi = 1 + 8 * (diff_patt(1)-1) - 8 + bifur(1);
            elseif patt_simul(diff_patt(1)-1) == 'Sleep'
                diff_modi = 1 + 8 * (diff_patt(1)-1);
            end
            % If diff_patt(1) == 1, patt_simul(diff_patt(1)-1) must be 'REST-S'
        else
            diff_modi = 1 + 8 * (diff_patt(1)-1);
        end
        % Calculate Dv which decide whether model is in the bistable region
        % or not
        Dv = -10.2 - (3.37 * 0.5) * ( const + coef_y * y1_2(diff_modi,2) + coef_x * y1_2(diff_modi,1) ) + v_vh * y1_2(diff_modi,6);
        % When model is in the bistable region
        if Dv <= 2.46
            % Change Vv,Vm as the value whenn the bifurcation occur
            V_0 = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);-4.66;-0.12;y1_2(diff_modi,6)];
            [t1_2(diff_modi:end,1),y1_2(diff_modi:end,:)] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day)+diff_modi-1:st_fi(day+1),1), V_0);
        else
            if diff_patt(1) == 1
                % Take the value of wake effort as max (or min) of the Vm during
                % 120th day of First 120 days-initialization
                V_0 = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);wake_effort_Vv;wake_effort_Vm;y1_2(diff_modi,6)];
                [t1_2(diff_modi:end,1),y1_2(diff_modi:end,:)] = ode15s(@(t,V) PCR_shift(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day)+diff_modi-1:st_fi(day+1),1), V_0);
                
                % Give wake effort untill
                % (a) When model enters the bistable region before 'Wake'
                % state of actigraphy data finished
                % (b) When 'Wake' time of actigraphy data finished
                
                % Calculate Dv which decide whether model is in the bistable region
                % or not
                Dv = -10.2 - (3.37 * 0.5) * ( const + coef_y * y1_2(diff_modi:real_sl_time(day)-st_fi(day)+1-8,2) + coef_x * y1_2(diff_modi:real_sl_time(day)-st_fi(day)+1-8,1) ) + v_vh * y1_2(diff_modi:real_sl_time(day)-st_fi(day)+1-8,6);
                % Find interval when model enters the bistable region
                % before 'Wake' state of actigraphy data finished
                circa_align = find(Dv <= 2.46);
                % Case (a)
                if isempty(circa_align) == 0
                    circa_align_start = circa_align(1);
                    % When stop using wake effort, change Vv, Vm as the value
                    % at the transition of wake to sleep occur
                    V_0 = [y1_2(diff_modi+circa_align_start-1,1);y1_2(diff_modi+circa_align_start-1,2);y1_2(diff_modi+circa_align_start-1,3);-4.66;-0.12;y1_2(diff_modi+circa_align_start-1,6)];
                    [t1_2(diff_modi+circa_align_start-1:end,1),y1_2(diff_modi+circa_align_start-1:end,:)] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day)+diff_modi+circa_align_start-2:st_fi(day+1),1), V_0);
                    % Case (b)
                else
                    V_0 = [y1_2(real_sl_time(day)-st_fi(day)+1-8,1);y1_2(real_sl_time(day)-st_fi(day)+1-8,2);y1_2(real_sl_time(day)-st_fi(day)+1-8,3);y1_2(real_sl_time(day)-st_fi(day)+1-8,4);y1_2(real_sl_time(day)-st_fi(day)+1-8,5);y1_2(real_sl_time(day)-st_fi(day)+1-8,6)];
                    [t1_2(real_sl_time(day)-st_fi(day)+1-8:end,1),y1_2(real_sl_time(day)-st_fi(day)+1-8:end,:)] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(real_sl_time(day)-8:st_fi(day+1),1), V_0);
                end
            elseif patt_simul(diff_patt(1)-1) == 'Wake'
                V_0 = [y1_2(diff_modi,1);y1_2(diff_modi,2);y1_2(diff_modi,3);wake_effort_Vv;wake_effort_Vm;y1_2(diff_modi,6)];
                [t1_2(diff_modi:end,1),y1_2(diff_modi:end,:)] = ode15s(@(t,V) PCR_shift(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day)+diff_modi-1:st_fi(day+1),1), V_0);
                
                % Give wake effort untill
                % (a) When model enters the bistable region before 'Wake'
                % state of actigraphy data finished
                % (b) When 'Wake' time of actigraphy data finished
                
                % Calculate Dv which decide whether model is in the bistable region
                % or not
                Dv = -10.2 - (3.37 * 0.5) * ( const + coef_y * y1_2(diff_modi:real_sl_time(day)-st_fi(day)+1-8,2) + coef_x * y1_2(diff_modi:real_sl_time(day)-st_fi(day)+1-8,1) ) + v_vh * y1_2(diff_modi:real_sl_time(day)-st_fi(day)+1-8,6);
                % Find interval when model enters the bistable region
                % before 'Wake' state of actigraphy data finished
                circa_align = find(Dv <= 2.46);
                % Case (a)
                if isempty(circa_align) == 0
                    circa_align_start = circa_align(1);
                    % When stop using wake effort, change Vv, Vm as the value
                    % at the transition of wake to sleep occur
                    V_0 = [y1_2(diff_modi+circa_align_start-1,1);y1_2(diff_modi+circa_align_start-1,2);y1_2(diff_modi+circa_align_start-1,3);y1_2(diff_modi+circa_align_start-1,4);y1_2(diff_modi+circa_align_start-1,5);y1_2(diff_modi+circa_align_start-1,6)];
                    [t1_2(diff_modi+circa_align_start-1:end,1),y1_2(diff_modi+circa_align_start-1:end,:)] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day)+diff_modi+circa_align_start-2:st_fi(day+1),1), V_0);
                    % Case (b)
                else
                    V_0 = [y1_2(real_sl_time(day)-st_fi(day)+1-8,1);y1_2(real_sl_time(day)-st_fi(day)+1-8,2);y1_2(real_sl_time(day)-st_fi(day)+1-8,3);y1_2(real_sl_time(day)-st_fi(day)+1-8,4);y1_2(real_sl_time(day)-st_fi(day)+1-8,5);y1_2(real_sl_time(day)-st_fi(day)+1-8,6)];
                    [t1_2(real_sl_time(day)-st_fi(day)+1-8:end,1),y1_2(real_sl_time(day)-st_fi(day)+1-8:end,:)] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(real_sl_time(day)-8:st_fi(day+1),1), V_0);
                end
            end
        end
    end
end

matching_data = [t1_2,y1_2];

end
