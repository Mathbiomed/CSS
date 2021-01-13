function dVdt = PCR(t,V,it,i, tau_c, mu, v_vh, coef_x, coef_y, const, gate)

Q_max = 100;                                                                % Sleep_wake parameter
theta = 10; sigma = 3;
v_vm = 2.1; v_mv = 1.8; v_vc = 3.37; 
A_m = 1.3; A_v = -10.2; 
tau_m = 10/(60*60); tau_v = 10/(60*60); chi = 45; 
Q_th = 1;
                        
kappa = (12 / pi); i_0 = 9500; p = 0.6; gamma = 0.23;                        % Circadian parameter
k = 0.55; beta = 0.013; f = 0.99669; alpha_0 = 0.16; b = 0.4; 

lambda = 60;

I = interp1(it,i,t);

dVdt = zeros(6,1);
dVdt(1) = (1 / kappa) * ( gamma * ( V(1) - (4 * (V(1).^3)/3) ) - V(2) .* ( ((24 / (f * tau_c))^2) + k * 19.9 * alpha_0 .* (1-V(3)) .* (1-b*V(1)) .* (1-b*V(2)) .* ((heaviside((Q_max ./ (1 + exp(-(V(5)-theta)/sigma)))-Q_th)+(heaviside((Q_max ./ (1 + exp(-(V(5)-theta)/sigma)))-Q_th)-1)*gate*(-0.03)).^p.*(I/i_0).^p) ) ) ;
dVdt(2) = (1 / kappa) * ( V(1) + 19.9 * alpha_0 * (1-V(3)) .* (1-b*V(1)) .* (1-b*V(2)) .* ((heaviside((Q_max ./ (1 + exp(-(V(5)-theta)/sigma)))-Q_th)+(heaviside((Q_max ./ (1 + exp(-(V(5)-theta)/sigma)))-Q_th)-1)*gate*(-0.03)).^p.*(I/i_0)^p));
dVdt(3) = lambda * (alpha_0 * ((heaviside((Q_max ./ (1 + exp(-(V(5)-theta)/sigma)))-Q_th)+(heaviside((Q_max ./ (1 + exp(-(V(5)-theta)/sigma)))-Q_th)-1)*gate*(-0.03)).^p.*(I/i_0)^p) .* (1-V(3)) - beta * V(3) );
dVdt(4) = (1 / tau_v) * ( -V(4) - v_vm * ( Q_max ./ (1+exp(-(V(5)-theta)/sigma)) ) + A_v - (v_vc * (1/2) * (const + coef_y * V(2) + coef_x * V(1) ) ) + (v_vh * V(6)) );
dVdt(5) = (1 / tau_m) * ( -V(5) - v_mv * ( Q_max ./ (1+exp(-(V(4)-theta)/sigma)) ) + A_m);
dVdt(6) = (1 / chi) * (-V(6) + mu * ( Q_max ./ (1+exp(-(V(5)-theta)/sigma)) ) );

end