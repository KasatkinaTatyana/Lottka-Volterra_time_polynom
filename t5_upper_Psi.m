function [F, params, coefs] = t5_upper_Psi(Y_0, Y_end, t_0, t_end, dt, constr2)

c_t = t5_calc_coefs(Y_0, Y_end, t_0, t_end);
c0 = c_t(1);
c1 = c_t(2);
c2 = c_t(3);
c3 = c_t(4);
c4 = c_t(5);
c5 = c_t(6);

t = t_0:dt:t_end;

Psi = c0 + c1/(t_end - t_0)*(t - t_0) + c2/(t_end - t_0)^2*(t - t_0).^2 + c3/(t_end - t_0)^3*(t - t_0).^3 + ...
      c4/(t_end - t_0)^4*(t - t_0).^4 + c5/(t_end - t_0)^5*(t - t_0).^5;

N_shift = 20;

N = (t_end - t_0)/dt + 1;

flag = 0;
% Ищем границы интервала, на котором нарушается ограничение
for i=1:N
    % ------- -- y > constr -----------------------------------
    if ((flag == 0)&&(Psi(i) >= constr2))
    % -----------------------------------------------------------------
        if ((i - N_shift) < 1)
            t_left = t_0;
        else
            t_left = t(i - N_shift);
        end
        flag = 1;
    end
    % --------- y < constr ------------------------------------
    if ((flag==1)&&(Psi(i) < constr2))
    % -----------------------------------------------------------------
        if ((i + N_shift) > N + 1)
            t_right = t_end;
        else
            t_right = t(i + N_shift);
        end
        break;
    end
end

t_1 = t_left:dt:t_right;     
N_d = 1000;
d_min = -1e7;
d_max = -1e6;
h_d = (d_min - d_max) / N_d;
for d = d_max : h_d : d_min
    Psi_1 = c0 + c1/(t_end - t_0)*(t_1 - t_0) + c2/(t_end - t_0)^2*(t_1 - t_0).^2 + c3/(t_end - t_0)^3*(t_1 - t_0).^3 +...
            c4/(t_end - t_0)^4*(t_1 - t_0).^4 + c5/(t_end - t_0)^5*(t_1 - t_0).^5 ...
            + d*(t_1 - t_left).^4.*(t_1 - t_right).^4;
            
    if (max(Psi_1) < constr2)
        break;
    end
end
dPsi_1 = c1/(t_end - t_0) + 2*c2/(t_end - t_0)^2*(t_1 - t_0) + 3*c3/(t_end - t_0)^3*(t_1 - t_0).^2 + ...
    4*c4/(t_end - t_0)^4*(t_1 - t_0).^3 + 5*c5/(t_end - t_0)^5*(t_1 - t_0).^4 ...
    + 4*d*((t_1 - t_left).^3.*(t_1 - t_right).^4 + (t_1 - t_left).^4.*(t_1 - t_right).^3);

ddPsi_1 = 2*c2/(t_end - t_0)^2 + 6*c3/(t_end - t_0)^3*(t_1 - t_0) + ...
          12*c4/(t_end - t_0)^4*(t_1 - t_0).^2 + 20*c5/(t_end - t_0)^5*(t_1 - t_0).^3 + ...
          4*d*(3*(t_1 - t_left).^2.*(t_1 - t_right).^4 + 4*(t_1 - t_left).^3.*(t_1 - t_right).^3 + ...
          4*(t_1 - t_left).^3.*(t_1 - t_right).^3 + 3*(t_1 - t_left).^4.*(t_1 - t_right).^2);

% Y_1_end = [Psi_1(end) dPsi_1(end) ddPsi_1(end)];
% 
% replace_part = t5_trajectory_synthesis(Y_1_end, Y_end, t_right, t_end, dt);
% coefs_2 = t5_calc_coefs(Y_1_end, Y_end, t_right, t_end);
% t_2 = t_right:dt:t_end;
% 
% Psi_2 = replace_part(1,:);
% dPsi_2 = replace_part(2,:);
% ddPsi_2 = 2*coefs_2(3)/(t_2(end)- t_2(1))^2 + 6*coefs_2(4)/(t_2(end) - t_2(1))^3*(t_2 - t_2(1)) + ...
%           12*coefs_2(5)/(t_2(end) - t_2(1))^4*(t_2 - t_2(1)).^2 + 20*coefs_2(6)/(t_2(end) - t_2(1))^5*(t_2 - t_2(1)).^3;

F = [t_1; Psi_1; dPsi_1; ddPsi_1];
% F2 = [t_2; Psi_2; dPsi_2; ddPsi_2];

params = [d; t_0; t_end; 0];

coefs = [c0 c1 c2 c3 c4 c5; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];