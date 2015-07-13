function F = control(t,y)

    global control_arr
    global count

    global Mas_u

    flag = false;

    % ��� ������� ���� ������ (10 ���������� � ������� ������� control_arr)
    % ���������� �������� �� ������

    for i=1:count
        t_0_i = control_arr(i,1);
        t_end_i = control_arr(i,2);

        if ((t < t_end_i)&&(t > t_0_i)) || (abs(t - t_0_i) < 1e-6)
            coefs = control_arr(i,3:8);
            c0 = coefs(1);
            c1 = coefs(2);
            c2 = coefs(3);
            c3 = coefs(4);
            c4 = coefs(5);
            c5 = coefs(6);
            d = control_arr(i,9);
            t_0_c = control_arr(i,11);
            t_end_c = control_arr(i,12);

            if (abs(control_arr(i,10) - 2) < 1e-6)   
                d0 = 4.375;
                d1 = -5.25;
                d2 = 1.875;

                Psi_t = c0 + c1/(t_end_c - t_0_c)*(t - t_0_c) + c2/(t_end_c - t_0_c)^2*(t - t_0_c)^2 + ...
                    c3/(t_end_c - t_0_c)^3*(t - t_0_c)^3 +...
                    c4/(t_end_c - t_0_c)^4*(t - t_0_c)^4 + c5/(t_end_c - t_0_c)^5*(t - t_0_c)^5 + ...
                    d*((t - t_0_i)/(t_end_i - t_0_i))^3*(d0 + ...
                    d1*(t - t_0_i)^2/(t_end_i - t_0_i)^2 + d2*(t - t_0_i)^4/(t_end_i - t_0_i)^4);

                dPsi_t = c1/(t_end_c - t_0_c) + 2*c2/(t_end_c - t_0_c)^2*(t - t_0_c) + 3*c3/(t_end_c - t_0_c)^3*(t - t_0_c)^2 + ...
                    4*c4/(t_end_c - t_0_c)^4*(t - t_0_c)^3 + 5*c5/(t_end_c - t_0_c)^5*(t - t_0_c)^4 + ...
                    d*(3*d0*(t - t_0_i)^2/(t_end_i - t_0_i)^3 + 5*d1*(t - t_0_i)^4/(t_end_i - t_0_i)^5 + ...
                    7*d2*(t - t_0_i)^6/(t_end_i - t_0_i)^7);

                ddPsi_t = 2*c2/(t_end_c - t_0_c)^2 + 6*c3/(t_end_c - t_0_c)^3*(t - t_0_c) + ...
                       12*c4/(t_end_c - t_0_c)^4*(t - t_0_c)^2 + 20*c5/(t_end_c - t_0_c)^5*(t - t_0_c)^3 + ... 
                       d*(6*d0*(t - t_0_i)/(t_end_i - t_0_i)^3 + 20*d1*(t - t_0_i)^3/(t_end_i - t_0_i)^5 + ...
                       42*d2*(t - t_0_i)^5/(t_end_i - t_0_i)^7);


            end
            if (abs(control_arr(i,10) - 1) < 1e-6)
                
                Psi_t = c0 + c1/(t_end_c - t_0_c)*(t - t_0_c) + c2/(t_end_c - t_0_c)^2*(t - t_0_c)^2 +...
                    c3/(t_end_c - t_0_c)^3*(t - t_0_c)^3 +...
                    c4/(t_end_c - t_0_c)^4*(t - t_0_c)^4 + c5/(t_end_c - t_0_c)^5*(t - t_0_c)^5 ...
                    + d*(t - t_0_i)^4*(t - t_end_i)^4;
                
                dPsi_t = c1/(t_end_c - t_0_c) + 2*c2/(t_end_c - t_0_c)^2*(t - t_0_c) +...
                    3*c3/(t_end_c - t_0_c)^3*(t - t_0_c)^2 + ...
                    4*c4/(t_end_c - t_0_c)^4*(t - t_0_c)^3 + 5*c5/(t_end_c - t_0_c)^5*(t - t_0_c)^4 ...
                    + 4*d*((t - t_0_i)^3*(t - t_end_i)^4 + (t - t_0_i)^4*(t - t_end_i)^3);
                
                ddPsi_t = 2*c2/(t_end_c - t_0_c)^2 + 6*c3/(t_end_c - t_0_c)^3*(t - t_0_c) + ...
                    12*c4/(t_end_c - t_0_c)^4*(t - t_0_c).^2 + 20*c5/(t_end_c - t_0_c)^5*(t - t_0_c).^3 + ...
                    4*d*(3*(t - t_0_i).^2.*(t - t_end_i).^4 + 4*(t - t_0_i).^3.*(t - t_end_i).^3 + ...
                    4*(t - t_0_i).^3.*(t - t_end_i).^3 + 3*(t - t_0_i).^4.*(t - t_end_i).^2);
            end
            if (abs(control_arr(i,10) - 0) < 1e-6) 
                
                Psi_t = c0 + c1/(t_end_c - t_0_c)*(t - t_0_c) + c2/(t_end_c - t_0_c)^2*(t - t_0_c)^2 + ...
                        c3/(t_end_c - t_0_c)^3*(t - t_0_c)^3 +...
                        c4/(t_end_c - t_0_c)^4*(t - t_0_c)^4 + c5/(t_end_c - t_0_c)^5*(t - t_0_c)^5;

                dPsi_t = c1/(t_end_c - t_0_c) + 2*c2/(t_end_c - t_0_c)^2*(t - t_0_c) + 3*c3/(t_end_c - t_0_c)^3*(t - t_0_c)^2 + ...
                         4*c4/(t_end_c - t_0_c)^4*(t - t_0_c)^3 + 5*c5/(t_end_c - t_0_c)^5*(t - t_0_c)^4;

                ddPsi_t = 2*c2/(t_end_c - t_0_c)^2 + 6*c3/(t_end_c - t_0_c)^3*(t - t_0_c) + ...
                          12*c4/(t_end_c - t_0_c)^4*(t - t_0_c)^2 + 20*c5/(t_end_c - t_0_c)^5*(t - t_0_c)^3;
            end 
            flag = true;
            break;
        end

    end

    % ��������� ����� ?
    if (flag == false)    
        coefs = control_arr(count,3:8);
        c0 = coefs(1);
        c1 = coefs(2);
        c2 = coefs(3);
        c3 = coefs(4);
        c4 = coefs(5);
        c5 = coefs(6);

        t_0_c = control_arr(count,11);
        t_end_c = control_arr(count,12);

        Psi_t = c0 + c1/(t_end_c - t_0_c)*(t - t_0_c) + c2/(t_end_c - t_0_c)^2*(t - t_0_c)^2 + ...
                c3/(t_end_c - t_0_c)^3*(t - t_0_c)^3 +...
                c4/(t_end_c - t_0_c)^4*(t - t_0_c)^4 + c5/(t_end_c - t_0_c)^5*(t - t_0_c)^5;

        dPsi_t = c1/(t_end_c - t_0_c) + 2*c2/(t_end_c - t_0_c)^2*(t - t_0_c) + 3*c3/(t_end_c - t_0_c)^3*(t - t_0_c)^2 + ...
                 4*c4/(t_end_c - t_0_c)^4*(t - t_0_c)^3 + 5*c5/(t_end_c - t_0_c)^5*(t - t_0_c)^4;

        ddPsi_t = 2*c2/(t_end_c - t_0_c)^2 + 6*c3/(t_end_c - t_0_c)^3*(t - t_0_c) + ...
                  12*c4/(t_end_c - t_0_c)^4*(t - t_0_c)^2 + 20*c5/(t_end_c - t_0_c)^5*(t - t_0_c)^3;
    end

    root1 = -0.1; root2 = -0.1;
    c1_st = -(root1 + root2);
    c0_st = root1*root2;

    u = (ddPsi_t - f_kan(y) - c1_st*(y(2) - dPsi_t) - c0_st*(y(1) - Psi_t)) / g_kan(y);


Mas_u = [Mas_u; t u];

F = u;