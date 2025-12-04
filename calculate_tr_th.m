% 计算通信时间阈值的函数
function tr_th = calculate_tr_th(err_r_th, C, V, d, T_s)

    y=qfuncinv(err_r_th);
    Lr_sqrt=(y/log(2)+sqrt((y/log(2)).^2+4*C*d/V))/(2*C/sqrt(V));
    L_r=ceil(Lr_sqrt.^2);
    tr_th = L_r*T_s;

end