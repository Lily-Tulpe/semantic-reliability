%% err_s=1-Q_1与err_s=1-Q对比
% 1个服务器，T_max = 10 
clear;

% 参数初始化
% 感知参数
P_s = 1;
delta = 0.0001;

% 通信参数
P_r = 1;
h_r = 0.5; 
d = 320;

% 计算参数
f = 3; % CPU频率
c = 0.15; % 任务量
a = 3.4955e-1; % 尺度参数
b = -0.0214; % 形状参数 
t_th = 5.7; % 计算阈值
F_W = 0.99;

T_s = 0.025;
sigma2 = 0.01;
T_max = 12;
L_max = floor(T_max/T_s);

% 错误率阈值
err_r_th = 0.01;
err_s_th = 0.01;
err_c_th = 0.01;

h_s_values = normrnd(0, sqrt(sigma2), 1000, 1);
ungamma_s = mean(h_s_values.^2);

[C, V] = calculate_para_c(P_r, h_r, sigma2);

tr_th = calculate_tr_th(err_r_th, C, V, d, T_s);
Lr_th = ceil(tr_th/T_s);

Lc_th = ceil((t_th+c/f)/T_s);

err_r = @(t_r) calculate_err_r(t_r, C, V, d, T_s);

err_c = @(t_c) calculate_err_c(a, b, t_c, c, f, t_th, F_W);

cst1 = @(t) err_r(t(1))- err_r_th;

cst3 = @(t) err_c(t(3)) - err_c_th;

cst_diff1 = @(t) calculate_fun_diff(cst1, t);

cst_diff3 = @(t) calculate_fun_diff(cst3, t);

cst_diff4 = @(t) [1; 1; 1];

err_s_unh_s = @(t_s) calculate_err_s_unh_s(t_s, P_s, ungamma_s, sigma2, delta, T_s);

obj_unh_s = @(t) err_r(t(1)) + err_s_unh_s(t(2)) + err_c(t(3));

cst2_unh_s = @(t) err_s_unh_s(t(2)) - err_s_th;

obj_diff_unh_s =@(t) calculate_fun_diff(obj_unh_s, t);

cst_diff2_unh_s = @(t) calculate_fun_diff(cst2_unh_s, t);

cst_diff_unh_s = {cst_diff1, cst_diff2_unh_s, cst_diff3, cst_diff4};

obj_best = NaN(2, length(h_s_values));
obj_best_h_s = zeros(2, 1);
obj_best_unh_s = zeros(2, 1);

x_ep = zeros(2, 3);
x_gl = zeros(2, 3);

for k = 1:1000

    h_s = h_s_values(k);  

    ts_th = calculate_ts_th(err_s_th, P_s, h_s, sigma2, delta, T_s);
    Ls_th = ceil(ts_th/T_s);

    if Ls_th + Lr_th + Lc_th >= L_max
        continue;
    end

    err_s = @(t_s) calculate_err_s(t_s, P_s, h_s, sigma2, delta, T_s);

    obj = @(t) err_r(t(1)) + err_s(t(2)) + err_c(t(3));

    cst2 = @(t) err_s(t(2)) - err_s_th;

    obj_diff =@(t) calculate_fun_diff(obj, t);

    cst_diff2 = @(t) calculate_fun_diff(cst2, t);

    cst_diff = {cst_diff1, cst_diff2, cst_diff3, cst_diff4};

    cst4 = @(t) sum(t)-T_max;

    cst = {cst1, cst2, cst3, cst4};

    max_abs = [T_max/3; T_max/3; T_max/3];
    opt_err = 1e-5;

    [x_ep(1, :), obj_best(1, k)] = ellipsoid_method(obj, obj_diff, cst, cst_diff, max_abs, opt_err, Ls_th, L_max, Lr_th, Lc_th, T_s, T_max);
    
    [x_gl(1, :), obj_best(2, k)] = global_search(Lr_th, Ls_th, Lc_th, T_s, T_max, obj);
  
end

obj_temp_ep = obj_best(1, :);
obj_temp_gl = obj_best(2, :);

obj_best_h_s(1, 1) = mean(obj_temp_ep, 'omitnan');
obj_best_h_s(2, 1) = mean(obj_temp_gl, 'omitnan');

cst_unh_s = {cst1, cst2_unh_s, cst3, cst4};

[ts_th_unh_s, ~, err_s_bi] = bisect(err_s_unh_s, 0, T_max, T_s, err_s_th);        
Ls_th_unh_s = ceil(ts_th_unh_s/T_s);

[x_ep(2, :), obj_best_unh_s(1, 1)] = ellipsoid_method(obj_unh_s, obj_diff_unh_s, cst_unh_s, cst_diff_unh_s, max_abs, opt_err, Ls_th_unh_s, L_max, Lr_th, Lc_th, T_s, T_max);

[x_gl(2, :), obj_best_unh_s(2, 1)] = global_search(Lr_th, Ls_th_unh_s, Lc_th, T_s, T_max, obj_unh_s);

disp('h_s is known:')
disp(obj_best_h_s);
disp('h_s is unknown:')
disp(obj_best_unh_s);

%%
function [x_best, obj_best] = ellipsoid_method(obj, obj_diff, cst, cst_diff, max_abs, opt_err, Ls_th, L_max, Lr_th, Lc_th, T_s, T_max)

    x_best = [];
    obj_best = inf;

    if Ls_th + Lr_th + Lc_th >= L_max
        return;
    end

    for j = 1:1000
    
        flag = 0;
    
        ts_th_rand = T_s*randi([Ls_th, L_max-Lr_th-Lc_th], 1, 1);
        tr_th_rand = T_s*randi([Lr_th, L_max-Ls_th-Lc_th], 1, 1);
        
        x0 = [tr_th_rand; ts_th_rand; T_max-tr_th_rand-ts_th_rand];
        
        for opt_i = 1 : length(cst)
            if cst{opt_i}(x0) >= 1e-10
                flag = 1;
                break; 
            end
        end

        if flag
            continue;
        end
    
        try
            [obj_temp, x_temp, ~] = ellipsoid_optimize(obj, obj_diff, cst, cst_diff, x0, max_abs, opt_err);
        catch ME
            % 如果优化器出错，记录并继续
    %         warning('ellipsoid_optimize error at iter %d: %s', i, ME.message);
            continue;
        end
    
        if obj_temp < obj_best
            obj_best = obj_temp;
            x_best = x_temp;
        end
    
    end

end

function [x_best, obj_best] = global_search(Lr_th, Ls_th, Lc_th, T_s, T_max, obj)

    obj_best = inf;
    for t_r = Lr_th*T_s:T_s:T_max-(Ls_th+Lc_th)*T_s
        for t_s = Ls_th*T_s:T_s:T_max-Lc_th*T_s-t_r
            t_c = T_max-t_r-t_s;             
            t = [t_r; t_s; t_c];
            obj_value = obj(t);
            if obj_value <= obj_best
                obj_best = obj_value;
                x_best = t;
            end      
        end
    end

end