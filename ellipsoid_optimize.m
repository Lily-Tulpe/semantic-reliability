function [obj_best, x_best, opt_loopi] = ellipsoid_optimize(obj, obj_diff, cst, cst_diff, x0, max_abs, opt_err)
%% Optimization tool based on ellipsoid method (Y)
% Especially for problems with complex form so that cvx cannot recognize.

% 1. Explanations:
% 1.1 INPUT:
%        obj            --  function handle of objective function (to be minimized)
%        obj_diff       --  function handle of derivative of objective function
%        cst, cst_diff  --  cell structure including function handles for constraints
%        x0             --  initial (feasible) input
%        max_abs        --  maximum distance in each dimension from optimal point to initial point x0
%        opt_err        --  maximum tolenrant error to the optimal
% 1.2 OUTPUT:
%        obj_best       --  optimal value (minimized)
%        x_best         --  optimal point x

% default output for error occurance

obj_best = []; % 定义的输出的最优目标值：一个值
x_best = []; % 定义输出的optimal point，得到optiaml value时的x取值

% 判断值的输入是否正确，如果有多个变量，那么这里是不是应该需要多个判断n个情况，比如obj中有n个变量，则需要判断
% size(obj,1)，size(obj,2)，size(obj,3)，...，size(obj,n)均满足 ~= 1的情况。
% 那么是否可以改写成这种形式，通过循环的方式来判断多个变量
% n = 10;
% for i=1:1:n
%    if size(obj,i) ~= 1 || size(obj_diff,i) ~= 1 || size(cst,i) ~= 1 || size(cst_diff,1) ~= 1
%        disp('Wrong objective format or wrong constraint format.')
%        return;
%    end
% end
if size(obj,1) ~= 1 || size(obj,2) ~= 1 || size(obj_diff,1) ~= 1 || size(obj_diff,2) ~= 1 ...
        || size(cst,1) ~= 1 || size(cst_diff,1) ~= 1 || size(cst,2) ~= size(cst_diff,2)
    disp('Wrong objective format or wrong constraint format.')
    return;
end

% reshape的作用：对矩阵进行处理，语法是reshape(A,m,n)或reshape(A,[m,n])，意思是将A的行列排成m行n列
% 另外，reshape是按照列来取数据的
max_abs = reshape( max_abs, size(max_abs,1) * size(max_abs,2), 1); % 将max_abs重新排列成size(max_abs,1) * size(max_abs,2)行，1列的形式
x0 = reshape( x0, size(x0,1) * size(x0,2), 1); % 将x0重新排列成size(x0,1) * size(x0,2)行1列的形式
% 这里同样是检查数据的输入是否符合规范
% 1、max_abs的数据形式应该和x0一致；
% 2、opt_err应该大于0；
if length(max_abs) ~= length(x0)
    disp('Wrong initial input or wrong range.');
    return;
elseif opt_err <= 0
    disp('Negative threshold error.');
    return;
end

% size是用来求矩阵大小的，size(A,n)
% 当n为1时，返回矩阵A的行数；当n为2时，返回矩阵A的列数
opt_cst_nr = size(cst, 2); % 返回矩阵cst的列数

opt_n = length(x0); % 向量x0的长度
opt_P = diag(opt_n * max_abs .^2); % diag：将矩阵对角化
opt_loopi = 0;
opt_objbest = inf;
opt_x = x0; % 记录每一次x(1)和x(2)的值
%opt_err=1e-6;

% Strictly feasible检查，即判断是否满足约束条件
for opt_i = 1 : opt_cst_nr
    if cst{opt_i}(opt_x) >= 1e-10
        disp('Initial point is not strictly feasible.');
        return;
    end
end

% 开始迭代，利用椭圆法进行求解
while 1 %optloop_i<5000 %%% note that for a large N a large number of iterations may be required.
    
    % mark the iteration -- objective iteration or constraint iteration
    opt_mark = 1;
    
    % 这里的情况相当于算法中x(k) is infeasible的情况
    % infeasible时，利用ellipsoid来根据constraint计算alpha
    for opt_i = 1 : opt_cst_nr
        if cst{opt_i}(opt_x) > 0
            opt_g = cst_diff{opt_i}(opt_x);
            opt_g = reshape(opt_g, opt_n, 1);
            opt_gn = opt_g / (opt_g' * opt_P * opt_g) ^0.5;
            opt_a = cst{opt_i}(opt_x) / (opt_g' * opt_P * opt_g) ^0.5;
            opt_mark = 0;
            break;
        end
    end
    
    % feasible 根据deep-cut ellipsoid来计算object
    % objective iteration
    if opt_mark == 1
        opt_obj = obj(opt_x);
        opt_g = obj_diff(opt_x);
        opt_g = reshape(opt_g, opt_n, 1);
        if opt_obj < opt_objbest
            opt_objbest = opt_obj;
            x_best = opt_x;
        end
        opt_gn = opt_g / (opt_g' * opt_P * opt_g) ^0.5;
        opt_a = (opt_obj - opt_objbest) / (opt_g' * opt_P * opt_g) ^0.5;
    end
    
    % stopping iteration checking
    if (opt_g' * opt_P * opt_g) ^0.5 <= opt_err && opt_mark == 1
%         disp( '椭球足够小!' )
        break; 
    end
    if opt_a > 1 && opt_mark == 0 
%         disp( '非凸!' )
        break;
%         opt_a = 0;
    end
    
    % 出现这种情况时，表示原问题不是一个凸问题
    if opt_a > 1
%         disp( 'Error! Convexity is not guaranteed!' )
       break;
    end
    
    % updating
    % 更新ellipsoid center和ellipsoid shape
    % opt_x表示ellipsoid center；opt_P表示ellipsoid shape
    opt_x = opt_x - (1 + opt_n * opt_a) / (opt_n + 1) * opt_P * opt_gn;
%     opt_x = max(opt_x, [1.875; 0.65; 5.75; 5.75; 5.75]);
    opt_P = opt_n ^2 / (opt_n ^2 - 1) * (1 - opt_a ^2) * (opt_P - 2 * (1 + opt_n * opt_a) / (opt_n + 1) / (opt_a + 1) * opt_P * (opt_gn * opt_gn') * opt_P);
    
    % opt_loopi记录迭代的次数
    opt_loopi = opt_loopi + 1;
end

obj_best = opt_objbest;

end

