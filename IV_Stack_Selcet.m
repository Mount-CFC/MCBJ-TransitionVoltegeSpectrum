%找出有突跃的traces并叠加
% 记录AO1 bias，AI0 sampling Voltage, current, Conductance LogG
% 最后计算TVS谱

% clc
clear 
close all
tic

%%%%   Parameters setting %%%%%
scan_voltge = 1;   %voltage scan range
% cur_density = 1;

min_cur = -1.5;      % Min limit of y axis in I-V plot
max_cur = 2.5;       % Max limit of y axis in I-V plot
min_logG = -4.5;     % Min limit of y axis in LogG-V plot
max_logG = -2.5;     % Max limit of y axis in LogG-V plot

%%%%%%% Selection Setting %%%%%%%%
%Average conductance should HIGHER than this value between -0.9~-0.8V
high_conductance = -3.5; 
%Average conductance should LOWER than this value between -0.3~-0.2V
low_conductance = 555;

%%%%%%% Plot Selection %%%%%%%%
plot_select = 0;   % If = 1, plot selected IV and LogG-V single traces



[filename,filepath]=uigetfile('*.tdms','Select data files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else 
    filename1{1}=filename;
end

num_files = length(filename1);
fprintf('Num of file(s): %d\n', num_files)
%读取电压、电流、电导
ForwardTraceBias = [];
ForwardTraceCurrent = [];
ForwardTraceLogG = [];
ReverseTraceBias = [];
ReverseTraceCurrent = [];
ReverseTraceLogG = [];

% 筛选
ForwardCurrent_slct = [];
ForwardBias_slct = [];
ReverseCurrent_slct = [];
ReverseBias_slct = [];

for n = 1:num_files
    struc=TDMS_readTDMSFile(filename1{n});
    data_bias=struc.data{1,3};              %第一行第3列，提取Bias
%     data_samplingV=struc.data{1,4}; %第一行第4列，提取sampling voltage
    data_Cur = struc.data{1,5};             % 第一行第5列，提取current
    data_logG = struc.data{1,6};            % 第一行第6列，提取log (G/G0)
    [ForwardTraceBias_temp,...
        ForwardTraceCurrent_temp,...
        ForwardTraceLogG_temp,...
        ReverseTraceBias_temp,...
        ReverseTraceCurrent_temp,...
        ReverseTraceLogG_temp] = CutIV(data_bias, data_Cur, data_logG, scan_voltge);

    ForwardTraceBias = [ForwardTraceBias ForwardTraceBias_temp];%在原有元胞后面添加新的元胞
    ForwardTraceCurrent = [ForwardTraceCurrent ForwardTraceCurrent_temp];
    ForwardTraceLogG = [ForwardTraceLogG ForwardTraceLogG_temp];
    ReverseTraceBias = [ReverseTraceBias ReverseTraceBias_temp];
    ReverseTraceCurrent = [ReverseTraceCurrent ReverseTraceCurrent_temp];
    ReverseTraceLogG = [ReverseTraceLogG ReverseTraceLogG_temp];
%     筛选函数
    [ForwardBias_slct_temp, ForwardCurrent_slct_temp] = SelectIV_diffANDhigh(ForwardTraceBias_temp,...
        ForwardTraceCurrent_temp,...
        low_conductance,...
        high_conductance);
    [ReverseBias_slct_temp, ReverseCurrent_slct_temp] = SelectIV_diffANDhigh(ReverseTraceBias_temp,...
        ReverseTraceCurrent_temp,...
        low_conductance,...
        high_conductance);
%     筛选后的电流电压
    ForwardCurrent_slct = [ForwardCurrent_slct  ForwardCurrent_slct_temp];
    ForwardBias_slct = [ForwardBias_slct  ForwardBias_slct_temp];
    ReverseCurrent_slct = [ReverseCurrent_slct ReverseCurrent_slct_temp];
    ReverseBias_slct = [ReverseBias_slct ReverseBias_slct_temp];
    clear ForwardTraceBias_temp ForwardTraceCurrent_temp ForwardTraceLogG_temp ReverseTraceBias_temp ReverseTraceCurrent_temp ReverseTraceLogG_temp

    clear ForwardBias_slct_temp ForwardBias_slct_temp



end

for i = 1:length(filename1)
    fprintf('File:%s\n', filename1{i})
end

% transform to current density
% for i = 1:length(ForwardTraceCurrent)
%     ForwardTraceCurrent{i} = ForwardTraceCurrent{i} ./ cur_density;
% end
% for j = 1:length(ReverseTraceCurrent)
%     ReverseTraceCurrent{j} = ReverseTraceCurrent{j} ./ cur_density;
% end
% for k = 1:length(ForwardCurrent_slct)
%     ForwardCurrent_slct{k} = ForwardCurrent_slct{k} ./ cur_density;
% end
% for m = 1:length(ReverseCurrent_slct)
%     ReverseCurrent_slct{m} = ReverseCurrent_slct{m} ./ cur_density;
% end


%% NO Selection
%I-V f+r
figure(1)
hist_all = plot_IV([ForwardTraceBias ReverseTraceBias], [ForwardTraceCurrent,ReverseTraceCurrent], -scan_voltge,scan_voltge,min_cur,max_cur,150,150);
title('F and R', 'Interpreter', 'tex','FontSize',15)
hold on
[XFittedAll,YFittedAll] = master_curve(hist_all,-scan_voltge,scan_voltge,min_cur,max_cur,150,150);%值与上面函数输入一致
plot(XFittedAll,YFittedAll,'-b','linewidth',1);
%%% 用于数据保存 %%%
% % 保存二维图
% save('HistIV_FR.txt','hist_all','-ascii')
% % type('Histcount.txt')
% % 方便导出拟合曲线，同时删除中间异常部分
XFittedIV_FR_part = XFittedAll(XFittedAll<-0.1 | XFittedAll>0.1)';
YFittedIV_FR_part = YFittedAll(XFittedAll<-0.1 | XFittedAll>0.1)';

%LogG-V f+r
figure(2)
hist = plot_IV([ForwardTraceBias ReverseTraceBias], [ForwardTraceLogG,ReverseTraceLogG], -scan_voltge,scan_voltge,min_logG,max_logG,150,150);
title('logG - F&R', 'Interpreter', 'tex','FontSize',15)
ylabel('Conductance / LogG', 'Interpreter', 'tex','FontSize',15)
hold on
[XFittedLogG,YFittedLogG] = master_curve(hist,-scan_voltge,scan_voltge,min_logG,max_logG,150,150);
plot(XFittedLogG,YFittedLogG,'-b','linewidth',1);
%%% 用于数据保存 %%%
% 保存二维图
% save('HistGV.txt','hist','-ascii')
% type('Histcount.txt')
% 方便导出拟合曲线，同时删除中间异常部分
% XFittedLogG_part = XFittedLogG(XFittedLogG<-0.1 | XFittedLogG>0.1)';
% YFittedLogG_part = YFittedLogG(XFittedLogG<-0.1 | XFittedLogG>0.1)';

%I-V Reverse
figure(3)
hist = plot_IV(ReverseTraceBias, ReverseTraceCurrent, -scan_voltge,scan_voltge,min_cur,max_cur,150,150);
title('Reverse', 'Interpreter', 'tex','FontSize',15)
hold on
[XFittedReverse,YFittedReverse] = master_curve(hist,-scan_voltge,scan_voltge,min_cur,max_cur,150,150);%值与上面函数输入一致
plot(XFittedReverse,YFittedReverse,'-b','linewidth',1);

%I-V Forward
figure(4)
hist = plot_IV(ForwardTraceBias , ForwardTraceCurrent, -scan_voltge,scan_voltge,min_cur,max_cur,150,150);
title('Forward', 'Interpreter', 'tex','FontSize',15)
hold on
[XFittedForward,YFittedForward] = master_curve(hist,-scan_voltge,scan_voltge,min_cur,max_cur,150,150);%值与上面函数输入一致
plot(XFittedForward,YFittedForward,'-b','linewidth',1);


%汇总
figure(5)
% 可绘制区域的位置和大小，指定为 [left bottom width height] 形式的向量。
% 此区域不包括图窗边框、标题栏、菜单栏和工具栏。
% help-figure属性
set(gcf, 'unit','centimeters','Position', [5,5,15,10])
subplot(2,2,1)
hist = plot_IV(ForwardTraceBias , ForwardTraceCurrent, -scan_voltge,scan_voltge, min_cur,max_cur,150,150);
title('Forward', 'Interpreter', 'tex','FontSize',15)
hold on
plot(XFittedForward,YFittedForward,'-b','linewidth',1);

subplot(2,2,2)
hist = plot_IV(ReverseTraceBias, ReverseTraceCurrent, -scan_voltge,scan_voltge,min_cur,max_cur,150,150);
title('Reverse', 'Interpreter', 'tex','FontSize',15)
hold on
plot(XFittedReverse,YFittedReverse,'-b','linewidth',1);


subplot(223)
hist = plot_IV([ForwardTraceBias ReverseTraceBias], [ForwardTraceCurrent,ReverseTraceCurrent], -scan_voltge,scan_voltge,min_cur,max_cur,150,150);
title('F and R', 'Interpreter', 'tex','FontSize',15)
hold on
plot(XFittedAll,YFittedAll,'-b','linewidth',1);

subplot(224)
hist = plot_IV([ForwardTraceBias ReverseTraceBias], [ForwardTraceLogG,ReverseTraceLogG], -scan_voltge,scan_voltge,min_logG,max_logG,150,150);
title('logG - F&R', 'Interpreter', 'tex','FontSize',15)
ylabel('Conductance / LogG', 'Interpreter', 'tex','FontSize',15)
hold on
plot(XFittedLogG,YFittedLogG,'-b','linewidth',1);


%% After selection
% 筛选后的曲线的热力图：正扫+反扫的子图
figure(6)
set(gcf, 'unit','centimeters','Position', [5,5,20,10])
%Forward
subplot(131)
hist = plot_IV(ForwardBias_slct , ForwardCurrent_slct, -scan_voltge,scan_voltge,min_cur,max_cur,150,150);
title('Forward selcet', 'Interpreter', 'tex','FontSize',15)
% hold on
% [XFittedSltF,YFittedSltF] = master_curve(hist,-1,1,-2,2,150,150);%值与上面函数输入一致
% plot(XFittedSltF,YFittedSltF,'-b','linewidth',1);
%画代表性的单条曲线
% hold on
% plot(ForwardBias_slct{7},ForwardCurrent_slct{7},'b')


%Reverse
subplot(132)
plot_IV(ReverseBias_slct , ReverseCurrent_slct, -scan_voltge,scan_voltge,min_cur,max_cur,150,150);
title('Reverse selcet', 'Interpreter', 'tex','FontSize',15)
% hold on
% [XFittedSltR,YFittedSltR] = master_curve(hist,-1,1,-2,2,150,150);%值与上面函数输入一致
% plot(XFittedSltR,YFittedSltR,'-b','linewidth',1);
% 画代表性的单条曲线
% hold on
% plot3 = plot(ReverseBias_slct{25},ReverseCurrent_slct{25},'b');
% plot3.Color(4) = 0.3;       %调节Color(4)这个参数可以设置不同的透明度


%Reverse and forward
subplot(133)
hist_IV_FR_slt = plot_IV([ForwardBias_slct ReverseBias_slct], [ForwardCurrent_slct ReverseCurrent_slct], -scan_voltge,scan_voltge,min_cur,max_cur,150,150);
title('F+R selcet', 'Interpreter', 'tex','FontSize',15)
hold on
[XFittedAll_slt,YFittedAll_slt] = master_curve(hist_IV_FR_slt,-scan_voltge,scan_voltge,min_cur,max_cur,150,150);%值与上面函数输入一致
plot(XFittedAll_slt,YFittedAll_slt,'-b','linewidth',1);
%%% 用于数据保存 %%%
% 保存二维图
% save('Hist_IV_FR_slt.txt','hist_IV_FR_slt','-ascii')
% type('Histcount.txt')
% IV 曲线还算正常，不需要筛选


fprintf(['Num of Forward trace:%d\t' ...
    'Selected:%d\n' ...
    'Num of Reverse trace:%d\t' ...
    'Selected:%d\t\n'], ...
    length(ForwardTraceCurrent), length(ForwardCurrent_slct), length(ReverseTraceCurrent), length(ReverseCurrent_slct))


%% 查看筛选后曲线的单条
if plot_select == 1
    k=1;
    row = 5;
    col = 8;
    for i=1:length(ForwardCurrent_slct)
        if k > row*col
            break
        end
        figure(20)
        t = subplot(row, col, k);
    %     plot(ForwardBias_slct{i},ForwardCurrent_slct{i})
    %     plot(ForwardCurrent_slct{i})
        plot(ForwardBias_slct{i}, ForwardCurrent_slct{i})
        ylim([min_cur,max_cur])
        title(num2str(i))
        k = k+1;

    end

    k=1;
    row = 5;
    col = 8;
    StartTrace = 1; %开始的条数序号
    % CurReal = (10 .^ currentSelected{i}) .* 1e-9;
    % logG = log10((CurReal ./ biasSelected{i}) ./ (77.6e-6));
    for i= StartTrace:length(ForwardBias_slct)
        if k > row*col
            break
        end
        CurReal = (10 .^ ForwardCurrent_slct{i}) .* 1e-9;
        logG = log10((CurReal ./ abs(ForwardBias_slct{i})) ./ (77.6e-6));
        figure(21)
        t = subplot(row, col, k);
    %     plot(ForwardBias_slct{i},ForwardCurrent_slct{i})
    %     plot(ForwardCurrent_slct{i})
        plot(ForwardBias_slct{i}, logG)
        ylim([min_logG,max_logG])
        title(num2str(i))
        k = k+1;

    %     clear CurReal logG

    end
end
%% TVS
% XFittedForward, YFittedForward 正扫拟合

TVSrange=0.1;
TVSrangetrans=1/TVSrange; %10
X_TVS=1 ./ XFittedForward; %1/V
% Y_TVS=log(abs(YfittedReal).* X_TVS .* X_TVS);
Y_TVS = log((10 .^ YFittedForward) .* X_TVS .* X_TVS);
for i=1:length(X_TVS)-1
    if X_TVS(i) >= -TVSrangetrans && X_TVS(i+1) < -TVSrangetrans
        sp=i;
    end
    if X_TVS(i) >= TVSrangetrans && X_TVS(i+1) < TVSrangetrans
        ep=i;
    end
end
figure(31);
set(gcf, 'unit','centimeters','Position', [10,2,17,17])
subplot(211)
% plot(X_TVS(1:sp),Y_TVS(1:sp),'-k','linewidth',1.3);
scatter(X_TVS(1:sp),Y_TVS(1:sp),'filled','black','s');
ylabel({'ln(\itI/V^2)'},'Interpreter','tex','FontSize',15);
xlabel({'1/V (V^{-1})'},'Interpreter','tex','FontSize',15);
title('Forward TVS','FontSize',15)
hold on;
% plot(X_TVS(ep:end),Y_TVS(ep:end),'-k','linewidth',1.3);
scatter(X_TVS(ep:end),Y_TVS(ep:end),'filled','black','s');

%反扫TVS
X_TVS_R = 1 ./ XFittedReverse;
Y_TVS_R = log((10 .^ YFittedReverse) .* X_TVS_R .* X_TVS_R);
for k=1:length(X_TVS)-1
    if X_TVS(k) >= -TVSrangetrans && X_TVS(k+1) < -TVSrangetrans
        sp_R=k;
    end
    if X_TVS(k) >= TVSrangetrans && X_TVS(k+1) < TVSrangetrans
        ep_R=k;
    end
end
subplot(212)
scatter(X_TVS_R(1:sp_R),Y_TVS_R(1:sp_R),'filled','black','s');
ylabel({'ln(\itI/V^2)'},'Interpreter','tex','FontSize',15);
xlabel({'1/V (V^{-1})'},'Interpreter','tex','FontSize',15);
title('Reverse TVS','FontSize',15)
hold on;
scatter(X_TVS_R(ep_R:end),Y_TVS_R(ep_R:end),'filled','black','s');


%筛选TVS
figure(55)
X_TVS_FR_slt = 1 ./ XFittedAll_slt;
Y_TVS_FR_slt = log((10 .^ YFittedAll_slt) .* X_TVS_FR_slt .* X_TVS_FR_slt);
for k=1:length(X_TVS)-1
    if X_TVS(k) >= -TVSrangetrans && X_TVS(k+1) < -TVSrangetrans
        sp_R=k;
    end
    if X_TVS(k) >= TVSrangetrans && X_TVS(k+1) < TVSrangetrans
        ep_R=k;
    end
end
scatter(X_TVS_FR_slt(1:sp_R),Y_TVS_FR_slt(1:sp_R),'filled','black','s');
ylabel({'ln(\itI/V^2)'},'Interpreter','tex','FontSize',15);
xlabel({'1/V (V^{-1})'},'Interpreter','tex','FontSize',15);
title('FR slt TVS','FontSize',15)
hold on;
scatter(X_TVS_FR_slt(ep_R:end),Y_TVS_FR_slt(ep_R:end),'filled','black','s');

%% 方便复制导出TVS数据
TVS_x_reverse_positive_bias = X_TVS_R(ep_R:end)';
TVS_y_reverse_positive_bias = Y_TVS_R(ep_R:end)';
TVS_x_reverse_negative_bias = X_TVS_R(1:sp_R)';
TVS_y_reverse_negative_bias = Y_TVS_R(1:sp_R)';

TVS_x_fr_slt_positive_bias = X_TVS_FR_slt(ep_R:end)';
TVS_y_fr_slt_positive_bias = Y_TVS_FR_slt(ep_R:end)';
TVS_x_fr_slt_negative_bias = X_TVS_FR_slt(1:sp_R)';
TVS_y_fr_slt_negative_bias = Y_TVS_FR_slt(1:sp_R)';



toc

%% 淳哥的TVS代码，结果一样的
% TVSrange=0.1;
% TVSrangetrans=1/TVSrange;
% X_TVS=1./XFittedForward;
% YfittedReal = 10.^ YFittedForward;
% Y_TVS=log(abs(YfittedReal).*X_TVS.*X_TVS);
% for i=1:length(X_TVS)-1
%     if X_TVS(i)>=-TVSrangetrans & X_TVS(i+1)<-TVSrangetrans
%         sp=i;
%     end
%     if X_TVS(i)>=TVSrangetrans & X_TVS(i+1)<TVSrangetrans
%         ep=i;
%     end
% end
% figure(50);
% plot(X_TVS(1:sp),Y_TVS(1:sp),'-k','linewidth',1.3);
% ylabel({'ln(I/V^2)'},'Interpreter','latex','FontSize',15);
% xlabel({'1/V (V^-1)'},'Interpreter','latex','FontSize',15);
% hold on;
% plot(X_TVS(ep:end),Y_TVS(ep:end),'-k','linewidth',1.3);


%% 筛选后的logG画图
logG_F = [];
logG_R = [];
for i=1:length(ForwardCurrent_slct)
    CurReal = (10 .^ ForwardCurrent_slct{i}) .* 1e-9;
    logG_F{i} = log10((CurReal ./ abs(ForwardBias_slct{i})) ./ (77.6e-6));
end
for i=1:length(ReverseCurrent_slct)
    CurReal = (10 .^ ReverseCurrent_slct{i}) .* 1e-9;
    logG_R{i} = log10((CurReal ./ abs(ReverseBias_slct{i})) ./ (77.6e-6));
end

%LogG-V f+r
figure(32)
hist_logG_slt = plot_IV([ForwardBias_slct ReverseBias_slct], [logG_F,logG_R], -scan_voltge,scan_voltge,min_logG,max_logG,150,150);
title('logG - selected', 'Interpreter', 'tex','FontSize',15)
ylabel('Conductance / LogG', 'Interpreter', 'tex','FontSize',15)
hold on
[XFittedLogG,YFittedLogG] = master_curve(hist_logG_slt,-scan_voltge,scan_voltge,min_logG,max_logG,150,150);
plot(XFittedLogG,YFittedLogG,'-b','linewidth',1);
%%% 用于数据保存 %%%
% 保存二维图
% save('HistGV.txt','hist','-ascii')
% type('Histcount.txt')
% 方便导出拟合曲线，同时删除中间异常部分
% XFittedLogG_part = XFittedLogG(XFittedLogG<-0.1 | XFittedLogG>0.1)';
% YFittedLogG_part = YFittedLogG(XFittedLogG<-0.1 | XFittedLogG>0.1)';








