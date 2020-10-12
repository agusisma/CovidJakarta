%% Research code by Agus Hasan

% Disclaimer: the analysis and results are strictly only for educational
% and research purposes and may contain errors.

% The effective (real-time) reproduction number (Rt) is calculated based on
% extended Kalman filter. I added a low-pass filter to reduce short term
% data fluctuation.

% The transmission index (TI) is calculated as Rt/Rmax

%==========================================================================
% COVIDMETER version 2.0 (3 day simulation average)
%==========================================================================

clear;
clc;
tic
%% load data
load DKI.txt; % load data: date | month | susceptible | active cases | recovered cases | deceased cases | daily confirmed cases

DATA = DKI;

%% Inputs
Tinf     = 12;                                % infectious time
Std_Tinf = 3;                                % standard deviation of infectious time
N   = sum(DATA(1,3:end));                    % number of population
CFR = DATA(end,6)/(sum(DATA(end,4:6)))/2;    % case fatality rate (with correction factor of 4)
tf  = length(DATA);                          % simulation time
tp  = 60;                                    % prediction time
td  = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf);
tdp = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf+tp);
dt  = 0.01;                                  % time steps
t   = dt:dt:tf;

%% Data matrix
C = [1 0 0 0 0 0;     % We have data of S, I, R, D, and C.
     0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 1 0 0
     0 0 0 0 1 0];
%% Parameters
sigma  = 1.96; %95 Confident Interval for infectious time

%% Noise
QF = diag([10 10 10 10 5 0.2]);   % process and measurement covariance matrices
RF = diag([100 10 10 5 1]);       % are considered as tuning parameters

%% For plotting
% Adding a low pass-filter to handle short-term data fluctuations 
windowSize = 300;
b          = (1/windowSize)*ones(1,windowSize);
a          = 1;
% Plotting Rt below 1
curve11 = 0*ones(1,tf);
curve22 = 1*ones(1,tf);
x2      = [td, fliplr(td)];
fs      = 48;
fs1     = 24;

%% Simulation
for j = 1:3
% Infection time
Ti     = Tinf-Std_Tinf*sigma+(j-1)*Std_Tinf*sigma;    % infection time with standard dev. 1 day
gamma  = (1-CFR)*(1/Ti);            % recovery rate
kappa  = CFR*1/Ti;                  % death rate

%% Initialization
xhat     = [N-1; 1; 0; 0; 1; 0];   % initial condition
Pplus    = 1000*eye(6);            % since we know excatly the initial conditions
xhatEff  = 0;
% for plotting
xArray       = [];
xhatArray    = [];
xhatEffArray = [];
% extended Kalman filter
for i=1:((tf-1)/dt)
     xhatArray    = [xhatArray xhat]; 
     xhatEffArray = [xhatEffArray xhatEff];      
     % assimilating the reported data
     y = [interp1(0:1:tf-1,DATA(:,3),t);
         interp1(0:1:tf-1,DATA(:,4),t);
         interp1(0:1:tf-1,DATA(:,5),t);
         interp1(0:1:tf-1,DATA(:,6),t);
         interp1(0:1:tf-1,DATA(:,7),t)];
     %y = y + sqrt(RF)*[randn randn randn randn]'*dt;
     % prediction
     xhat(1) = xhat(1)-(gamma+kappa)*xhat(6)*xhat(1)*xhat(2)*dt/N;
     xhat(2) = xhat(2)+(gamma+kappa)*xhat(6)*xhat(1)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt;
     xhat(3) = xhat(3)+gamma*xhat(2)*dt;
     xhat(4) = xhat(4)+kappa*xhat(2)*dt;
     xhat(5) = xhat(5)+(gamma+kappa)*xhat(2)*dt-xhat(5)*dt;
     xhat(6) = xhat(6);
    % calculating the Jacobian matrix
    FX    = [1-(gamma+kappa)*xhat(6)*xhat(2)*dt/N -(gamma+kappa)*xhat(6)*xhat(1)*dt/N 0 0 0 -(gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             (gamma+kappa)*xhat(6)*xhat(2)*dt/N 1+(gamma+kappa)*xhat(6)*xhat(1)*dt/N-(gamma+kappa)*dt 0 0 0 (gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             0 gamma*dt 1 0 0 0;
             0 kappa*dt 0 1 0 0;
             0 (gamma+kappa)*dt 0 0 1-dt 0;
             0 0 0 0 0 1];
    Pmin  = FX*Pplus*FX'+QF;
    % update 
    KF    = Pmin*C'*inv(C*Pmin*C'+RF);  % Kalman gain
    xhat  = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(6)-KF*C)*Pmin;
    xhat(6) = max(0,xhat(6));           % the reproduction number cannot be negative
    xhatEff = (xhat(1)/N)*xhat(6);      % calculating the effective repsoduction number
end

%% Plotting

xhatArray(6,:) = filter(b,a,xhatEffArray);

xhatSArray  = [];
xhatS       = xhatArray(1,tf);
xhatIArray  = [];
xhatI       = xhatArray(2,tf);
xhatRArray  = [];
xhatR       = xhatArray(3,tf);
xhatDArray  = [];
xhatD       = xhatArray(4,tf);
xhatHArray  = [];
xhatH       = xhatArray(5,tf);
xhatRtArray = [];
xhatRt      = xhatArray(6,tf);
for i=1:tf-1
    xhatSArray  = [xhatSArray xhatS];
    xhatS       = xhatArray(1,(1/dt)*i);
    xhatIArray  = [xhatIArray xhatI];
    xhatI       = xhatArray(2,(1/dt)*i);
    xhatRArray  = [xhatRArray xhatR];
    xhatR       = xhatArray(3,(1/dt)*i);
    xhatDArray  = [xhatDArray xhatD];
    xhatD       = xhatArray(4,(1/dt)*i);
    xhatHArray  = [xhatHArray xhatH];
    xhatH       = xhatArray(5,(1/dt)*i);
    xhatRtArray = [xhatRtArray xhatRt];
    xhatRt      = xhatArray(6,(1/dt)*i);
end

xhatSArray  = [xhatSArray xhatS];
xhatIArray  = [xhatIArray xhatI];
xhatRArray  = [xhatRArray xhatR];
xhatDArray  = [xhatDArray xhatD];
xhatHArray  = [xhatHArray xhatH];
xhatRtArray = [xhatRtArray xhatRt];

M(j,:) = smooth(xhatRtArray);

end

for l = 1:tf
    curve2(l)      = max(M(:,l));
    xhatRtArray(l) = mean(M(:,l));
    curve1(l)      = min(M(:,l));
end

    %% Forecasting

    mRt  = mean(xhatRtArray(end-14:end));
    R0   = max(xhatRtArray);
    XRT = [1-(xhatRt/R0) xhatRt(end)/R0];
    
    CI = (mean(xhatRtArray(end-14:end))/R0);
    if XRT(end) > 1
        XRT = [0 1];
        CI = 1;
    end
    cc  = [.1 .5 1.];
    [val,id1] = min(abs(CI-cc));
    cc(id1) = CI;
    
   
    for m = 1:3

        contact = cc(m);
        xp   = [DATA(end,3); DATA(end,4); DATA(end,5); DATA(end,6); DATA(end,7)]; % initial condition
        xpArray     = [];
        
        cRt = contact*R0;
        RtBArray = [];
        for j=1:tp/dt
            RtB = mRt-((mRt-cRt)/(tp/dt))*j;
            RtBArray = [RtBArray RtB];
        end
        Rt = RtBArray;
        for i=1:tp/dt
            xpArray = [xpArray xp];
            
            xp(1) = xp(1)-(gamma+kappa)*Rt(i)*xp(1)*xp(2)*dt/N;
            xp(2) = xp(2)+(gamma+kappa)*Rt(i)*xp(1)*xp(2)*dt/N-(gamma+kappa)*xp(2)*dt;
            xp(3) = xp(3)+gamma*xp(2)*dt;
            xp(4) = xp(4)+kappa*xp(2)*dt;
            xp(5) = xp(5)+(gamma+kappa)*xp(2)*dt-xp(5)*dt;
        end
        
        xpSArray  = [];
        xpS       = xpArray(1,tf);
        xpIArray  = [];
        xpI       = xpArray(2,tf);
        xpRArray  = [];
        xpR       = xpArray(3,tf);
        xpDArray  = [];
        xpD       = xpArray(4,tf);
        xpHArray  = [];
        xpH       = xpArray(5,tf);
        for i=1:tp
            xpSArray  = [xpSArray xpS];
            xpS       = xpArray(1,(1/dt)*i);
            xpIArray  = [xpIArray xpI];
            xpI       = xpArray(2,(1/dt)*i);
            xpRArray  = [xpRArray xpR];
            xpR       = xpArray(3,(1/dt)*i);
            xpDArray  = [xpDArray xpD];
            xpD       = xpArray(4,(1/dt)*i);
            xpHArray  = [xpHArray xpH];
            xpH       = xpArray(5,(1/dt)*i);
        end

        xIpredic(m,:) = [xhatIArray xpIArray];
        xRpredic(m,:) = [xhatRArray xpRArray];
        xDpredic(m,:) = [xhatDArray xpDArray];
        xHpredic(m,:) = [xhatHArray xpHArray];        
    end

toc

figure(1)
subplot(2,2,1)
plot(tdp,xIpredic(3,:),':y','LineWidth',6)
hold on;
plot(tdp,xIpredic(2,:),':b','LineWidth',6)
hold on;
plot(tdp,xIpredic(1,:),':g','LineWidth',6)
hold on;
plot(td,xhatIArray,':r','LineWidth',6)
set(gca,'color','none','FontSize',fs1)
plot([datetime(2020,DATA(end,2),DATA(end,1)) datetime(2020,DATA(end,2),DATA(end,1))],[0 1e+9],'--k','LineWidth',4)
%ylim([0 max(xIpredic(3,:))])
ylim([0 10000])
xlim([min(tdp) max(tdp)])
title('Active Case')
grid on
grid minor
xtickangle(45)
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs1-4)
subplot(2,2,2)
plot(tdp,xRpredic(3,:),':y','LineWidth',6)
hold on;
plot(tdp,xRpredic(2,:),':b','LineWidth',6)
hold on;
plot(tdp,xRpredic(1,:),':g','LineWidth',6)
hold on;
plot(td,xhatRArray,':r','LineWidth',6)
set(gca,'color','none','FontSize',fs1)
plot([datetime(2020,DATA(end,2),DATA(end,1)) datetime(2020,DATA(end,2),DATA(end,1))],[0 1e+9],'--k','LineWidth',4)
ylim([0 max(xRpredic(3,:))])
xlim([min(tdp) max(tdp)])
title('Recovered Case')
grid on
grid minor
xtickangle(45)
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs1-4)
subplot(2,2,3)
plot(tdp,xDpredic(3,:),':y','LineWidth',6)
hold on;
plot(tdp,xDpredic(2,:),':b','LineWidth',6)
hold on;
plot(tdp,xDpredic(1,:),':g','LineWidth',6)
hold on;
plot(td,xhatDArray,':r','LineWidth',6)
set(gca,'color','none','FontSize',fs1)
plot([datetime(2020,DATA(end,2),DATA(end,1)) datetime(2020,DATA(end,2),DATA(end,1))],[0 1e+9],'--k','LineWidth',4)
ylim([0 max(xDpredic(3,:))])
xlim([min(tdp) max(tdp)])
title('Deceased Case')
grid on
grid minor
xtickangle(45)
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs1-4)
subplot(2,2,4)
plot(tdp,xIpredic(3,:)+xRpredic(3,:)+xDpredic(3,:),':y','LineWidth',6)
hold on;
plot(tdp,xIpredic(2,:)+xRpredic(2,:)+xDpredic(2,:),':b','LineWidth',6)
hold on;
plot(tdp,xIpredic(1,:)+xRpredic(1,:)+xDpredic(1,:),':g','LineWidth',6)
hold on;
plot(td,xhatIArray+xhatRArray+xhatDArray,':r','LineWidth',6)
set(gca,'color','none','FontSize',fs1)
plot([datetime(2020,DATA(end,2),DATA(end,1)) datetime(2020,DATA(end,2),DATA(end,1))],[0 1e+9],'--k','LineWidth',4)
ylim([0 max(xIpredic(3,:)+xRpredic(3,:)+xDpredic(3,:))])
xlim([min(tdp) max(tdp)])
title('Total Case')
grid on
grid minor
xtickangle(45)
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs1-4)