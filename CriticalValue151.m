%%Critical Value Calculation for Table 1.5.1
Confi_Level     = 0.95;%[0.90 0.95 0.99];
Normal_CV = nan(size(Confi_Level));
N       = 100;
KS      = 10:10:50;
theta1  = 0;
theta2  = 1;
sigma   = 1;
% for i = 1:length(Confi_Level)   
%     syms x
%     X = solve(normcdf(x,0,1)==Confi_Level(i));
%     Normal_CV(i) = abs(eval(X));
% end
Normal_CV =3.06;
%% Normal with constant known variance
sigma1  = @(lambda,theta1,theta2)lambda*(1-lambda)*(theta1-theta2)^2;
sigma2  = @(lambda,sigma)lambda*(1-lambda)*sigma^2;
sigma3  = @(theta1,theta2,sigma)(theta1-theta2)^2/sigma^2;
sigma4  = @(sigma)sigma^2;
mus     = @(theta1,theta2,sigma,ks,n)(1/2/sigma^2)*(ks*(n-ks)/n)*(theta1-theta2)^2;

Tt      = @(t,sigma)t/sigma^2; 
Ht      = @(t,sigma)t^2*sigma^2/2;
dA      = @(theta,sigma)theta/sigma^2;
tao1    = dA(theta1,sigma);
tao2    = dA(theta2,sigma);
delta   = (tao1-tao2);

%% 
Times   = 5000;
Data    = nan(N,1);
for K = 1%:length(KS)
    kstar  = KS(K);
    lambda = kstar/N;
    %% asympotic value
    Zn1 = Normal_CV*sqrt(4*N*sigma1(lambda,theta1,theta2))+2*mus(theta1,theta2,sigma,kstar,N);
    Zn2 = Normal_CV*sqrt(4*N*sigma2(lambda,sigma))+2*mus(theta1,theta2,sigma,kstar,N);
    Zn3 = Normal_CV*sqrt(4*kstar*sigma3(theta1,theta2,sigma))+2*mus(theta1,theta2,sigma,kstar,N);
    Zn4 = Normal_CV*sqrt(4*kstar*delta^2*sigma4(sigma))+2*mus(theta1,theta2,sigma,kstar,N);
    %% simulation
    LogL    = nan(Times,1);
    for T = 1:Times
        Data(1:kstar)       = normrnd(theta1,sigma,kstar,1);
        Data(kstar+1:end)   = normrnd(theta2,sigma,N-kstar,1);
%         Data                = normrnd(theta1,sigma,N,1);
%         LogL(T)             = 2*(kstar*Ht(mean(Tt(Data(1:kstar),sigma)),sigma)+(N-kstar)*Ht(mean(Tt(Data(kstar+1:end),sigma)),sigma)-N*Ht(mean(Tt(Data,sigma)),sigma));
        logl = [];
        for i = 1:N-1            
            logl = [logl,2*(i*Ht(mean(Tt(Data(1:i),sigma)),sigma)+(N-i)*Ht(mean(Tt(Data(i+1:end),sigma)),sigma)-N*Ht(mean(Tt(Data,sigma)),sigma))];
        end
        LogL(T) = max(logl);
    end   
    Critical_Value_Sim = quantile(LogL,Confi_Level);
    Power = nan(Times,5);
    for T = 1:Times
        Data(1:kstar)       = normrnd(theta1,sigma,kstar,1);
        Data(kstar+1:end)   = normrnd(theta2,sigma,N-kstar,1);
        logl = [];
        for i = 1:N-1            
            logl = [logl,2*(i*Ht(mean(Tt(Data(1:i),sigma)),sigma)+(N-i)*Ht(mean(Tt(Data(i+1:end),sigma)),sigma)-N*Ht(mean(Tt(Data,sigma)),sigma))];
        end
        Statistic  = max(logl);
%         Statistic  = 2*(kstar*Ht(mean(Tt(Data(1:kstar),sigma)),sigma)+(N-kstar)*Ht(mean(Tt(Data(kstar+1:end),sigma)),sigma)-N*Ht(mean(Tt(Data,sigma)),sigma));
        Power(T,1) = Statistic>Zn1;
        Power(T,2) = Statistic>Zn2;
        Power(T,3) = Statistic>Zn3;
        Power(T,4) = Statistic>Zn4;
        Power(T,5) = Statistic>Critical_Value_Sim;
    end
    sum(Power)/Times
    
    
end