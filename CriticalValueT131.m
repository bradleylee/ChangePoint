%%Critical Value Calculation for Table 1.3.1

%% exponetial 
N               = [20,50,100,500];
Confi_Level     = [0.90 0.95 0.99];
len_N           = length(N);
len_CL          = length(Confi_Level);
sim_times       = 10000;
%%
CriticalValue   = nan(len_N*len_CL,1);
for i = 1:len_N
    n       = N(i);
    LogL    = nan(sim_times,1);
    for times = 1:sim_times
        data    = exprnd(1,n,1);
        logl    = [];
        for k = 1:n-1
            xk      = mean(data(1:k));
            xks     = mean(data(k+1:end));
            xn      = mean(data);
            logl    = [logl,2*(k*log(1/xk)+(n-k)*log(1/xks)-n*log(1/xn))];
        end
        LogL(times) = sqrt(max(logl));
    end
    CriticalValue((i-1)*len_CL+1:i*len_CL) = quantile(LogL,Confi_Level);    
end

%% Possion
CriticalValue   = nan(len_N*len_CL,1);
for i = 1:len_N
    n       = N(i);
    LogL    = nan(sim_times,1);
    for times = 1:sim_times
        data    = poissrnd(1,n,1);
        logl    = [];
        for k = 1:n-1
            xk      = mean(data(1:k));
            xks     = mean(data(k+1:end));
            xn      = mean(data);
            logl    = [logl,2*(k*xk*log(xk)+(n-k)*xks*log(xks)-n*xn*log(xn))];
        end
        LogL(times) = sqrt(max(logl));
    end
    CriticalValue((i-1)*len_CL+1:i*len_CL) = quantile(LogL,Confi_Level);    
end
%% Normal with known sigma^2
CriticalValue   = nan(len_N*len_CL,1);
for i = 1:len_N
    n       = N(i);
    LogL    = nan(sim_times,1);
    for times = 1:sim_times
        data    = normrnd(0,1,n,1);
        sigma   = 1;
        logl    = [];
        for k = 1:n-1
            sk      = sum(data(1:k));
            sn      = sum(data);
            logl    = [logl,(n/sigma/k/(n-k))*(sk-k*sn/n)^2];
        end
        LogL(times) = sqrt(max(logl));
    end
    CriticalValue((i-1)*len_CL+1:i*len_CL) = quantile(LogL,Confi_Level);    
end

%% Normal with unknown sigma^2

CriticalValue   = nan(len_N*len_CL,1);
for i = 1:len_N
    n       = N(i);
    LogL    = nan(sim_times,1);
    for times = 1:sim_times
        data    = normrnd(0,rand(1),n,1);
        logl    = [];
        for k = 1:n-1
            sigmak  = (var(data(1:k))+var(data(k+1:end)))/n;
            sk      = sum(data(1:k));       
            sn      = sum(data);
            sigman  = sigmak + (n/k/(n-k))*(sk-k*sn/n)^2;
            logl    = [logl,n*(log(sigman)-log(sigmak))];
        end
        LogL(times) = sqrt(max(logl));
    end
    CriticalValue((i-1)*len_CL+1:i*len_CL) = quantile(LogL,Confi_Level);    
end

%% Normal with change sigma^2

CriticalValue   = nan(len_N*len_CL,1);
for i = 1:len_N
    n       = N(i);
    LogL    = nan(sim_times,1);
    for times = 1:sim_times
        data    = normrnd(0,rand(1),n,1);
        logl    = [];
        for k = 1:n-1
            sigmak  = (var(data(1:k))+var(data(k+1:end)))/n;
            sk      = sum(data(1:k));       
            sn      = sum(data);
            sigman  = sigmak + (n/k/(n-k))*(sk-k*sn/n)^2;
            sigmaks = var(data(k+1:end))/(n-k);
            logl    = [logl,n*log(sigman)-k*log(sigmak)-(n-k)*log(sigmaks)];
        end
        LogL(times) = sqrt(max(logl));
    end
    CriticalValue((i-1)*len_CL+1:i*len_CL) = quantile(LogL,Confi_Level);    
end