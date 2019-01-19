%% Programs for Limit Theorems in Change-Point Analysis by Scorgo and Horvath 1997

%Critical Value Calculation by Theroem 1.3.1 


A   = @(x)sqrt(2*log(x));
D   = @(d,x)2*log(x)+0.5*d*log(log(x))-log(gamma(d/2));
N   = [20,50,100,500,1000]; % Table 1.3.1 column 1
Confi_Level     = [0.90 0.95 0.99]; % Table 1.3.1 column 2: 1-alpha
len_N           = length(N);
len_CL          = length(Confi_Level);
CriticalValue   = nan(len_N*len_CL,1);
for i = 1:len_N
    n = N(i);
    for j = 1:len_CL
        t                               = -log(log(Confi_Level(j))/(-2));
        CriticalValue(len_CL*(i-1)+j)   = (t+D(1,log(n)))/A(log(n));    %Table 1.3.1 column 3
    end
end

