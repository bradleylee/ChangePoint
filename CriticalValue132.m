%%Critical Value Calculation by Theroem 1.3.2
B   = @(h,l)log((1-h)*(1-l)/h/l);
D   = @(d,x)(x^d)*exp(-x^2/2)/2^(d/2)/gamma(d/2);
N   = [20,50,100,500,1000];
Confi_Level = [0.90 0.95 0.99];
len_N = length(N);
len_CL = length(Confi_Level);
CriticalValue = nan(len_N*len_CL,1);
for i = 1:len_N
    n = N(i);
    h = log(n)^1.5/n;
    l = h;
    for j = 1:len_CL
        syms x
        x = solve(D(1,x)*(B(l,h)-B(l,h)/(x^2)+4/(x^2))== 1-Confi_Level(j));
        CriticalValue(len_CL*(i-1)+j)   = eval(x); %Table 1.3.1 column 4 u_hat
    end
end
