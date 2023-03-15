function [f,cov] = benchmark_func(x, func_num)
[objF,covG,covH]=CEC2017(x,func_num);
deta=0.0001;
%% 计算与最优值的误差
[np,~]=size(objF);
for i=1:np
    f(i,1)=objF(i,1);
    cov(i,1)=sum(max(0,covG(i,:)))+sum(max((abs(covH(i,:))-deta),0));
end
end