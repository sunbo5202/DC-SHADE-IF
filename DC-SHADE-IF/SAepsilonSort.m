function [temp_fit, sorted_index]=SAepsilonSort(fit,cov,direct,epsilon,t)
%% 采用SAeplision排序，第一次无epsilon，采用SF
if t==1
    markc=find(cov==0);
    markf=find(cov~=0);
else
    markc=find(cov<=epsilon);
    markf=find(cov>epsilon);
end
%% 可行解与不可行解分成两组排序
[~,index_markc]=sort(fit(markc));
orderc=markc(index_markc);
[~,index_markf]=sort(cov(markf));
orderf=markf(index_markf);
%% 不可行解在约束违反相同的情况下适应度好的优先
i=1;
while i<length(orderf)
    index_orderf=find(cov(orderf)==cov(orderf(i)));
    %index_orderf=orderf(index_orderf);
    length_io=length(index_orderf);
    if length_io>1
        [~,index_of]=sort(fit(orderf(index_orderf)));
        orderf(i:i+length_io-1)=orderf(index_orderf(index_of));
        i=i+length_io;
    else
        i=i+1;
    end
end
sorted_index=[orderc;orderf];
temp_fit=fit(sorted_index);
switch direct
    case 'ascend'
        
    case 'descend'
        sorted_index=flipud(sorted_index);
        temp_fit=flipud(temp_fit);
end