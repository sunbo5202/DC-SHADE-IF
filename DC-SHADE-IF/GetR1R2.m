function [r1,r2]=GetR1R2(neibor,pop_size)
r1=ones(1,pop_size);r2=ones(1,pop_size);
for i=1:pop_size
    index1=find(neibor(i,:)==1);
    order1=randperm(length(index1));
    r1(i)=index1(order1(1));r2(i)=index1(order1(2));
end
end