function p=GetReconnectPro(pop,sorted_index,pop_size,DI,fitness,nfes,max_nfes)
%% 计算FDC
for i=1:pop_size
    distance(i,1)=sqrt(sum((pop(i,:)-pop(sorted_index(1),:)).^2));
end
cfd=(sum((fitness-mean(fitness)).*(distance-mean(fitness))))/pop_size;
sfd=std(fitness)*std(distance);
if cfd==sfd
    fdc=1;
elseif cfd==0
    fdc=0;
else
    fdc=cfd/sfd;
end
x=2*(1-(nfes/max_nfes));
p=(exp(x)-exp(-x))/(exp(x)+exp(-x));

if fdc<=-0.15
    p=p*DI;%减小
    
elseif fdc>0.15
    p=p*(1+DI);%增大
end