function neibor= BuildNeibor(neibor_size,reconnect_p,min_degree,pop_size,sorted_index)
%% 构建最近邻网络的邻域
neibor=zeros(pop_size,pop_size);
for i=neibor_size/2+1:pop_size-neibor_size/2
    neibor(sorted_index(i),sorted_index(i+1:i+neibor_size/2))=1;
    neibor(sorted_index(i),sorted_index(i-1:-1:i-neibor_size/2))=1;
end
for i=1:neibor_size/2
    neibor(sorted_index(i),sorted_index(i+1:i+neibor_size/2))=1;
    for j=i-1:-1:i-neibor_size/2
        if j<1
            neibor(sorted_index(i),sorted_index(pop_size+j))=1;
        else
            neibor(sorted_index(i),sorted_index(j))=1;
        end
    end
end
for i=pop_size-neibor_size/2+1:pop_size
    neibor(sorted_index(i),sorted_index(i-1:-1:i-neibor_size/2))=1;
    for j=i+1:i+neibor_size/2
        if j>pop_size
            neibor(sorted_index(i),sorted_index(j-pop_size))=1;
        else
            neibor(sorted_index(i),sorted_index(j))=1;
        end
    end
end
%% 边随机重连，构建小世界网络邻域
for i=1:pop_size
    lock=[i];
    for j=1:pop_size
        if (neibor(i,j)==1)&&(ismember(j,lock)==0)&&(rand<reconnect_p)
            index0=find(neibor(i,:)==0);
            index0(find(index0==i))=[];
            in=randi(length(index0));
            lock(end+1)=index0(in);
            neibor(i,j)=0;neibor(j,i)=0;
            neibor(i,index0(in))=1;neibor(index0(in),i)=1;
        end
    end
end
%% 补度
for i=1:pop_size
    degree=length(find(neibor(i,:)==1));
    if degree<min_degree
        index0=find(neibor(i,:)==0);
        in=randperm(length(index0));
        neibor(i,index0(in(1:min_degree-degree)))=1;neibor(index0(in(1:min_degree-degree)),i)=1;
    end
end