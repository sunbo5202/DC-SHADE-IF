%% 较上一版本，约束处理改为两阶段，第一阶段搜索不可行域
%% 最终版本，下一步调参，待调参数：
%% switch_rate：0.1-0.5 阶段转换最高代数控制
%% DI：0.1-0.9 小世界网络重连概率扰动强度控制
%% neibor_rate：0.1-0.9 小世界网络邻域大小控制 %%确定为0.3，下一步调整DI
%% DI确定为0.2，下一步调整switch_rate
clc;
clear;
format long e;
format compact;
global initial_flag
problem_size = 10;
n = problem_size;
max_nfes = 20000 * problem_size;
rand('seed', sum(100 * clock));
for func = 1:28
    t1=cputime;
    %% 用于问题框架的参数
    lu=GetLU(func,problem_size);
    problem = func;
    %% Record the best results
    outcome = [];
    fprintf('\n-------------------------------------------------------\n')
    fprintf('Function = %d, Dimension size = %d\n', func, problem_size)
    %% ========================================
    for run_id = 1 : 25
        %%  parameter settings for SHADE
        p_best_rate = 0.1;
        arc_rate = 2;
        memory_size = problem_size;
        pop_size = 100;
        %% =============================
        %% 用于分阶段的参数
        phase=1;switch_rate=0.4;
        %% =================
        %% 用于SAepsilon的参数
        N=pop_size;
        Tmax=max_nfes/pop_size;
        Tc=0.2*Tmax;
        theta=0.3*N;
        cp=5;
        Th1=100;
        Th2=2;
        %% ============================
        %% Initialize the main population
        popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
        pop = popold; % the old population becomes the current population
        %fitness = benchmark_func(pop, problem);
        initial_flag = 0;
        [fit_p,cov_p] = benchmark_func(pop, problem);
        if phase==1
            fitness=cov_p;
        else
            fitness=fit_p;
        end
        t=1;
        %% =============
        nfes = 0;
        bsf_fit_var = 1e+30;%%要改
        bsf_solution = zeros(1, problem_size);
        
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        for i = 1 : pop_size
            nfes = nfes + 1;
            
            if fitness(i) < bsf_fit_var
                bsf_fit_var = fitness(i);
                bsf_solution = pop(i, :);
            end
            if nfes > max_nfes; break; end
        end
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        
        memory_sf = 0.5 .* ones(memory_size, 1);
        memory_cr = 0.5 .* ones(memory_size, 1);
        memory_pos = 1;
        
        archive.NP = arc_rate * pop_size; % the maximum size of the archive
        archive.pop = zeros(0, problem_size); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions
        G=1;
       
        %% main loop
        while nfes < max_nfes
            %% 
            pop = popold; % the old population becomes the current population
            if t==1
                [temp_fit, sorted_index] =SAepsilonSort(fitness,cov_p,'ascend',0,t);
            else
                [temp_fit, sorted_index] =SAepsilonSort(fitness,cov_p,'ascend',epsilon(end),t);
            end        
            mem_rand_index = ceil(memory_size * rand(pop_size, 1));
            mu_sf = memory_sf(mem_rand_index);
            mu_cr = memory_cr(mem_rand_index);
            
            %% for generating crossover rate
            cr = normrnd(mu_cr, 0.1);
            term_pos = find(mu_cr == -1);
            cr(term_pos) = 0;
            cr = min(cr, 1);
            cr = max(cr, 0);
            
            %% for generating scaling factor
            sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
            pos = find(sf <= 0);
            while ~ isempty(pos)
                sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
                pos = find(sf <= 0);
            end
            sf = min(sf, 1);
            pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
            randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions

            %% 根据FDC扰动
            DI=0.2;%扰动强度disturbance intensity
            p(G)=GetReconnectPro(pop,sorted_index,pop_size,DI,fitness,nfes,max_nfes);
            %% 构建小世界网络邻域
            reconnect_p=p(G);%重连概率
            neibor_rate=0.3;
            neibor_size=neibor_rate*pop_size;%邻域大小
            min_degree=2;%最小度
            neibor= BuildNeibor(neibor_size,reconnect_p,min_degree,pop_size,sorted_index);
            
            %% r1,r2从邻域中选择
            [r1,r2]=GetR1R2(neibor,pop_size);
            %% 变异
            vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - pop(r2, :));
            vi = boundConstraint(vi, pop, lu);
            
            mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop(mask);
            
            [fit_c,cov_c] = benchmark_func(ui, problem);
            if phase==1
                children_fitness=cov_c;
            else
                children_fitness=fit_c;
            end
            %% 更新epsilon
            if phase==2
                if t==1
                    phi_order_valueP=sort(cov_p,'descend');
                    epsilon0=phi_order_valueP(ceil(theta));
                    epsilon1(t)=epsilon0;
                else
                    if epsilon0>Th1
                        phi_order_valueC=sort(cov_c,'descend');
                        epsilon2(t)=phi_order_valueC(ceil(theta));
                    else
                        epsilon2(t)=epsilon0;
                    end
                    if (epsilon2(t)>Th2) && (epsilon2(t)<epsilon1(t-1))
                        epsilon1(t)=epsilon2(t);
                    else
                        epsilon1(t)=epsilon1(t-1);
                    end
                end
                if t<=Tc
                    epsilon(t)=epsilon1(t)*((1-(t/Tc))^cp);
                else
                    epsilon(t)=0;
                end
                t=t+1;
            end
            %% =================================
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            for i = 1 : pop_size
                nfes = nfes + 1;
                
                if children_fitness(i) < bsf_fit_var
                    bsf_fit_var = children_fitness(i);
                    bsf_solution = ui(i, :);
                end
                if nfes > max_nfes; break; end
            end
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            if phase==2
                dif = abs((fitness+cov_p) - (children_fitness+cov_c));
            else
            	dif=abs(fitness-children_fitness);
            end
            
            %% I == 0: the parent is better; I == 1: the offspring is better
            if phase==1
            	I = (fitness > children_fitness);
            else
                for i=1:pop_size
                    if ((cov_p(i,1)<epsilon(end))&&(cov_c(i,1)<epsilon(end))) || (cov_p(i,1)==cov_c(i,1))
                        if children_fitness(i,1)<fitness(i,1)
                            I(i,1)=1;
                        else
                            I(i,1)=0;
                        end
                    else
                        if cov_c(i,1)<=cov_p(i,1)
                            I(i,1)=1;
                        else
                            I(i,1)=0;
                        end
                    end
                end
            end
            goodCR = cr(I == 1);
            goodF = sf(I == 1);
            dif_val = dif(I == 1);
            archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
            fitness(I==1)=children_fitness(I==1);
            cov_p(I==1)=cov_c(I==1);
            
            popold = pop;
            popold(I == 1, :) = ui(I == 1, :);
            
            num_success_params = numel(goodCR);
            
            if num_success_params > 0
                sum_dif = sum(dif_val);
                dif_val = dif_val / sum_dif;
                
                %% for updating the memory of scaling factor
                memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
                
                %% for updating the memory of crossover rate
                if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
                    memory_cr(memory_pos)  = -1;
                else
                    memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
                end
                
                memory_pos = memory_pos + 1;
                if memory_pos > memory_size;  memory_pos = 1; end
            end
            %% 阶段转换
            if phase==1
                phase1_iter_result(G)=bsf_fit_var;
                if (length(find(fitness==0))/pop_size>=0.5) || G>((max_nfes/pop_size)*switch_rate)
                    phase=2;
                    fitness=fit_p;
                    switchG=G;
                end
            end
            G=G+1;
        end
        [~,index_fin]=SAepsilonSort(fitness,cov_p,'ascend',0,1);
        bsf_fit_var=fitness(index_fin(1));bsf_cov_var=cov_p(index_fin(1));
        fprintf('%d th run, bsf_fit_var = %1.8e, bsf_cov_var = %1.8e\n', run_id , bsf_fit_var,bsf_cov_var)
        SRF04(func,run_id)=bsf_fit_var;
        SRC04(func,run_id)=bsf_cov_var;
        shadeSwitchG(func,run_id)=switchG;
        outcome = [outcome bsf_fit_var];
    end %% end 1 run
    %%  sort(outcome)
    
    fprintf('\n')
    fprintf('mean error value = %f, std = %f\n', mean(outcome), std(outcome))
    fprintf('running time= %f\n', cputime-t1)
end %% end 1 function run
save('ParametersSR0_4.mat','SRF04','SRC04')