function lu=GetLU(func,problem_size)
Xmax = 100 + zeros(28,1);
Xmax(4) = 10;
Xmax(5) = 10;
Xmax(6) = 20;
Xmax(7) = 50;
Xmax(9) = 10;
Xmax(19) = 50;
Xmax(28) = 50;
 lu = [-Xmax(func) * ones(1, problem_size); Xmax(func) * ones(1, problem_size)];
end