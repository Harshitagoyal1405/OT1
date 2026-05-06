clc
clear all;
a = [1 2 2 ; 2 4 1 ];
b =[10;15];
c = [5 2 3];
n = size(a,2)
s = eye(size(a,1))
A = [a s b]
cost = zeros(1,size(A,2))
cost(1:n) = c

bv = n+1:1:size(A,2)-1
zjcj = cost(bv)*A-cost
zcj=[A;zjcj]
simptable=array2table(zcj)
simptable.Properties.VariableNames(1:size(zcj,2))={'x1','x2','x3','s1','s2','sol'}

run=true
while(run)
    zc=zjcj(1:end-1)
if any(zc<0)  
    fprintf('the given soln is not optimal')
    [min_val pvt_col] = min(zc)
    if all(A(:,pvt_col)<=0)
        eror('the given lpp unbounded')
    else
        sol= A(:,end)
        col= A(:,pvt_col)
    for i = 1:size(A,1)
        if col(i)>0
            ratio(i) = sol(i)/col(i)
        else
            ratio(i)= inf
        end
    end
    [minratio pvt_row] = min(ratio)
    end
    pvt_key=A(pvt_row,pvt_col)
    bV(pvt_row)= pvt_col
    A(pvt_row, :)= (A(pvt_row,:))/pvt_key
    for i=1:size(A,1)
        if i~=pvt_row
            A(i,:)= A(i,:)-A(i,pvt_col)* A(pvt_row,:)
        end
    end
    zjcj = zjcj-zjcj(pvt_col)*A(pvt_row,:)
    zcj=[A;zjcj]
    array2table(zcj, 'VariableNames',{'X1','X2','X3','S1','S2','SOL'},'RowNames',{'BV1', 'BV2', 'Zj-Cj'})
else
  fprintf("final optimal table has been obtained and optimal value is %f",zcj(end))
  run= false
end
end