function [result] = PartialTrace(n1, n2)
%PartialTrace4 此处显示有关此函数的摘要
%   此处显示详细说明
result=0;
temp=1;
I=eye(2);
s0=[1;0];
s1=[0;1];
base={};
%子系统n2的基底获取
for i=1:log2(size(n1,1))
    if ~ismember(i,n2)
        if size(base,2)==0
            base{1}=kron(temp,I);
        else
            for j=1:size(base,2)
                base{j}=kron(base{j},I);
            end
        end
    else
        if size(base,2)==0
            base{1}=kron(temp,s0);
            base{2}=kron(temp,s1);
        else
            nowsize=size(base,2);
           for j=1:nowsize
               base=[base,kron(base{j},s0)];
               base=[base,kron(base{j},s1)];
           end
           base(1:nowsize)=[];          
%            for j=1:nowsize
%                base(j)=[];  
%            end  
        end    
    end    
end
%偏迹运算
for i=1:size(base,2)
    result=result+base{i}'*n1*base{i};
end

end

