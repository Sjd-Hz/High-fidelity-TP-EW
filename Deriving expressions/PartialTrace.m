function [result] = PartialTrace(n1, n2)
%PartialTrace: this is a function of the partial trace operation on the quantum subsystem
%Input: n1 is the density matrix of the whole quantum system
%       n2 is a array which contains the sequence numbers that make up the quantum subsystem
%Output: result is the density matrix after the partial trace operation
%Example: n1 is a 8*8 density matrix, n2=[1,3], this respresents that the
%partial trace on the quantum subsystem consisting of the first and third
%quantum bit
result=0;
temp=1;
I=eye(2);
s0=[1;0];
s1=[0;1];
base={};
%obtain the basis corresponding to the subsystem (n2)
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
%partial trace operation
for i=1:size(base,2)
    result=result+base{i}'*n1*base{i};
end

end

