S=load('mycielskian13.mat');
B=S.Problem.A;
tic
C=B.*(B*B);
e=ones(size(B,1),1);
c3=(C*e);
val=sum(c3,'all');
toc