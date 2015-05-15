function spW=ComputeDecomposition(W)
n=size(W,1);
W(abs(W)>0)=-1;
W(eye(n)>0)=-sum(W,2);
W(eye(n)>0)=max(W(eye(n)>0),1);
spW=chordalDecomposition(W);
end