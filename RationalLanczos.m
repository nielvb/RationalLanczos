function [V,W,T,S,Stil,Ttil] = RatLan_paper(A,v,w,n,B,L,Lambda,Beta)
% RATLAN computes the biorthogonal bases of rational Krylov subspaces with 
%  given poles and the oblique projections on these subspaces as a pencil.
%
% This function constructs biorthogonal bases V and W (W*V=I) for two
% rational Krylov subspaces Krat and Lrat. M* denotes the conjugate 
% transpose of some matrix M.
% The subspaces are Krat(A,v,L,B) = span{v,(L(1)A- B(1)I)^{-1}Av,(L(2)A- B(2)I)^{-1}Av}
% and Lrat(A*,w,Beta,Lambda) = span{w,(Beta(1)A*- Lambda(1)I)^{-1}A*v,(Beta(2)A- Lambda(2)I)^{-1}Av}
% In other words B./L are poles of Krat and Lambda./beta of Lrat.
%
% The oblique projection of A onto Krat and orthogonal to Lrat is represented by
% the tridiagonal matrix-pair (T,S) such that W*AVS=T.
% The oblique projection of A* onto Lrat and orthogonal to Krat is represented by
% the tridiagonal matrix-pair (Stil,Ttil) such that V*A*WTtil=Stil.

% INPUT:
%       - A = some (m x m)-matrix 
%       - n = maximum size of subpaces
%       - L, B = variables determining pole i of Krat by their ratio (B(i)/L(i))
%       - Beta, Lambda = variables determining pole i of Lrat by their ratio (Lambda(i)/Beta(i))
%       - v = starting vector for constructing the subspace Krat
%       - w = starting vector for constructing the subspace Lrat
% Output: with k= min(n+1,m)
%       - V = (m x k)-matrix representing the basis for Krat
%       - W = (m x k)-matrix representing the basis for Lrat
%       - T,S = two tridiagonal matrices of size (k x n) representing the
%       projection of A onto Krat orthogonal to Lrat, i.e., W*AVS=T.
%       - Stil,Ttil = two tridiagonal matrices of size (k x n) representing the
%       projection of A* onto Lrat orthogonal to Krat, i.e., V*A*WTtil=Stil.

I = eye(size(A));
vbar_1 = (L(1)*A-B(1)*I)\v;
wbar_1 = (Beta(1)*A'-Lambda(1)*I)\w;

temp1 = dot(w,vbar_1)/dot(w,A*vbar_1);
temp2 = dot(v,wbar_1)/dot(A*v,wbar_1);

vhat_2 = -temp1*A*vbar_1 + vbar_1;
what_2 = -temp2*A'*wbar_1 + wbar_1;



[sigma,gamma_1,r,d_1] = NormalizeInit(vhat_2,what_2)

l_2 = L(1)/r;
b_2 = B(1)/r
beta_2 = Beta(1)/sigma;
lambda_2 = Lambda(1)/sigma;
c_1 = d_1*temp1;
delta_1 = gamma_1*temp2;


v_2 = vhat_2*r*d_1;
w_2 = what_2*sigma*gamma_1;

V = [v,v_2];
W = [w,w_2];



T(1,1) = d_1;
T(2,1) = b_2;
S(1,1) = c_1;
S(2,1) = l_2;
Stil(1,1) = gamma_1;
Stil(2,1) = lambda_2;
Ttil(1,1) = delta_1;
Ttil(2,1) = beta_2;


for i = 2:n
    if(i==2)
        betarand = randn(1); % some random value for beta
        lambdarand = Lambda(n);%randn(1); % some random value for lambda
        vtilde = 1/betarand' * ((L(i)*A-B(i)*I)\((betarand'*A-lambdarand*I)*V(:,i-1)));
        lrand = randn(1);
        brand = B(n);%randn(1);
        wtilde = 1/lrand' * ((Beta(i)*A'-Lambda(i)*I)\((lrand'*A'-brand*I)*W(:,i-1)));
    else
        if(Beta(i-2)==0)
            vtilde = ((L(i)*A-B(i)*I)\((-Lambda(i-2)*I)*V(:,i-1)));
        else
            vtilde = 1/Beta(i-2)' * ((L(i)*A-B(i)*I)\((Beta(i-2)'*A-Lambda(i-2)*I)*V(:,i-1)));
        end
        if(L(i-2)==0)
            wtilde = ((Beta(i)*A'-Lambda(i)*I)\((-B(i-2)*I)*W(:,i-1)));
        else
            wtilde = 1/L(i-2)' * ((Beta(i)*A'-Lambda(i)*I)\((L(i-2)'*A'-B(i-2)*I)*W(:,i-1)));
        end
    end
    
        vbar = (L(i)*A-B(i)*I)\V(:,i);
        wbar = (Beta(i)*A'-Lambda(i)*I)\W(:,i);
        
        tempv1 = (dot(W(:,i-1),vtilde)*dot(W(:,i),A*vbar) - dot(W(:,i),vtilde)*dot(W(:,i-1),A*vbar))/...
            (dot(W(:,i-1),vbar)*dot(W(:,i),A*vbar)-dot(W(:,i),vbar)*dot(W(:,i-1),A*vbar));
        tempv2 = tempv1 * dot(W(:,i),vbar)/dot(W(:,i),A*vbar) - dot(W(:,i),vtilde)/dot(W(:,i),A*vbar);
        tempw1 = (dot(V(:,i-1),wtilde)*dot(A*V(:,i),wbar)-dot(V(:,i),wtilde)*dot(A*V(:,i-1),wbar))/...
            (dot(V(:,i-1),wbar)*dot(A*V(:,i),wbar)-dot(V(:,i),wbar)*dot(A*V(:,i-1),wbar));
        tempw2 = tempw1*dot(V(:,i),wbar)/dot(A*V(:,i),wbar)-dot(V(:,i),wtilde)/dot(A*V(:,i),wbar);
        
        vhat = -tempv2*A*vbar+tempv1*vbar-vtilde;
        what = -tempw2*A'*wbar+tempw1*wbar-wtilde;
        
        
        [alpha,u] = NormalizeIter(vhat,what);
        V = [V,u*vhat];
        W = [W,alpha*what];
        d_i = u*tempv1;
        c_i = u*tempv2;
        gamma_i = alpha*tempw1;
        delta_i = alpha*tempw2;
        
        if(i==2)
            T(i-1,i) = ((lambdarand/betarand)*u')';
            Stil(i-1,i) = ((brand/lrand)*alpha')';
            S(i-1,i) = u;
            Ttil(i-1,i) = alpha;
        else
            if(Beta(i-2)==0)
                T(i-1,i) = (Lambda(i-2)*u')';
                S(i-1,i) = 0;
            else
                T(i-1,i) = ((Lambda(i-2)/Beta(i-2))*u')';
                S(i-1,i) = u;
            end
            if(L(i-2)==0)
                Stil(i-1,i) = (B(i-2)*alpha')';
                Ttil(i-1,i) = 0;
            else
                Stil(i-1,i) = ((B(i-2)/L(i-2))*alpha')';
                Ttil(i-1,i) = alpha;
            end
        end
        
        T(i,i) = d_i;
        T(i+1,i) = B(i);
        S(i,i) = c_i;
        S(i+1,i) = L(i);
        
        Stil(i,i) = gamma_i;
        Stil(i+1,i) = Lambda(i);
        Ttil(i,i) = delta_i;
        Ttil(i+1,i) = Beta(i);
end
