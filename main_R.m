clear

R = 3;
ct = 50000;
Rvec = [1: 0.5: 3];
for n = 1: length(Rvec)
    R = Rvec(n);
    parfor i = 1 : ct 
%%% method 2: using fmincon as the alternative 
            h11 = abs(complex(randn(1,1)/sqrt(2),randn(1,1)/sqrt(2)))^2;
            h22 = abs(complex(randn(1,1)/sqrt(2),randn(1,1)/sqrt(2)))^2;
            h2 = max(0.001,min(h11,h22))/2;
            h1 = max(h11,h22);

            nonlcon = @mycons;%(x,N,tm,m,P,hall);
            options = optimoptions('fmincon','Display', 'off','MaxFunctionEvaluations', 300000); %display off
            x0 = [1  ]';
            A = []; % No other constraints
            b = [];
            Aeq = [];
            beq = [];
            lb = [];
            ub = []; 
            x = fmincon(@(x) x+exp(3*R/2)/h2*sqrt(1+x*h2)/(1+x*h1)-exp(R/2)/h2*sqrt(x*h2+1)...
                ,x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,h1,h2,R),options);
            y = exp(3*R/2)/h2*sqrt(1+x*h2)/(1+x*h1)-exp(R/2)/h2*sqrt(1+x*h2);
            mu = 1/(x*h2+1); lam = sqrt(exp(R)/h2^2/mu);
            P = [x y; lam-1/mu/h2 lam-1/h2];
            Poma =  [(exp(R)-1)/h1 (exp(R)-1)/h2]  ;
  
            %U2: OMA U1: HNOMA
            Pnew=zeros(2,2);
            Pnew(1,1) = (exp(R)-1)/h1;
            mux = 1/(Pnew(1,1)*h2+1); lamx = sqrt(exp(R)/h2^2/mux);
            Pnew(2,1) = lamx - 1/h2 ;
            Pnew(2,2) = lamx - 1/mux/h2;

            Pomaitr(i) = sum(Poma);
            Pminconitr(i) = sum(sum(P));
            Panaitr(i) = sum(sum(Pnew));
             
    end
    Pomaall(n) = mean(Pomaitr);
    Pminconall(n) =  mean(Pminconitr);
    Panaall(n) = mean(Panaitr);
end
plot(Rvec,Pomaall,Rvec,Panaall,Rvec,Pminconall)  

function [c,ceq] = mycons(x,h1,h2,R)
 
c(1) =  x*h1+1-exp(R);
c(2) = -x(1); 
ceq = [];
 
end 