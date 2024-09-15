clear
close all

R = 3;
%%% method 2: using fmincon as the alternative 
h11 = abs(complex(randn(1,1)/sqrt(2),randn(1,1)/sqrt(2)))^2;
h22 = abs(complex(randn(1,1)/sqrt(2),randn(1,1)/sqrt(2)))^2;
h1 = max(h11,h22);
h2 = min(h11,h22);

Poma =  [(exp(R)-1)/h1 (exp(R)-1)/h2]  ; 
           
stepz = Poma(1)/300;            
z1vec = [0:stepz: 2*Poma(1)];
z2vec = [0:stepz: 2*Poma(1)]; 
Pex = zeros(length(z1vec),length(z2vec)); 
for i = 1 : length(z1vec)
    for k = 1 : length(z2vec)
        Pomaplan(i,k) = Poma(1);
        z1 = z1vec(i);
        z2 = z2vec(k);
        if log(1+z1*h1)+log(1+ exp(-R/2)*h2*z2/sqrt(z1*h2+1))>=R & z1*h2+1<=exp(R)
            Pex(i,k) = z1+z2;
        else
            Pex(i,k) = inf;
        end
    end
end
[id1,id2] = find(Pex==min(Pex(:)));
%[z1vec(id1) z2vec(id1)]
z1 = z1vec(id1); z2 = z2vec(id2);
[min(min(Pex)) Poma(1)]

mesh(Pex)

hold
alpha 0.5

mesh(Pomaplan)
 