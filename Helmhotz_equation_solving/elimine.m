function [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu)

Nbpt = size(AA,1);

tilde_AA = AA;
tilde_LL = LL;

for i=1:Nbpt
    if Refneu(i)==1 || Refneu(i)==2
        tilde_LL(i) = 0;
        tilde_AA(i,:) = 0;
        tilde_AA(:,i) = 0;
        tilde_AA(i,i) = 1;
    end
end
