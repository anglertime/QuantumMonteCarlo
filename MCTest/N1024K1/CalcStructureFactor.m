function SkHist = CalcStructureFactor(coords,maxN)
% Calculating the structure factor S(k)
global Box
N = length(coords);
n = 1:maxN; 
SkHist(1:maxN, 1:maxN, 1:maxN) = 0; 
legal_kx = 2*pi*n/Box(1);
legal_ky = 2*pi*n/Box(2);
legal_kz = 2*pi*n/Box(3);
for kx = 1:maxN 
   for ky = 1:maxN 
        for kz = 1:maxN 
             k=[legal_kx(kx),legal_ky(ky),legal_kz(kz)];
             a=rhoK(coords,k); 
             SkHist(kx,ky,kz) = a*conj(a)/N;
        end 
   end 
end

function rhok = rhoK(coords,k) 
% Calculating the Rho_K
rhok=0;
for ptcl1=1:length(coords)
    rhok=rhok+exp(i*dot(k,coords(ptcl1,:))); 
end
   


