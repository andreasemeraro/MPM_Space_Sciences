
%calcolo costellazione ottimale
clc
clear all
close all
global R mu eps j

%dati marte
R = 3390; %km
mu = 42828; %km^3/s^2

%costanti
eps = 5*pi/180;

min_sat = 30; %min 30 per avere pct=80% 
max_sat = 100;

%input
j = 1;
N = min_sat:1:max_sat;
P = 5:2:size(N,2)-1;
inc = pi/180*(90:-1:70);
h = 300:100:2000;
F = 1;   % phasing=F*2*pi/N
%utiliziamo costellazione Walker

%inizializzazione
J1 = zeros(size(N,2),size(h,2),size(P,2),size(inc,2));
J3 = zeros(size(N,2),size(h,2),size(P,2),size(inc,2));
J5 = zeros(size(N,2),size(h,2),size(P,2),size(inc,2));
J7 = zeros(size(N,2),size(h,2),size(P,2),size(inc,2));

opp_lb = sqrt(mu)/pi.*2/(1000+R)^1.5;
opp_ub = sqrt(mu)/pi.*(N(end)*(P(end)-1))/(300+R)^1.5;


%cicli for per calcolare le properties della costellazione e riempire le
%matrici J1 J3 J5 J7. Al fondo dello script sono segnate le funzioni usate
for n=1:size(N,2)
    
    for k=1:size(h,2)
        
        for p=1:size(P,2)
            
            for i=1:size(inc,2)
            
                cov = funcov(N(n),h(k));
            
                dv_eol = funeol(N(n),h(k));
                
        
                 %se N/P non da resto zero allora il risultato non è valido
                 %e riempiamo le matrici con un -1 per segnalarlo.
                if mod(N(n),P(p)) == 0 
         
                    opp = funopp(N(n),P(p),h(k));
                
                    psi = funmissd(N(n),P(p),F,inc(i));
                  
                    bld = P(p);
            
                else
                    opp = -1;
                
                    psi = -1;
                
                    bld = -1;
            end
            
           
           
           if opp == -1 || psi == -1 || bld == -1
               J1(n,k,p,i) = -1;
               J3(n,k,p,i) = -1;
               J5(n,k,p,i) = -1;
               J7(n,k,p,i) = -1;
               
           else
              J1(n,k,p,i) = cov;    
              J3(n,k,p,i) = (log10(pi/sqrt(mu)*opp)-log10(pi/sqrt(mu)*opp_lb))/(log10(pi/sqrt(mu)*opp_ub)-log10(pi/sqrt(mu)*opp_lb))+(pi/2-psi)/(pi/2);
              J5(n,k,p,i) = bld;
              J7(n,k,p,i) = log10(dv_eol);
           
           end
           
           
        end
          
     end 
    end
end
 
%calcoliamo il max e min di ogni matrice e le riscaliamo per avere gli
%stessi ordini di grandezza
lim_J1 = [max(J1,[],'all') min(J1(J1(:)>0),[],'all')];
lim_J3 = [max(J3,[],'all') min(J3(J3(:)>0),[],'all')];
lim_J5 = [max(J5,[],'all') min(J5(J5(:)>0),[],'all')];
lim_J7 = [max(J7,[],'all') min(J7(J7(:)>0),[],'all')];

J1 = (J1-lim_J1(2))/(lim_J1(1)-lim_J1(2));
J3 = (J3-lim_J3(2))/(lim_J3(1)-lim_J3(2));
J5 = (J5-lim_J5(2))/(lim_J5(1)-lim_J5(2));
J7 = (J7-lim_J7(2))/(lim_J7(1)-lim_J7(2));

%calcoliamo la matrice J come media delle altre matrici
J = (J1 + J3 + J5 + J7)/4;

%con altri cicli for andiamo a ricercare il valor minimo di J

A=1; B=1; C=1; D=1;
for a=1:size(N,2)
    
    for b=1:size(h,2)
       
        for c=1:size(P,2)
            
            for d=1:size(inc,2)
          
            if (J(a,b,c,d) < J(A,B,C,D)) && (J(a,b,c,d)>0)
                
                A = a;
                B = b;
                C = c;
                D = d;
                
            end                  
                      
            
        end
    end
end
end

J(A,B,C,D)

sprintf("number of satellites, N: %d",N(A))
sprintf("constellation altitude, h: %d km",h(B))
sprintf("number of orbital planes, P: %d",P(C))
sprintf("inclination of the orbits, inc: %d",inc(D)*180/pi)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cov] = funcov(N,h)
%UNTITLED2 Summary of this function goes here
%
global R eps j

%theta = 180-50-90-eps
theta = acos(R/(R+h)*cos(eps))-eps;

cov = N*(1-cos(theta))./(2*j);


end

function [opp] = funopp(N,P,h)
%UNTITLED2 Summary of this function goes here
%
global mu R

opp = sqrt(mu)/pi*N*(P-1)/(h+R)^(3/2);

end

% function [psi] = funmissd(N,P,F)
% %UNTITLED2 Summary of this function goes here
% %
% global mu R
% 
% n = 1:1:(N/P);
% 
% m = 1:1:(P-1);
% 
% psi=2*pi;
% 
% for n=1:size(n)
%     for m=1:size(m)
%         
%         psi_new = acos(cos(n*F/2)^2 - sin(n*F/2)^2*cos(m*N/P));
%         
%         if psi_new < psi
%             psi = psi_new;
%         end
%          
% end
% end 
% 
% end
% 
function [psi] = funmissd(N,P,F,i)
%UNTITLED2 Summary of this function goes here
%
global mu R

n = 1:1:(N/P);

m = 1:1:(P-1);

psi=2*pi;

for n=1:size(n)
    for m=1:size(m)
        
        alpha = n*F/2 + atan(tan(m*N/P/2)*cos(i));
        cosbeta =cos(i)^2 +sin(i)^2*cos(m*N/P);
        psi_new = acos(cos(alpha)^2 - sin(alpha)^2*cosbeta);
        
        if psi_new < psi
            psi = psi_new;
        end
         
end
end 

end

function [dv_eol] = funeol(N,h)

global R mu


a_n =  h+R;
r_pd = 100+ R;

dv_eol = N*(sqrt(mu/a_n)-sqrt(2*mu/a_n-(2*mu/(a_n+r_pd))));

end

% function [dv_sk] = funsk(N,h)
% 
% global R mu
% 
% 
% a_n = h + R;
% a_f = h + R - 50;
% 
% dv_sk = N*(-sqrt(mu/a_f)+sqrt(2*mu/a_f-(2*mu/(a_n+a_f)))+sqrt(mu/a_n)-sqrt(2*mu/a_n - (2*mu/(a_n+a_f))));
% 
% end
