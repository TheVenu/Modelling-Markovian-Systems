#Calculating the ground state then save to R_gstate and use in DiagEvolOperator
#The Ak is calculated according to wigner mapping tensored with identities in Aknew
global sigma_x=[0 1; 1 0];
global sigma_z=[1 0;0 -1];
global sigma_y=[0 -i;i 0];
t_step=0.01;
J=1;
a=0;        #start point
b=2;        #end point
n=(abs(b-a))/t_step+1;       

N=2;
He=0;Hg=0;
L=0;

delta=0.1;
alpha=0.5;
theta=pi/4;



function Aknew=aknew(k,lambda,delta)
  sigma_x=[0 1; 1 0];
  sp=[0 1; 0 0];
  sm=[0 0; 1 0];
  sigma_plus=1/2*[-1 1; -1 1];
  sigma_minus=1/2*[-1 -1; 1 1];
  N=2;
  Aknew=0;
  s=cos(thetae(k,lambda,delta)/2)*sigma_plus-i*sin(thetae(k,lambda,delta)/2)*sigma_minus;
  j=1;
  while(j<=N)  
    temp=1;
    f=1;
    while(f<j) 
    temp=kron(temp,sigma_x);
    f++;
    endwhile
    temp=kron(temp,s);
    
    f=j+1;
    while (f<=N)
    temp=kron(temp,eye(2));
    f++;
    endwhile
      
    Aknew+=exp(-i*k*j)*temp/sqrt(2^N);
    j++;
  endwhile
 
endfunction



function q=alp(k,lambda,delta)
  q=(thetae(k,lambda,0)-thetae(k,lambda,delta))/2;
endfunction


function t=thetae(k,lambda,delta)
  a=1;
  t=atan(-sin(a*k)/(cos(a*k)-(lambda+delta)));
endfunction


function e=ep(k,lambda,delta)
  a=1;
  J=1;
  e=2*J*sqrt(1+(lambda+delta)*(lambda+delta)-2*(lambda+delta)*cos(k*a)); #CHANGE HERE
endfunction




lambda=1;
delta=0.01;
k=4;

phie_z=[0 1];
for i=1:2^(N)-2
  if(mod(i,2)==1)
  phie_z(end+1)=0;
  else
  phie_z(end+1)=1;
  endif
  endfor


temp=eye(2^(N),2^(N));

for k=1:N
  Akplus=aknew(k,lambda,delta);
  Akminus=aknew(-k,lambda,delta);
  temp*=(cos(alp(k,lambda,delta))*eye(2^(N),2^(N))-i*sin(alp(k,lambda,delta))*Akplus'*Akminus');
endfor
phig_z=temp*phie_z';
save('R_gstate.m','phig_z');



