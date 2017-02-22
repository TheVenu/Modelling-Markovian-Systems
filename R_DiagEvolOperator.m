#Evolution operators Ue and Ug are directly calculated in their diagonal form.
#A vector storing the energy values are stored into temp. Then temp is used to 
#create a diagonal matrix. 

load('R_gstate.m','phig_z');
global sigma_x=[0 1; 1 0];
global sigma_z=[1 0;0 -1];
global sigma_y=[0 -i;i 0];
t_step=0.01;
J=1;
a=0;        #start point
b=2;        #end point
n=(abs(b-a))/t_step+1;       
count=0;
N=2;
He=0;Hg=0;
L=0;

lambda=0.1;
delta=0.1;
function e=ep(k,lambda,delta) #energy function 
  a=1;
  J=1;
  N=2;
  e=2*J*sqrt(1+(lambda+delta)*(lambda+delta)-2*(lambda+delta)*cos(k*a)); #CHANGE HERE
endfunction

function Ue=Ue(lambda,delta) #unitary evolution matrix for excited state
  N=2;t_step=0.01;
  
  for k=1:2^(N-1)
  temp(end+1)=exp(-i*t_step*ep(k,lambda,delta));  
  temp(end+1)=exp(-i*t_step*ep(k,lambda,delta));  
  
  endfor
  Ue=diag(temp);
endfunction


function Ug=Ug(lambda,delta) #unitary evolution matrix for excited state
  N=2;t_step=0.01;
  for k=1:2^(N-1)
  temp(end+1)=exp(-i*t_step*ep(k,lambda,0));  
  temp(end+1)=exp(-i*t_step*ep(k,lambda,0));  
  
  endfor
  Ug=diag(temp);
endfunction

#initial state. This is where phi_g should be introduced
phi_z=phig_z;

phie_z=[0 1];
for i=1:2^(N)-2
  if(mod(i,2)==1)
  phie_z(end+1)=0;
  else
  phie_z(end+1)=1;
  endif
  endfor


for r=1:n
    phi_g=phi_z;
    phi_e=phi_z;
    ug=Ug(a+(r-1)*t_step,0);
    ue=Ue(a+(r-1)*t_step,delta);
  for q=1:n
    zz(q,r)=abs(phi_g'*phi_e)^2;  
    phi_g=ug*phi_g;
    phi_e=ue*phi_e;
    endfor 
endfor 

#to plot the graph
x=linspace(a,b,n);
x=y=linspace(a,b,n); #x=lambda y=time
[xx,yy]=meshgrid(x,y);
figure;
mesh(xx,yy,zz);
xlabel("lambda");
ylabel("time");
title("LosEcho vs time Corrected Hamiltonian");

