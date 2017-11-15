% This is a function to load the particle position data from 
% ksph and plot it in 2d. It's only a few lines

function ppos(n)
clear A
if(n<=9) p1='0000';
elseif(n<=99) p1='000';
elseif(n<=999) p1='00';
elseif(n<=9999) p1='0';
else p1='';
endif
fname=strcat('r',p1,num2str(n),'.out');
A=load(fname);

plot(A(:,1),A(:,2),'k.');

endfunction


