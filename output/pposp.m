% This is a function to load the particle position data from 
% ksph and plot it in 2d. It's only a few lines

function pposp(n)
clear A
if(n<=9) p1='0000';
elseif(n<=99) p1='000';
elseif(n<=999) p1='00';
elseif(n<=9999) p1='0';
else p1='';
endif
fname=strcat('r',p1,num2str(n),'.out');
A=load(fname);


fname=strcat('p',p1,num2str(n),'.out');
B=load(fname);
C=B(:,1);
colormap('jet');
scatter(A(:,1),A(:,2),15,C,'filled')
colorbar

endfunction


