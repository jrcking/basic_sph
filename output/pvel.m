function pvel(n)
clear A
if(n<=9) p1='0000';
elseif(n<=99) p1='000';
elseif(n<=999) p1='00';
elseif(n<=9999) p1='0';
else p1='';
endif
fname=strcat('r',p1,num2str(n),'.out');
A=load(fname);
fname=strcat('v',p1,num2str(n),'.out');
B=load(fname);
subplot(2,1,1),plot(A(:,2),B(:,1),'k.');
xlabel('distance from wall');ylabel('streamwise velocity');
subplot(2,1,2),plot(A(:,2),B(:,2),'k.');
xlabel('distance from wall');ylabel('cross-stream velocity');
endfunction
