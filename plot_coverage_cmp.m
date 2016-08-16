figure(10);
d_x=[2 502 1002 1502 2002 2502 3002 3502 3875 4050];
d=[1 1 0.99 0.98 0.97 0.92 0.88 0.72 0.2 0];
c(1)=coverage_rec(2);c(2)=coverage_rec(502);c(3)=coverage_rec(1002);c(4)=coverage_rec(1502);c(5)=coverage_rec(2002);
c(6)=coverage_rec(2502);c(7)=coverage_rec(3002);c(8)=coverage_rec(3502);c(9)=coverage_rec(4002);c(10)=coverage_rec(4502);
c(11)=coverage_rec(last_round+1);
c=c/64;
c_x=[2 502 1002 1502 2002 2502 3002 3502 4002 4502 last_round+1];
e_x=[2 100 250 300 400 500 550]; % LEACH
e=[1 1 0.98 0.9 0.6 0.2 0]; 
f_x=[2 500 1000 1100 1150 1200]; %PEGASIS
f=[1 1 1 0.97 0.2 0];

plot(c_x,c,d_x,d,f_x,f,e_x,e);

xlabel('Number of rounds','fontsize',14);
ylabel('Sensing coverage ratio','fontsize',14);
