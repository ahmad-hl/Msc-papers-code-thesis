clear
xm=100;
ym=100;
sink.x=0.5*xm;
sink.y=0.5*ym;
n=100

%Optimal Election Probability of a node
%to become cluster head
p=0.05;
P=0.05;
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.1;
%alpha
a=1;

INFINITY = 999999999999999;
rmax=3000;
do=sqrt(Efs/Emp);
Et=0;

figure(1);

for i=1:1:n
    % X Axis
    S1(i).xd=rand(1,1)*xm;
    S2(i).xd=S1(i).xd;
    S3(i).xd=S1(i).xd;
    XR1(i)=S1(i).xd;
    XR2(i)=S2(i).xd;
    XR3(i)=S3(i).xd;
    
    % Y Axis
    S1(i).yd=rand(1,1)*ym;
    S2(i).yd=S1(i).yd;
    S3(i).yd=S1(i).yd;
    YR1(i)=S1(i).yd;
    S1(i).G=0;
    YR2(i)=S2(i).yd;
    S2(i).G=0;
    YR3(i)=S3(i).yd;
    S3(i).G=0;

    %initially there are no cluster heads only nodes   
    temp_rnd0 = i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S1(i).E=Eo;
        S2(i).E=Eo;
        S3(i).E=Eo;
        S2(i).ENERGY=0;
        plot(S2(i).xd,S2(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1) 
        S1(i).E=Eo*(1+a);
        S3(i).E=Eo*(1+a);
        S2(i).E=Eo*(1+a);
        S2(i).ENERGY=1;
        plot(S2(i).xd,S2(i).yd,'+');
        hold on;
    end
    E3(i)= S3(i).E;
    Et=Et+E3(i);
    
    %initially there are no cluster heads only nodes
    S1(i).type='N';
    S2(i).type='N';
    S3(i).type='N';
end

S1(n+1).xd=sink.x;
S1(n+1).yd=sink.y;
S2(n+1).xd=sink.x;
S2(n+1).yd=sink.y;
S3(n+1).xd=sink.x;
S3(n+1).yd=sink.y;
plot(S2(n+1).xd,S2(n+1).yd,'x');
hold on;  

%%***** LEACH (Low Energy Adaptive Clustering Hierarchy) **********%%

figure(1);

countCHs1=0;
cluster1=1;
flag_first_dead1=0;
flag_teenth_dead1=0;
flag_all_dead1=0;
dead1=0;
first_dead1=0;
teenth_dead1=0;
all_dead1=0;
allive1=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS1=0;
packets_TO_CH1=0;
throughput1=0;

for r=0:1:rmax     
    r
    
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S1(i).G=0;
        S1(i).cl=0;
    end
  end
dead1=0;
for i=1:1:n
    if (S1(i).E<=0)
        dead1=dead1+1; 
        if (dead1==1)
           if(flag_first_dead1==0)
              first_dead1=r;
              flag_first_dead1=1;
           end
        end
        if(dead1==0.1*n)
           if(flag_teenth_dead1==0)
              teenth_dead1=r;
              flag_teenth_dead1=1;
           end
        end
        if(dead1==n)
           if(flag_all_dead1==0)
              all_dead1=r;
              flag_all_dead1=1;
           end
        end
    end
    if S1(i).E>0
        S1(i).type='N';
    end
end
STATISTICS.DEAD1(r+1)=dead1;
STATISTICS.ALLIVE1(r+1)=allive1-dead1;

countCHs1=0;
cluster1=1;
for i=1:1:n
   if(S1(i).E>0)
   temp_rand=rand;     
   if ( (S1(i).G)<=0)  
        if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs1=countCHs1+1;
            packets_TO_BS1=packets_TO_BS1+1;
            throughput1 = throughput1+1;
             S1(i).type='C';
            S1(i).G=round(1/p)-1;
            C1(cluster1).xd=S1(i).xd;
            C1(cluster1).yd=S1(i).yd;
           distance=sqrt( (S1(i).xd-(S1(n+1).xd) )^2 + (S1(i).yd-(S1(n+1).yd) )^2 );
            C1(cluster1).distance=distance;
            C1(cluster1).id=i;
            X1(cluster1)=S1(i).xd;
            Y1(cluster1)=S1(i).yd;
            cluster1=cluster1+1;
           
           distance;
            if (distance>do)
                S1(i).E=S1(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S1(i).E=S1(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
   end
 end 
end
STATISTICS.CLUSTERHEADS1(r+1)=countCHs1;
for c=1:1:cluster1-1
    x1(c)=0;
end
y1=0;
z1=0;
for i=1:1:n
   if ( S1(i).type=='N' && S1(i).E>0 )
     if(cluster1-1>=1)
       min_dis=Inf;
       min_dis_cluster=0;
       for c=1:1:cluster1-1
           temp=min(min_dis,sqrt( (S1(i).xd-C1(c).xd)^2 + (S1(i).yd-C1(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
               x1(c)=x1(c)+1;
           end
       end
       
            min_dis;
            if (min_dis>do)
                S1(i).E=S1(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S1(i).E=S1(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        
            S1(C1(min_dis_cluster).id).E = S1(C1(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
            packets_TO_CH1=packets_TO_CH1+1; 
            throughput1 = throughput1+1;
 
        S1(i).min_dis=min_dis;
        S1(i).min_dis_cluster=min_dis_cluster;
    else
        y1=y1+1;
        min_dis=sqrt( (S1(i).xd-S1(n+1).xd)^2 + (S1(i).yd-S1(n+1).yd)^2 );
            if (min_dis>do)
                S1(i).E=S1(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S1(i).E=S1(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS1=packets_TO_BS1+1;
            throughput1 = throughput1+1;
     end
  end
end

STATISTICS.PACKETS_TO_CH1(r+1)=packets_TO_CH1;
STATISTICS.PACKETS_TO_BS1(r+1)=packets_TO_BS1;
STATISTICS.THROUGHPUT1(r+1)=throughput1;
end


%%***** DEEC (Distributed Energy Efficient Clustering) **********%%

d1=0.765*xm/2;
K=sqrt(0.5*n*do/pi)*xm/d1^2;
d2=xm/sqrt(2*pi*K);
Er=4000*(2*n*ETX+n*EDA+K*Emp*d1^4+n*Efs*d2^2);

countCHs3=0;
cluster3=1;
flag_first_dead3=0;
flag_teenth_dead3=0;
flag_all_dead3=0;

dead3=0;
first_dead3=0;
teenth_dead3=0;
all_dead3=0;

allive3=n;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS3=0;
packets_TO_CH3=0;
throughput3=0;

for r=0:1:rmax     
    r
    
  if(mod(r, round(1/P) )==0)
    for i=1:1:n
        S3(i).G=0;
        S3(i).cl=0;
    end
  end
Ea=Et*(1-r/rmax)/n;
dead3=0;
for i=1:1:n
    if (S3(i).E<=0)
        dead3=dead3+1; 
        if (dead3==1)
           if(flag_first_dead3==0)
              first_dead3=r;
              flag_first_dead3=1;
           end
        end
        if(dead3==0.1*n)
           if(flag_teenth_dead3==0)
              teenth_dead3=r;
              flag_teenth_dead3=1;
           end
        end
        if(dead3==n)
           if(flag_all_dead3==0)
              all_dead3=r;
              flag_all_dead3=1;
           end
        end
    end
    if S3(i).E>0
        S3(i).type='N';
    end
end
STATISTICS.DEAD3(r+1)=dead3;
STATISTICS.ALLIVE3(r+1)=allive3-dead3;

countCHs3=0;
cluster3=1;
for i=1:1:n
 if Ea>0
 p(i)=P*n*S3(i).E*E3(i)/(Et*Ea);
 if(S3(i).E>0)
   temp_rand=rand;     
   if ( (S3(i).G)<=0)  
        if(temp_rand<= (p(i)/(1-p(i)*mod(r,round(1/p(i))))))
            countCHs3=countCHs3+1;
            packets_TO_BS3=packets_TO_BS3+1;
            throughput3 = throughput3+1;
            S3(i).type='C';
            S3(i).G=round(1/p(i))-1;
            C3(cluster3).xd=S3(i).xd;
            C3(cluster3).yd=S3(i).yd;
           distance=sqrt( (S3(i).xd-(S3(n+1).xd) )^2 + (S3(i).yd-(S3(n+1).yd) )^2 );
            C3(cluster3).distance=distance;
            C3(cluster3).id=i;
            X3(cluster3)=S3(i).xd;
            Y3(cluster3)=S3(i).yd;
            cluster3=cluster3+1;
           distance;
            if (distance>do)
                S3(i).E=S3(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S3(i).E=S3(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    
   end
   
 end 
 end
end
STATISTICS.CLUSTERHEADS3(r+1)=countCHs3;

for c=1:1:cluster3-1
    x3(c)=0;
end

for i=1:1:n
   if ( S3(i).type=='N' && S3(i).E>0 )
     if(cluster3-1>=1)
       min_dis=Inf;
       min_dis_cluster=0;
       for c=1:1:cluster3-1
           temp=min(min_dis,sqrt( (S3(i).xd-C3(c).xd)^2 + (S3(i).yd-C3(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
               x3(c)=x3(c)+1;
           end
       end
 
            min_dis;
            if (min_dis>do)
                S3(i).E=S3(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S3(i).E=S3(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
       
            S3(C3(min_dis_cluster).id).E = S3(C3(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
            packets_TO_CH3=packets_TO_CH3+1;
            throughput3 = throughput3+1;
            
      
        S3(i).min_dis=min_dis;
        S3(i).min_dis_cluster=min_dis_cluster;
     else
            min_dis=sqrt( (S3(i).xd-S3(n+1).xd)^2 + (S3(i).yd-S3(n+1).yd)^2 );
            if (min_dis>do)
                S3(i).E=S3(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S3(i).E=S3(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
           packets_TO_BS3=packets_TO_BS3+1;
           throughput3 = throughput3+1;
   end
  end
end

 STATISTICS.PACKETS_TO_CH3(r+1)=packets_TO_CH3;
 STATISTICS.PACKETS_TO_BS3(r+1)=packets_TO_BS3;
 STATISTICS.THROUGHPUT3(r+1)=throughput3;
end

%%***** SEP (Stable Election Protocol) **********%%
figure(1);
%counter for CHs
countCHs2=0;
%counter for CHs per round
rcountCHs2=0;
cluster2=1;
countCHs2;
rcountCHs2=rcountCHs2+countCHs2;
flag_first_dead2=0;
first_dead2=0;
flag_all_dead2=0;
all_dead2=0;
teenth_dead2=0;
flag_teenth_dead2=0;
allive2=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS2=0;
packets_TO_CH2=0;
throughput2=0;

for r=0:1:rmax
    r

  %Election Probability for Normal Nodes
  pnrm=( P/ (1+a*m) );
  %Election Probability for Advanced Nodes
  padv= ( P*(1+a)/(1+a*m) );
    
  %Operation for heterogeneous epoch
  if(mod(r, round(1/pnrm) )==0)
    for i=1:1:n
        S2(i).G=0;
        S2(i).cl=0;
    end
  end

 %Operations for sub-epochs
 if(mod(r, round(1/padv) )==0)
    for i=1:1:n
        if(S2(i).ENERGY==1)
            S2(i).G=0;
            S2(i).cl=0;
        end
    end
  end

 
hold off;

%Number of dead nodes
dead2=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

figure(2);

for i=1:1:n
    %checking if there is a dead node
    if (S2(i).E<=0)
        plot(S2(i).xd,S2(i).yd,'red .');
        dead2=dead2+1;
        if(S2(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S2(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        if (dead2==1)
           if(flag_first_dead2==0)
              first_dead2=r;
              flag_first_dead2=1;
           end
        end
        if(dead2==0.1*n)
           if(flag_teenth_dead2==0)
              teenth_dead2=r;
              flag_teenth_dead2=1;
           end
        end
        if(dead2==n)
           if(flag_all_dead2==0)
              all_dead2=r;
              flag_all_dead2=1;
           end
        end
        hold on;    
    end
    if S2(i).E>0
        S2(i).type='N';
        if (S2(i).ENERGY==0)  
        plot(S2(i).xd,S2(i).yd,'o');
        end
        if (S2(i).ENERGY==1)  
        plot(S2(i).xd,S2(i).yd,'+');
        end
        hold on;
    end
end
plot(S2(n+1).xd,S2(n+1).yd,'x');


STATISTICS.DEAD2(r+1)=dead2;
STATISTICS.DEAD_N(r+1)=dead_n;
STATISTICS.ALLIVE2(r+1)=allive2-dead2;
STATISTICS.DEAD_A(r+1)=dead_a;

countCHs2=0;
cluster2=1;
for i=1:1:n
   if(S2(i).E>0)
   temp_rand=rand;     
   if ( (S2(i).G)<=0)

 %Election of Cluster Heads for normal nodes
 if( ( S2(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )

            countCHs2=countCHs2+1;
            packets_TO_BS2=packets_TO_BS2+1;
            throughput2 = throughput2+1;
            S2(i).type='C';
            S2(i).G=100;
            C(cluster2).xd=S2(i).xd;
            C(cluster2).yd=S2(i).yd;
            %plot(S2(i).xd,S2(i).yd,'k*');
            
            distance=sqrt( (S2(i).xd-(S2(n+1).xd) )^2 + (S2(i).yd-(S2(n+1).yd) )^2 );
            C(cluster2).distance=distance;
            C(cluster2).id=i;
            X(cluster2)=S2(i).xd;
            Y(cluster2)=S2(i).yd;
            cluster2=cluster2+1;
            
            %Calculation of Energy dissipated
            distance;
            if (distance>do)
                S2(i).E=S2(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S2(i).E=S2(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    


 %Election of Cluster Heads for Advanced nodes
 if( ( S2(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
        
            countCHs2=countCHs2+1;
            packets_TO_BS2=packets_TO_BS2+1;
            throughput2 = throughput2+1;
            S2(i).type='C';
            S2(i).G=100;
            C(cluster2).xd=S2(i).xd;
            C(cluster2).yd=S2(i).yd;
            %plot(S2(i).xd,S2(i).yd,'k*');
            
            distance=sqrt( (S2(i).xd-(S2(n+1).xd) )^2 + (S2(i).yd-(S2(n+1).yd) )^2 );
            C(cluster2).distance=distance;
            C(cluster2).id=i;
            X(cluster2)=S2(i).xd;
            Y(cluster2)=S2(i).yd;
            cluster2=cluster2+1;
            
            %Calculation of Energy dissipated
            distance;
            if (distance>do)
                S2(i).E=S2(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S2(i).E=S2(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    
    end
   end
    STATISTICS.PACKETS_TO_BS2(r+1)=packets_TO_BS2;
end


STATISTICS.CLUSTERHEADS2(r+1)=cluster2-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S2(i).type=='N' && S2(i).E>0 )
     if(cluster2-1>=1)
       min_dis=sqrt( (S2(i).xd-S2(n+1).xd)^2 + (S2(i).yd-S2(n+1).yd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster2-1
           temp=min(min_dis,sqrt( (S2(i).xd-C(c).xd)^2 + (S2(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %Energy dissipated
        if(min_dis>0)
            S2(C(min_dis_cluster).id).E = S2(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
         packets_TO_CH2=packets_TO_CH2+1; 
         throughput2 = throughput2+1;
        end

       S2(i).min_dis=min_dis;
       S2(i).min_dis_cluster=min_dis_cluster;
           
   end
   end
 STATISTICS.PACKETS_TO_CH2(r+1)=packets_TO_CH2;
end
hold on;

countCHs2;
rcountCHs2=rcountCHs2+countCHs2;

STATISTICS.THROUGHPUT2(r+1)=throughput2;


%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems

%[vx,vy]=voronoi(X,Y);
%plot(X,Y,'r*',vx,vy,'b-');
% hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);

end


%Plot Result
figure(3);
r=0:rmax;
plot(r,STATISTICS.DEAD1,'-k',r,STATISTICS.DEAD2,'-r',r,STATISTICS.DEAD3,'-b');
title('\bf LEACH VS SEP VS DEEC');
legend('LEACH','SEP','DEEC');
xlabel('Rounds');
ylabel('Network Life Time');
figure(4);
plot(r,STATISTICS.ALLIVE1,'-k',r,STATISTICS.ALLIVE2,'-r',r,STATISTICS.ALLIVE3,'-b');
title('\bf LEACH VS SEP VS DEEC');
legend('LEACH','SEP','DEEC');
xlabel('Rounds');
ylabel('Network Life Time');
figure(5);
title('\bf LEACH VS SEP VS DEEC');
legend('LEACH','SEP','DEEC');
plot(r,STATISTICS.PACKETS_TO_BS1,'-k',r,STATISTICS.PACKETS_TO_BS2,'-r',r,STATISTICS.PACKETS_TO_BS3,'-b');
title('\bf LEACH VS SEP VS DEEC');
legend('LEACH','SEP','DEEC');
xlabel('Rounds');
ylabel('Packets To BS');
figure(6);
plot(r,STATISTICS.PACKETS_TO_CH1,'-k',r,STATISTICS.PACKETS_TO_CH2,'-r',r,STATISTICS.PACKETS_TO_CH3,'-b');
title('\bf LEACH VS SEP VS DEEC');
legend('LEACH','SEP','DEEC');
xlabel('Rounds');
ylabel('Packets To CHs');
figure(7);
plot(r,STATISTICS.THROUGHPUT1,'-k',r,STATISTICS.THROUGHPUT2,'-r',r,STATISTICS.THROUGHPUT3,'-b');
title('\bf LEACH VS SEP VS DEEC');
legend('LEACH','SEP','DEEC');
xlabel('Rounds');
ylabel('Throughput');
figure(8);
subplot(3,1,1);
plot(r,STATISTICS.CLUSTERHEADS1,'-k');
title('\bf LEACH VS SEP VS DEEC');
legend('LEACH');
xlabel('Rounds');
ylabel('Cluster Heads');
subplot(3,1,2);
plot(r,STATISTICS.CLUSTERHEADS2,'-r');
legend('SEP');
xlabel('Rounds');
ylabel('Cluster Heads');
subplot(3,1,3);
plot(r,STATISTICS.CLUSTERHEADS3,'-b');
legend('DEEC');
xlabel('Rounds');
ylabel('Cluster Heads');
figure(9);
xlabel(first_dead1);
ylabel(first_dead2);
title(first_dead3);
figure(10);
xlabel(all_dead1);
ylabel(all_dead2);
title(all_dead3);