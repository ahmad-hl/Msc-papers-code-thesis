clear
xm=100;
ym=100;
sink.x=xm*0.5;
sink.y=ym*0.5;
n=100;

%Optimal Election Probability of a node
%to become cluster head
p=0.1;
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

INFINITY = 999999999999999;
rmax=2000;
do=sqrt(Efs/Emp);
Et=0;

figure(1);
for i=1:1:n
    % X Axis
    S1(i).xd=rand(1,1)*xm;
    S2(i).xd=S1(i).xd;
    XR1(i)=S1(i).xd;
    XR2(i)=S2(i).xd;
    % Y Axis
    S1(i).yd=rand(1,1)*ym;
    S2(i).yd=S1(i).yd;
    S3(i).yd=S1(i).yd;
    YR1(i)=S1(i).yd;
    S1(i).G=0;
    YR2(i)=S2(i).yd;
    S2(i).G=0;
    
    %initially there are no cluster heads only nodes   
    S1(i).E=Eo;
    S2(i).E=S1(i).E;
    E2(i)=S2(i).E;
    if(sink.x>xm)
        if(sink.y>ym)
            axis([0 sink.x+10 0 sink.y+10])
        else
            axis([0 sink.x+10 0 ym])
        end
    else if(sink.y>ym)
            axis([0 xm 0 sink.y+10])
        else
            axis([0 xm 0 ym])
        end
    end
    plot(S2(i).xd,S2(i).yd,'o');
    hold on;
    Et=Et+E2(i);
    
    %initially there are no cluster heads only nodes
    S1(i).type='N';
    S2(i).type='N';
    S2(i).type='N';
end

if(sink.x>xm)
        if(sink.y>ym)
            axis([0 sink.x+10 0 sink.y+10])
        else
            axis([0 sink.x+10 0 ym])
        end
    else if(sink.y>ym)
            axis([0 xm 0 sink.y+10])
        else
            axis([0 xm 0 ym])
        end
end
    
S1(n+1).xd=sink.x;
S1(n+1).yd=sink.y;
S2(n+1).xd=sink.x;
S2(n+1).yd=sink.y;
plot(S2(n+1).xd,S2(n+1).yd,'X');
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
 
        S1(i).min_dis=min_dis;
        S1(i).min_dis_cluster=min_dis_cluster;
     else
        min_dis=sqrt( (S1(i).xd-S1(n+1).xd)^2 + (S1(i).yd-S1(n+1).yd)^2 );
            if (min_dis>do)
                S1(i).E=S1(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S1(i).E=S1(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS1=packets_TO_BS1+1;
     end
  end
end
STATISTICS.PACKETS_TO_CH1(r+1)=packets_TO_CH1;
STATISTICS.PACKETS_TO_BS1(r+1)=packets_TO_BS1;
 
%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems
%[vx,vy]=voronoi(X1,Y1);
%plot(X1,Y1,'r*',vx,vy,'b-');
%hold on;
%voronoi(X1,Y1);
%axis([0 xm 0 ym]);
end


%%***** DE LEACH (Distributed Energy For LEACH) **********%%
d1=0.765*xm/2;
K=sqrt(0.5*n*do/pi)*xm/d1^2;
d2=xm/sqrt(2*pi*K);
Er=4000*(2*n*ETX+n*EDA+K*Emp*d1^4+n*Efs*d2^2);

countCHs2=0;
cluster2=1;
flag_first_dead2=0;
flag_teenth_dead2=0;
flag_all_dead2=0;
dead2=0;
first_dead2=0;
teenth_dead2=0;
all_dead2=0;
allive2=n;
Etmp = Et;
Et
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS2=0;
packets_TO_CH2=0;
    
for r=0:1:rmax     
    r

  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S2(i).G=0;
        S2(i).cl=0;
    end
  end
  
  if(Etmp>0)
    Ea=Etmp*(1-r/rmax)/n;
    Etmp = Etmp - Er;
  end
  
dead2=0;
for i=1:1:n
    if (S2(i).E<=0)
        dead2=dead2+1; 
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
    end
    if S2(i).E>0
        S2(i).type='N';
    end
end
STATISTICS.DEAD2(r+1)=dead2;
STATISTICS.ALLIVE2(r+1)=allive2-dead2;

countCHs2=0;
cluster2=1;
for i=1:1:n
    %if(S2(n+1).xd - S2(i).xd <= 0 &&  S2(n+1).yd - S2(i).xd <=0)
    distance=sqrt( (S2(i).xd-(S2(n+1).xd) )^2 + (S2(i).yd-(S2(n+1).yd) )^2 );
    if (distance<=do)
        if(Ea>0)
            if(S2(i).E>0)
                temp_rand=rand;     
                if ( (S2(i).G)<=0)  
                    if(temp_rand<= (p/(1-p*mod(r,round(1/p)))) && S2(i).E>Ea)
                        countCHs2=countCHs2+1;
                        packets_TO_BS2=packets_TO_BS2+1;
                        S2(i).type='C';
                        S2(i).G=round(1/p)-1;
                        C2(cluster2).xd=S2(i).xd;
                        C2(cluster2).yd=S2(i).yd;
                        distance=sqrt( (S2(i).xd-(S2(n+1).xd) )^2 + (S2(i).yd-(S2(n+1).yd) )^2 );
                        C2(cluster2).distance=distance;
                        C2(cluster2).id=i;
                        X2(cluster2)=S2(i).xd;
                        Y2(cluster2)=S2(i).yd;
                        cluster2=cluster2+1;
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
        end
        S2(n+1).xd
    else
        if(S2(i).E>0)
            temp_rand=rand;     
            if ( (S2(i).G)<=0)    
                if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                    countCHs2=countCHs2+1;
                    packets_TO_BS2=packets_TO_BS2+1;
                    S2(i).type='C';
                    S2(i).G=round(1/p)-1;
                    C2(cluster2).xd=S2(i).xd;
                    C2(cluster2).yd=S2(i).yd;
                    distance=sqrt( (S2(i).xd-(S2(n+1).xd) )^2 + (S2(i).yd-(S2(n+1).yd) )^2 );
                    C2(cluster2).distance=distance;
                    C2(cluster2).id=i;
                    X2(cluster2)=S2(i).xd;
                    Y2(cluster2)=S2(i).yd;
                    cluster2=cluster2+1;
           
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
    end
end

STATISTICS.CLUSTERHEADS2(r+1)=countCHs2;

for i=1:1:n
   if ( S2(i).type=='N' && S2(i).E>0 )
     if(cluster2-1>=1)
       min_dis=Inf;
       min_dis_cluster=0;
       for c=1:1:cluster2-1
           temp=min(min_dis,sqrt( (S2(i).xd-C2(c).xd)^2 + (S2(i).yd-C2(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
 
            min_dis;
            if (min_dis>do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
       
            S2(C2(min_dis_cluster).id).E = S2(C2(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
            packets_TO_CH2=packets_TO_CH2+1;
      
        S2(i).min_dis=min_dis;
        S2(i).min_dis_cluster=min_dis_cluster;
     else
            min_dis=sqrt( (S2(i).xd-S2(n+1).xd)^2 + (S2(i).yd-S2(n+1).yd)^2 );
            if (min_dis>do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
           packets_TO_BS2=packets_TO_BS2+1;
   end
  end
end
 STATISTICS.PACKETS_TO_CH2(r+1)=packets_TO_CH2;
 STATISTICS.PACKETS_TO_BS2(r+1)=packets_TO_BS2;
 
 

end

%Plot Result
figure(2);
r=0:rmax;
plot(r,STATISTICS.DEAD1,'-k',r,STATISTICS.DEAD2,'-.r');
title('\bf (LEACH) VS (BS Distance Adaptive LEACH)');
legend('LEACH','BS-DA LEACH');
xlabel('Rounds');
ylabel('Network Life Time');
figure(4);
plot(r,STATISTICS.ALLIVE1,'-k',r,STATISTICS.ALLIVE2,'-.r');
title('\bf (LEACH) VS (BS Distance Adaptive LEACH)');
legend('LEACH','BS-DA LEACH');
xlabel('Rounds');
ylabel('Network Life Time');
figure(5);
title('\bf (LEACH) VS (BS Distance Adaptive LEACH)');
legend('LEACH','BS BS-DA LEACH');
plot(r,STATISTICS.PACKETS_TO_BS1,'-k',r,STATISTICS.PACKETS_TO_BS2,'-.r');
title('\bf (LEACH) VS (BS Distance Adaptive LEACH)');
legend('LEACH','BS-DA LEACH');
xlabel('Rounds');
ylabel('Packets To BS');
figure(6);
plot(r,STATISTICS.PACKETS_TO_CH1,'-k',r,STATISTICS.PACKETS_TO_CH2,'-.r');
title('\bf (LEACH) VS (BS Distance Adaptive LEACH)');
legend('LEACH','BS-DA LEACH');
xlabel('Rounds');
ylabel('Packets To CHs');
figure(7);

subplot(2,1,1);
plot(r,STATISTICS.CLUSTERHEADS1,'-k');
title('\bf (LEACH) VS (BS Distance Adaptive LEACH)');
legend('LEACH');
xlabel('Rounds');
ylabel('Cluster Heads');
subplot(2,1,2);
plot(r,STATISTICS.CLUSTERHEADS2,'-.r');
title('\bf (LEACH) VS (BS Distance Adaptive LEACH)');
legend('BS-DA LEACH');
xlabel('Rounds');
ylabel('Cluster Heads');
