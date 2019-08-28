%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Author: Rodolfo Aguirre
%%The University of Texas at El Paso
%%raguirre4@utep.edu
%%915-637-2822
%%THIS PROGRAM CALCULATES THE POSITION OF GRAIN BOUNDARIES OVER
%%TIME, CALCULATES THE VELOCITY, THEN PLOTS THE VELOCITY OVER TEMPERATURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    
%%CLEAN WORKSPACE, CLEAN COMMAND WINDOW AND CLOSE ANY GRAPH
clc 
clear all
%clearvars -except VelkT 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for enerfol = 1:3
% INITIALIZE COMPOSITIONS
wz = 0;
zb = 0;
us = 0;

Boltz = 8.6173303e-5;

%numfiles = 1266000; %Sigma3
%numfiles = 125000; %Sigma 7 800 angs my simulations
numfiles = 2500000; %Sigma 11
%numfiles = 1250000; %Sigma 3 no Anneal

step = 10;
ast = 0;
count = 0; 
cou = 0;
co = 0;
cow = 0;
cowbell = 0;
kk = 0;
% set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Location of the simulation files at different temperatures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_list = {'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1250', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1300', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1350', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1400', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1450', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1500', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1550', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1600', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1650', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1700', ...
'/Volumes/My Passport/LONGtimeScaleSimulationFiles/Sigma7-111-E0.76T1750'};
%%
 figure(3)
 set(gcf,'color','w');

for v = 1:length(folder_list)
    pathname = folder_list{v};
    cd(pathname);
    % v
% if v == 1 
%     numfiles =350000;
% else if v == 2
%         numfiles = 425000;
%     else if v == 3
%             numfiles = 350000;
%         else if v == 4
%                 numfiles = 400000;
%             else if v == 5 
%                     numfiles = 400000;
%                 else if v == 6
%                         numfiles = 375000;
%                     else if v == 7
%                             numfiles = 375000;
%                         else if v == 8
%                                 numfiles = 325000;
%                             else if v == 9
%                                     numfiles = 325000;
%                                 else if v == 10
%                                         numfiles = 375000;
%                                     else if v == 11
%                                             numfiles = 375000;
%                                         else if v == 12
%                                                 numfiles = 350000;
%                                             else if v == 13
%                                                     numfiles = 375000;
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end


%fid = importdata('C:\Users\Brandon\Dropbox\GB Mobility\GB Mobility Simulations\GB_mobility_Rudy_Copy\Sigma11-311\E1.0NewAI800Y\T1800\r.3000', ' ', 9);
%cd ('C:\Users\Brandon\Dropbox\GB Mobility\GB Mobility Simulations\GB_mobility_Rudy_Copy\Sigma11-311\E1.0withAnneal\T1800');

%for k = 40000:1000:numfiles %Sigma 3 in my simulations
for k = 0:25000:numfiles %Sigma 3 in my simulations
%for k = 25000:25000:numfiles %Sigma 3 no Anneal and Xiaowang's simulations
%for k = 40000:5000:numfiles
% kk = kk + 1;
kk = 1;
    myfilename = sprintf('r.%d', k);
    mydata{kk} = importdata(myfilename, ' ', 9);
  
%close all
%whos mydata
y = mydata{kk}.data(:,4);
%y = -1*y + 1196.35503; %FOR Sigma 3 - 111 double the size (~1600 Angs)
%y = -1*y + 1192.8; %FOR Sigma 11 (~1600 Angs)
%y = -1*y + 1212.031423; %FOR Sigma 7 (~1600 Angs)

%y = -1*y + 612.7005; %ESTE ES EL Y QUE ESTABA USANDO PARA TODOS LAS
%GRAFICAS DE ANTES DE QUE CAMBIARA A LAS SIMULACIONES MAS LARGAS EN Y
%y = -1*y + 612.7005; %for sigma11 (~800 Angs)my simulations
%y = -1*y + 607.29453; %for sigma3 (~800 Angs)my simulations
%y = -1*y + 606.267; %for sigma7 (~800 Angs)my simulations


%NEW LONGsimulations

y = -1*y + 1212.14; %FOR Sigma 7 (~1600 Angs)

uop = mydata{kk}.data(:,6);
nop = mydata{kk}.data(:,7);

matrix(:,:) = [y uop];

% figure(1)
% set(gcf,'color','w');
% hold on
% plot(matrix(:,1),matrix(:,2),'k.');
% ylabel('unscaled order parameter')
% xlabel('Distance (A)')
%   xlim([-1200 1200])
%   ylim([0 3.5])
% drawnow

%FOR SIGMA 7 and 11 ONLY
indices = find((matrix(:,2)<1.7));
matrix(indices,:) = [];
indices2 = find((matrix(:,2)>1.8));
matrix(indices2,:) = [];
indices3 = find((matrix(:,1)>450));
matrix(indices3,:) = [];

% %FOR SIGMA 3 ONLY
% indices = find((matrix(:,2)<1));
% matrix(indices,:) = [];
% indices2 = find((matrix(:,2)>1.5));
% matrix(indices2,:) = [];
% indices3 = find((matrix(:,1)>250));
% matrix(indices3,:) = [];

%FOR SIGMA 3 ONLY Y SIZE = ~1600
% indices = find((matrix(:,2)<1));
% matrix(indices,:) = [];
% indices2 = find((matrix(:,2)>1.5));
% matrix(indices2,:) = [];
% indices3 = find((matrix(:,1)>450));
% matrix(indices3,:) = [];



count = count + 1;


% figure(2)
% set(gcf,'color','w'); 
% %hold on 
% plot(matrix(:,1),matrix(:,2),'b.');
% ylabel('unscaled order parameter')
% xlabel('Distance (A)')
%  xlim([-1200 1200])
%  ylim([0 3.5])
%  drawnow
 
%COUNT THE NUMBER OF POINTS PER DELTA Y 
len = length(matrix);
for w = -400:0.1:400     %for w = -250:0.1:250  for Y = ~800   or   -400:0.1:400 for Y = ~1600   the resolution depends on the step of the for loop
    cowbell = cowbell + 1;
    for q = 1:len
        if matrix(q,1)< w+1 && matrix(q,1)>w
            cow = cow + 1;
        end
    end
    m(cowbell,:) = [w cow]; 

    cow = 0;
end
cowbell = 0;

%m = sortrows(m);
    
gb(:,:,count) = [m(:,1) m(:,2)];

pas = find(gb(:,2,count)>100,1,'first');
pas2 = find(gb(:,2,count)>110,1,'last');

posR = gb(pas,1,count);
posL = gb(pas2,1,count);

time = k*(0.004e-12)/(1e-9);

test = linspace(0,700,100);



% figure(10)
% set(gcf,'color','w');
% 
% %  hold on 
%   %plot(gb(:,1,count),gb(:,2,count),'k.');
%      plot(gb(:,1,count),gb(:,2,count),'k.',posR,test,'r.',posL,test,'r.'); 
%      plot(gb(:,1,count),gb(:,2,count),'k.',posL,test,'r.',posR,test,'r.');
%      
% %      plot(posR,test,'r.',posL,test,'r.'); 
% %      plot(posL,test,'r.',posR,test,'r.');
%   xlim([-500 500])
%   ylim([0 1500])
%   ylabel('Counts')
%   xlabel('Distance (A)')
%  
%     %drawnow


figure(3)
set(gcf,'color','w');

hold on

if v == 1   
plot(time,posR,'k.',time,posL,'k.','MarkerSize',15);
end 
if v == 2
  plot(time,posR,'k+',time,posL,'k+');
end
if v == 3
    plot(time,posR,'kx',time,posL,'kx');
end
if v == 4
    plot(time,posR,'ks',time,posL,'ks');
end
if v == 5
    plot(time,posR,'k^',time,posL,'k^');
end
if v == 6
    plot(time,posR,'kd',time,posL,'kd');
end
if v == 7
    plot(time,posR,'k<',time,posL,'k<');
end
if v == 8
    plot(time,posR,'k>',time,posL,'k>');
end
if v == 9
    plot(time,posR,'kp',time,posL,'kp');
end
if v == 10
    plot(time,posR,'kh',time,posL,'kh');
end
if v == 11
    plot(time,posR,'ko',time,posL,'ko');
end
% if v == 12
%     plot(time, posR,'kv',time,posL,'kv');
% end
% if v == 13
%     plot(time, posR,'k.',time,posL,'k.','MarkerSize',5);
% end

xlabel('Time (ns)','Interpreter','latex','FontSize',20)
ylabel('Distance (\AA)','Interpreter','latex','FontSize',20)
ylim([-400 400])
xlim([0 10])
 drawnow

vel(count,:,v) = [time posL posR];

    

c = 0;
c2 = 0;
index = 0;
index2 = 0;
posR = 0;
posL = 0;
cou = 0; 
co = 0;
indices = 0;
indices2 = 0;
indices3 = 0;
y = [];
uop = [];
matrix = [];
di = [];
cow = 0; 
m = [];
mydata{:,:} = [];
end

%% PLOT LOCAL VELOCITY OF EVERY TWO POINTS


% numofpoints = length(vel(:,1,v));
% 
% for vas = 1:numofpoints-1
%     tiime = vel(vas,1,v);
%     veLL = (abs(vel(vas+1,2,v)) - abs(vel(vas,2,v)))/(vel(vas+1,1,v) - vel(vas,1,v));
%     veRR = abs(vel(vas+1,3,v) - vel(vas,3,v))/(vel(vas+1,1,v) - vel(vas,1,v));
%     aveveL = (veLL + veRR)/2;
%    
%     localVel(vas,:) = [tiime aveveL];
%     %localVel(vas,:) = [tiime veLL veRR];
% end



% set(gcf,'color','w');
% figure(4)
% hold on
% if v == 1
% plot(localVel(:,1),localVel(:,2),'k.','MarkerSize',15);
% % plot(localVel(:,1),localVel(:,2),'k.',localVel(:,1),localVel(:,3),'k.','MarkerSize',15);
% end 
% if v == 2
%  plot(localVel(:,1),localVel(:,2),'k+');
% %  plot(localVel(:,1),localVel(:,2),'k+',localVel(:,1),localVel(:,3),'k+');
% end
% if v == 3
%     plot(localVel(:,1),localVel(:,2),'kx');
% %    plot(localVel(:,1),localVel(:,2),'kx',localVel(:,1),localVel(:,3),'kx');
% end
% if v == 4
%      plot(localVel(:,1),localVel(:,2),'ks');
% %     plot(localVel(:,1),localVel(:,2),'ks',localVel(:,1),localVel(:,3),'ks');
% end
% if v == 5
%     plot(localVel(:,1),localVel(:,2),'k^');
% %    plot(localVel(:,1),localVel(:,2),'k^',localVel(:,1),localVel(:,3),'k^');
% end
% if v == 6
%      plot(localVel(:,1),localVel(:,2),'kd');
% %    plot(localVel(:,1),localVel(:,2),'kd',localVel(:,1),localVel(:,3),'kd');
% end
% if v == 7
%   plot(localVel(:,1),localVel(:,2),'k<');   
% %     plot(localVel(:,1),localVel(:,2),'k<',localVel(:,1),localVel(:,3),'k<');
%  
% end
% if v == 8
%     plot(localVel(:,1),localVel(:,2),'k>');
% %     plot(localVel(:,1),localVel(:,2),'k>',localVel(:,1),localVel(:,3),'k>');
% end
% xlabel('Time (ns)','Interpreter','latex','FontSize',20)
% ylabel('Velocity (\AA/ns)','Interpreter','latex','FontSize',20)
% ylim([-700 700])
%  drawnow




 numofpoints = 0;
 %localVel(:,:) = [];
 count = 0;
end
%% PLOT DISTANCE VS TIME SEPARATE FROM PREVIOUS CODE
close all


figure(100)
set(gcf,'color','w');
for U = 1:11
    
hold on

if U == 1   
plot(vel(:,1,U),vel(:,2,U),'k.',vel(:,1,U),vel(:,3,U),'k.','MarkerSize',15);
end 
if U == 2
  plot(vel(:,1,U),vel(:,2,U),'k+',vel(:,1,U),vel(:,3,U),'k+');
end
if U == 3
    plot(vel(:,1,U),vel(:,2,U),'kx',vel(:,1,U),vel(:,3,U),'kx');
end
if U == 4
    plot(vel(:,1,U),vel(:,2,U),'ks',vel(:,1,U),vel(:,3,U),'ks');
end
if U == 5
    plot(vel(:,1,U),vel(:,2,U),'k^',vel(:,1,U),vel(:,3,U),'k^');
end
if U == 6
    plot(vel(:,1,U),vel(:,2,U),'kd',vel(:,1,U),vel(:,3,U),'kd');
end
if U == 7
    plot(vel(:,1,U),vel(:,2,U),'k<',vel(:,1,U),vel(:,3,U),'k<');
end
if U == 8
    plot(vel(:,1,U),vel(:,2,U),'k>',vel(:,1,U),vel(:,3,U),'k>');
end
if U == 9
    plot(vel(:,1,U),vel(:,2,U),'kp',vel(:,1,U),vel(:,3,U),'kp');
end
if U == 10
    plot(vel(:,1,U),vel(:,2,U),'kh',vel(:,1,U),vel(:,3,U),'Kh');
end
if U == 11
    plot(vel(:,1,U),vel(:,2,U),'ko',vel(:,1,U),vel(:,3,U),'ko');
end


xlabel('Time (ns)','Interpreter','latex','FontSize',20)
ylabel('Distance (\AA)','Interpreter','latex','FontSize',20)
%ylim([0 75])
%xlim([0 0.5])
%drawnow


end

U = 0;





%% LINEAR REGRESION TO CALCULATE SLOPE (VELOCITY) OF EACH GB

figure(100)

for v = 1:length(folder_list)
    n = length(vel);
    spL = sum(vel(:,1,v).*vel(:,2,v));
    spR = sum(vel(:,1,v).*vel(:,3,v));
    sxs = sum(vel(:,1,v).^2);
    sx = sum(vel(:,1,v));
    sL = sum(vel(:,2,v));
    sR =  sum(vel(:,3,v));
    x = sx./n;
    yy1 = sR./n;
    yy2 = sL./n;
    
    slopeR = (n.*spR - sx.*sR)./(n.*sxs - (sx).^2);
    slopeL = (n.*spL - sx.*sL)./(n.*sxs - (sx).^2);
    constR = yy1 - slopeR*x;
    constL = yy2 - slopeL*x;
    
    fittime = linspace(0.2,0.8,1000);
    fitR = slopeR*fittime + constR;
    fitL = slopeL*fittime + constL;
    
    slopes11(v,:,3) = [slopeR slopeL];
    
    plot(fittime,fitR,'k',fittime,fitL,'k');
    hold on 
    
    if v == 1
        T = 1200;
    else if v == 2
            T = 1600;
        else if v == 3
                T = 1700;
            else if v ==4
                    T = 1800;
                else if v == 5
                        T = 1900;
                    else if v == 6
                            T = 2000;
%                         else if v == 7
%                                 T = 2100;
%                             else if v == 8
%                                     T = 2200;
%                                  else if v == 9
%                                          T = 1650;
%                                      else if v == 10
%                                              T = 1700;
%                                          else if v == 11
%                                                  T = 1750;
%                                              else if v == 12
%                                                      T = 1800;
%                                                  else if v == 13
%                                                          T = 1850;
%                                                      end
%                                                  end
%                                              end
%                                          end
%                                      end
%                                 end
%                             end
                        end
                    end
                end
            end
        end
    end
    
    

    VelkT(v,:,2) = [1/(Boltz*T) (abs(slopeR) + abs(slopeL))/2];





%vel = [];
 tiime = 0;
 veLL = 0;
 veRR = 0;
 vas = 0;
 mydata = {};
 kk = 0;
end
%end 


%% LINEAR REGRESSION TO CALCULATE THE SLOPE OF THE AVERAGE VELOCITY BETWEEN GU AND GD

figure(5)
%hold on
set(gcf,'color','w');

for tum = 1:1
losq = log(VelkT(:,2,tum));
num = length(VelkT(:,2,tum));
sop = sum(VelkT(:,1,tum).*losq);
soxs = sum(VelkT(:,1,tum).^2);
sox = sum(VelkT(:,1,tum));
soy = sum(losq);
xx = sox/num;
yy = soy/num;

slp = (num.*sop - sox.*soy)./(num.*soxs - (sox).^2);
constant = yy - slp*xx;

fitty = linspace(5,10,1000);
fits = slp*fitty + constant;
%slp
semilogy(VelkT(:,1,1),VelkT(:,2,1),'ks');
%semilogy(VelkT(:,1,1),VelkT(:,2,1),'ks',VelkT(:,1,2),VelkT(:,2,2),'k*',VelkT(:,1,3),VelkT(:,2,3),'k^');
hold on
% plot(fitty,fits,'k');
end

% plot(VelkT(:,1,1),log(VelkT(:,2,1)),'ks');
% plot(VelkT(:,1,2),log(VelkT(:,2,2)),'k*');
% plot(VelkT(:,1,3),log(VelkT(:,2,3)),'k^');

%semilogy(VelkT(:,1,1),VelkT(:,2,1),'ks',VelkT(:,1,2),VelkT(:,2,2),'k*',VelkT(:,1,3),VelkT(:,2,3),'k^');
%semilogy(VelkT(:,1,1),VelkT(:,2,1),'ks');
ylim([-100 1000]);
%xlim([5 10]);
hold on
%plot(fitty,fits,'k');
%semilogy(VelkT(:,1,2),VelkT(:,2,2),'b*');
%semilogy(VelkT(:,1,3),VelkT(:,2,3),'r*');
ylabel('ln(v)(ln\AA/ns)','Interpreter','latex','FontSize',20)
%ylabel('velocity(\AA/ns)','Interpreter','latex','FontSize',20)
xlabel('1/kT','Interpreter','latex','FontSize',20)


 %% CODE TO EXTRACT DATA OF EACH TEMPERATURE TO TXT FILES
 shas = 1200;

for o=1:11
    shas = shas + 50;
     
%      if (shas > 1200 && shas < 1600)
%          continue
%      end
    
    A = vel(:,:,o);
    save(sprintf('BothGBs_E0.76_T_%d.txt', shas), 'A', '-ASCII', '-double');
end
 

%% CODE TO PLOT THE INSTANT VELOCITY (VELOCITY PER EVERY TWO POINTS)
     %v = 13;
    
% songui = length(vel);
% 
% figure(10)
% set(gcf,'color','w');
% 
% 
% 
% 
% 
% for v = 5:5    
%     for ok = 1:songui - 1
%     
%         
%         instanvel = (vel(ok+1,3,v) - vel(ok,3,v))/(vel(ok+1,1,v) - vel(ok,1,v));   
%         instanvel2 = (vel(ok+1,2,v) - vel(ok,2,v))/(vel(ok+1,1,v) - vel(ok,1,v));
%         
%         
%         [time instavel]
%         
%        
% hold on
% if v == 1
%     plot(vel(ok+1,1,v),instanvel,'k.','MarkerSize',15);
%     plot(vel(ok+1,1,v),instanvel2,'r.','MarkerSize',15)
% %plot(time,posR,'k.',time,posL,'k.','MarkerSize',15);
% end 
% if v == 2
%     plot(vel(ok+1,1,v),instanvel,'k+');
%     plot(vel(ok+1,1,v),instanvel2,'r+');
%   %plot(time,posR,'k+',time,posL,'k+');
% end
% if v == 3
%     plot(vel(ok+1,1,v),instanvel,'kx');
%      plot(vel(ok+1,1,v),instanvel2,'rx');
%    % plot(time,posR,'kx',time,posL,'kx');
% end
% if v == 4
%     plot(vel(ok+1,1,v),instanvel,'ks');
%      plot(vel(ok+1,1,v),instanvel2,'rs');
%    % plot(time,posR,'ks',time,posL,'ks');
% end
% if v == 5
%     plot(vel(ok+1,1,v),instanvel,'k^');
%     plot(vel(ok+1,1,v),instanvel2,'r^');
%    % plot(time,posR,'k^',time,posL,'k^');
% end
% if v == 6
%     plot(vel(ok+1,1,v),instanvel,'kd');
%      plot(vel(ok+1,1,v),instanvel2,'rd');
%     %plot(time,posR,'kd',time,posL,'kd');
% end
% if v == 7
%     plot(vel(ok+1,1,v),instanvel,'k<');
%      plot(vel(ok+1,1,v),instanvel2,'r<');
%    % plot(time,posR,'k<',time,posL,'k<');
% end
% if v == 8
%     plot(vel(ok+1,1,v),instanvel,'k>');
%      plot(vel(ok+1,1,v),instanvel2,'r>');
%     %plot(time,posR,'k>',time,posL,'k>');
% end
% if v == 9
%     plot(vel(ok+1,1,v),instanvel,'kp');
%      plot(vel(ok+1,1,v),instanvel2,'rp');
%     %plot(time,posR,'kp',time,posL,'kp');
% end
% if v == 10
%     plot(vel(ok+1,1,v),instanvel,'kh');
%      plot(vel(ok+1,1,v),instanvel2,'rh');
%     %plot(time,posR,'kh',time,posL,'kh');
% end
% if v == 11
%     plot(vel(ok+1,1,v),instanvel,'ko');
%      plot(vel(ok+1,1,v),instanvel2,'ro');
%     %plot(time,posR,'ko',time,posL,'ko');
% end
% if v == 12
%     plot(vel(ok+1,1,v),instanvel,'kv');
%      plot(vel(ok+1,1,v),instanvel2,'rv');
%    % plot(time, posR,'kv',time,posL,'kv');
% end
% if v == 13
%     plot(vel(ok+1,1,v),instanvel,'k.','MarkerSize',5);
%     plot(vel(ok+1,1,v),instanvel,'r.','MarkerSize',5);
%     %plot(time, posR,'k.',time,posL,'k.','MarkerSize',5);
% end
% 
% xlabel('Time (ns)','Interpreter','latex','FontSize',20)
% ylabel('Velocity (\AA/ns)','Interpreter','latex','FontSize',20)
% ylim([-400 200])
% xlim([0 10])
%  drawnow

      
 
%     end
%     
% 
% end



%%
%  figure(12)
%  set(gcf,'color','w');
%  
% 
% plot(act(:,1),act(:,2),'-ks');
% 
% ylabel('Ea - E*c','Interpreter','latex','FontSize',20)
% xlabel('Grain Boundary Type (\Sigma)','FontSize',15)
% 
% 
% ylim([0 12]);
% xlim([0 12]);



%% CODE TO SHIFT THE GRAPHS FOR NON-LINEAR FITTING
%%WE WANT TO FIT THE FOLLOWING EQUATION: y = a0(1-exp(-a1x))
% 
for p = 2:6
vel(:,1,p) = vel(:,1,p) - abs(vel(1,1,p));
vel(:,3,p) = vel(:,3,p) + abs(vel(1,3,p));
vel(:,2,p) = vel(:,2,p) - vel(1,2,p);
vel(:,2,p) = -1*vel(:,2,p);
end


%% CODE FOR FITTING NON-LINEAR CURVE TO THE DATA
% 
for L = 1:6

%INITIAL GUESSES POSITIVE
azero = 65 ; %a0
aone = 1.5  ; %a1

%INITIAL GUESSES NEGATIVE
azero2 = -110; %a0 Initial guess for Down Grain Boundary
aone2 = 1.5;   %a1 Initial guess for Down Grain Boundary 



for t = 1:100
    
    
      for r = 1:length(vel)


         %PARTIAL DERIVATIVES POSITIVE
         pone = 1 - exp(-1*aone*vel(r,1,L));
         ptwo = azero*vel(r,1,L)*exp(-1*aone*vel(r,1,L));

         zzero(r,:) = [pone,ptwo];

         fofx = azero*(1 - exp(-1*aone*vel(r,1,L)));

         D(r,:) = [vel(r,3,L) - fofx];
         
         %PARTIAL DERIVATIVES NEGATIVE
         pone2 = 1 - exp(-1*aone2*vel(r,1,L));
         ptwo2 = azero2*vel(r,1,L)*exp(-1*aone2*vel(r,1,L));

         zzero2(r,:) = [pone2,ptwo2];

         fofx2 = azero2*(1 - exp(-1*aone2*vel(r,1,L)));

         D2(r,:) = [vel(r,2,L) - fofx2];        

      end

      %POSITIVE
      
zot = transpose(zzero);
mone = zot*zzero;
inmone = inv(mone);

mtwo = zot*D;

deltaA = inmone*mtwo;

azero = azero + deltaA(1,1);
aone = aone + deltaA(2,1);

   %NEGATIVE
   zot2 = transpose(zzero2);
mone2 = zot2*zzero2;
inmone2 = inv(mone2);

mtwo2 = zot2*D2;

deltaA2 = inmone2*mtwo2;

azero2 = azero2 + deltaA2(1,1);
aone2 = aone2 + deltaA2(2,1);

end
%as(L,:) = [azero aone azero2 aone2];

   coeffi(L,:) = [azero aone];
     coeffi2(L,:) = [azero2 aone2];
     
figure(32)
set(gcf,'color','w');
hold on
oneline = linspace(0,0.35,1000);
fofxx = coeffi(L,1)*(1 - exp(-1*coeffi(L,2)*oneline));
fofxx2 = coeffi2(L,1)*(1 - exp(-1*coeffi2(L,2)*oneline));

%plot(vel(:,1,L),vel(:,3,L),'k+',oneline,fofxx);

if L == 1   

    plot(vel(:,1,L),vel(:,3,L),'k.',oneline,fofxx);
    plot(vel(:,1,L),vel(:,2,L),'k.',oneline,fofxx2);
%plot(oneline,fofxx,oneline,fofxx2);

end 
if L == 2
     plot(vel(:,1,L),vel(:,3,L),'k+',oneline,fofxx);
     plot(vel(:,1,L),vel(:,2,L),'k+',oneline,fofxx2);
%plot(oneline,fofxx,oneline,fofxx2);
end
if L == 3
    plot(vel(:,1,L),vel(:,3,L),'kx',oneline,fofxx);
    plot(vel(:,1,L),vel(:,2,L),'kx',oneline,fofxx2);
 % plot(oneline,fofxx,oneline,fofxx2);
end
if L == 4
     plot(vel(:,1,L),vel(:,3,L),'ks',oneline,fofxx);
     plot(vel(:,1,L),vel(:,2,L),'ks',oneline,fofxx2);
 %plot(oneline,fofxx,oneline,fofxx2);
end
if L == 5
    plot(vel(:,1,L),vel(:,3,L),'k^',oneline,fofxx);
    plot(vel(:,1,L),vel(:,2,L),'k^',oneline,fofxx2);
 %plot(oneline,fofxx,oneline,fofxx2);
end
if L == 6
     plot(vel(:,1,L),vel(:,3,L),'kd',oneline,fofxx);
     plot(vel(:,1,L),vel(:,2,L),'kd',oneline,fofxx2);
 %plot(oneline,fofxx,oneline,fofxx2);
end

   end
ylim([-200 200])
xlim([0 0.35])

xlabel('Time (ns)','Interpreter','latex','FontSize',20)
ylabel('Distance (\AA)','Interpreter','latex','FontSize',20)

%% CODE FOR FITTING NON-LINEAR CURVE TO THE DATA with NEW TERMS (REGRESSION)
%close all
% a = 0;
% figure(100)
% 
%     if (T > 1200 && T < 1600)
%         continue
%     end
%     a = a + 1;


close all 
clc
clearvars -except vel

A = 4000;
Ea = 1.92;
c = 2.02;
E = 0.75;
k = 8.6173303e-5;
aone = 5;
x = linspace(0,0.35,1000);
a = 10^-3; 
T = 1600;




 for L = 2:2
    
for t = 1:2
    
    
      for r = 1:length(vel)        
          
  %DEFINITIONS
         kT = k*T;
         R = -1*(Ea - c*E)/kT;
         Ti = -1*aone*vel(r,1,L);
         %Ti = -1*aone*x;
  %PARTIAL DERIVATIVES POSITIVE
         pone = (exp(R)/aone) - (exp(R+Ti)/aone);
         ptwo = (a/aone^2) - ((A*exp(R))/aone^2) + (A*(-1*Ti + 1)*exp(R+Ti)/aone^2) ...
             - (A*(-1*Ti + 1)*exp(Ti)/aone^2);
         pthree = (A*E*exp(R)/(aone*kT)) - (A*E*exp(R+Ti)/(aone*kT));
         pfour = (-1*A*exp(R)/(aone*kT)) + (A*exp(R+Ti)/(aone*kT));

         zzero(r,:) = [pone,ptwo,pthree,pfour];

         fofx = (a*vel(r,1,L)) - (a/aone) + (A*exp(R)/aone) - (A*exp(R+Ti)/aone) + ...
            (a*exp(Ti)/aone);
          
        %fofx = (a*x) - (a/aone) + (A*exp(R)/aone) - (A*exp(R+Ti)/aone) + ...
        %(a*exp(Ti)/aone);         
                  
        D(r,:) = vel(r,3,L) - fofx;
     
      end
%POSITIVE
      
zot = transpose(zzero);
mone = zot*zzero;
inmone = inv(mone);

mtwo = zot*D;

deltaA = inmone*mtwo;

A = A + deltaA(1,1);
aone = aone + deltaA(2,1);
c = c + deltaA(3,1);
Ea = Ea + deltaA(4,1);

end

figure(32)
set(gcf,'color','w');
hold on
oneline = linspace(0,0.35,1000);

if L == 2   
plot(vel(:,1,L),vel(:,3,L),'k.',x,fofx,'k.');
xlabel('Time (ns)','Interpreter','latex','FontSize',20)
ylabel('Distance (\AA)','Interpreter','latex','FontSize',20)
end 
 end

 
 
 