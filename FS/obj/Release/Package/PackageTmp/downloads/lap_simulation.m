function lap_simulation()



%Outputs easy to read results    
format short g

%RPM/Torque Curves from DELFT 2007 WR450 Torque Curves
rpm = [500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4090 4155 4226 4284 4349 4405 4472 4531 4594 4657 4723 4786 4844 4916 4980 5046 5112 5179 5246 5313 5386 5459 5527 5603 5681 5760 5839 5925 6002 6088 6167 6246 6329 6410 6486 6565 6648 6732 6814 6897 6977 7064 7142 7234 7322 7402 7489 7573 7657 7739 7824 7902 7988 8065 8144 8220 8299 8377 8446 8521 8595 8667 8743 8813 8880 8952 9017 9084 9142 9210 9270 9332 9391 9452 9509 9567 9620 9683 9739 9794 9840 9894 9941 9987 10034 10057 10096 10106 10111 10101];
tor = [2.0054 2.3754 2.7454 3.1154 3.4854 3.8554 4.2254 4.5954 4.9654 5.3354 5.7054 6.0754 6.4454 6.8154 7.1854 7.5554 7.9254 8.2954 8.6654 9.0354 9.4054 9.7754 10.1454 10.5154 10.8854 11.2554 11.6254 11.9954 12.3654 12.7354 13.1054 13.4754 13.8454 14.2154 14.5854 14.9554 15.244 18.722 22.274 23.754 24.642 24.79 24.79 24.864 24.938 25.012 25.086 25.308 25.53 25.678 25.9 26.122 26.344 26.418 26.566 26.714 27.01 27.454 27.972 28.638 29.674 30.266 31.006 31.894 32.264 32.56 32.634 32.56 32.412 32.19 32.116 32.116 32.264 32.56 32.782 32.856 33.152 33.340 33.892 34.04 34.114 34.114 34.114 34.114 34.04 33.744 33.596 33.374 33.004 32.708 32.486 32.042 31.672 31.228 30.932 30.562 30.266 29.97 29.674 29.23 28.934 28.49 28.046 27.676 27.084 26.566 25.9 25.382 25.086 24.642 24.272 23.902 23.606 23.458 23.384 23.014 22.422 21.386 20.572 19.832 19.462 18.056 15.244 13.394 9.028 4.144];











%Calculating BSFC:
%BSFC=r/P
%P=power produced, P=ax*vehicle velocity
%r is the fuel consumption rate(g-s^-1)
%BSFC is a constant for most engines
%r*total time = fuel used.
%http://en.wikipedia.org/wiki/Brake_specific_fuel_consumption 
%Link states that the typical BSFC for gasoline engines is 0.45 to 0.37
%Taking the average of 0.45 and 0.37 yields a BSFC of .41
%BSFC is in units of lb/hp*h
%ax=ft/s2*ft/s

BSFC=.41/3600;












    
%Setting future variables to zero.
corner_total=0;
prev_corners=0;


%This whole section is commented out because it was deemend unnecessary,
%however in future versions of this program it could prove useful.  The
%intent here was to allow the user to manually input the number of corners
%and straight sections.  


%Input data including the number of straight sections, number of corners
%following those sections and their radii




%Input Section for User
%straight_sections=input('THE NUMBER OF STRAIGHT SECTIONS:  ');

%for i=1:straight_sections
%    dist(i)=input('STRAIGHT SECTION LENGTH(ft):  ');
%    number_of_corners(i)=input('INPUT THE NUMBER OF ASSOCIATED CORNERS AFTER STRAIGHT SECTION:  ');
    
%    for j=1:number_of_corners(i)
        
%        corner_radius(j)=input('INPUT THE RADII OF EACH CORNER(ft):  ');
%    end
    
%    j=0;
    
%    for k=(prev_corners+1):(number_of_corners(i)+prev_corners)
        
        
        
%        j=j+1;    
        
%        corner_total(k)=corner_radius(j);
        
%   end
        

 
%    prev_corners=number_of_corners(i)+prev_corners;
%    clear corner_radius
    
    
%end


%FSAE West 2007 Endurance Track Map
%Endurance Event is 13.66 miles, total distance/track distance will equal
%number of laps.


number_of_laps=12.4991;


%Straight Section Distances
straight_sections=18;
dist(1)=280;
dist(2)=50;
dist(3)=182;
dist(4)=109;
dist(5)=323;
dist(6)=440;
dist(7)=127;
dist(8)=94;
dist(9)=55;
dist(10)=415;
dist(11)=205;
dist(12)=62;
dist(13)=353;
dist(14)=16;
dist(15)=64.6;
dist(16)=281;
dist(17)=84.1;
dist(18)=28.8;


%Associated Corners
number_of_corners(1)=2;
corner_radius(1)=51.43;
cornering_dist(1)= 156.79;
corner_radius(2)=18.40;
cornering_dist(2)= 40.27;

number_of_corners(2)=1;
corner_radius(3)=79.09;
cornering_dist(3)= 75.64;

number_of_corners(3)=1;
corner_radius(4)=87.66;
cornering_dist(4)= 56.15;

number_of_corners(4)=1;
corner_radius(5)=82.51;
cornering_dist(5)= 206.37;

number_of_corners(5)=6;
corner_radius(6)=76.86;
cornering_dist(6)= 242.80;
corner_radius(7)=47.89;
cornering_dist(7)= 82.41;
corner_radius(8)=47.89;
cornering_dist(8)= 80.65;
corner_radius(9)=43.54;
cornering_dist(9)= 59.58;
corner_radius(10)=47.89;
cornering_dist(10)= 72.07;
corner_radius(11)=47.89;
cornering_dist(11)= 39.90;

number_of_corners(6)=2;
corner_radius(12)=74.51;
cornering_dist(12)= 38.63;
corner_radius(13)=140.29;
cornering_dist(13)= 81.29;

number_of_corners(7)=2;
corner_radius(14)=140.29;
cornering_dist(14)= 81.83;
corner_radius(15)=162.34;
cornering_dist(15)= 105.12;

number_of_corners(8)=2;
corner_radius(16)=74.29;
cornering_dist(16)= 124.2;
corner_radius(17)=80.51;
cornering_dist(17)= 236.36;

number_of_corners(9)=2;
corner_radius(18)=143.66;
cornering_dist(18)= 74.47;
corner_radius(19)=50.46;
cornering_dist(19)= 119.36;

number_of_corners(10)=2;
corner_radius(20)=68.40;
cornering_dist(20)= 82.14;
corner_radius(21)=59.71;
cornering_dist(21)= 70.1;

number_of_corners(11)=1;
corner_radius(22)=9;
cornering_dist(22)= 14.13;

number_of_corners(12)=1;
corner_radius(23)=9;
cornering_dist(23)= 11.14;

number_of_corners(13)=1;
corner_radius(24)=9;
cornering_dist(24)= 14.14;

number_of_corners(14)=1;
corner_radius(25)=28.11;
cornering_dist(25)= 44.16;

number_of_corners(15)=2;
corner_radius(26)=36.23;
cornering_dist(26)= 50.86;
corner_radius(27)=66.51;
cornering_dist(27)= 59.21;

number_of_corners(16)=1;
corner_radius(28)=88;
cornering_dist(28)= 46.08;

number_of_corners(17)=2;
corner_radius(29)=70.29;
cornering_dist(29)= 86.48;
corner_radius(30)=25.66;
cornering_dist(30)= 54.41;

number_of_corners(18)=2;
corner_radius(31)=35.94;
cornering_dist(31)= 56.46;
corner_radius(32)=32.91;
cornering_dist(32)= 51.7;

%Saves all of the corner radii into one array.
for i=1:32
    corner_total(i)=corner_radius(i);
end





  








%Car Parameters
vehicle_weight = 530;
c_d = 1;
frontal_a = 9.24021;
rho = 0.00236;
f_0 = 0.01417;
f_s = 0.01;
drivetrain_eff = .85;
g = 32.17405;
pi = 3.1416;
r=10.25/12;
iteration_time=0.001;
wheelbase=61;
track_width=50;
track_a=track_width/2;

%coeff_friction=1.5;
hcg=12;
weight_distribution=.5;  % Percent to the Front
c_length=wheelbase*weight_distribution;
b_length=wheelbase*(1-weight_distribution);


%Transmission and shifting parameters

nt_1 = 29/12;   %1st Gear Ratio
nt_2 = 26/15;   %2nd Gear Ratio
nt_3 = 21/16;   %3rd Gear Ratio
nt_4 = 21/20;   %4th Gear Ratio
nt_5 = 21/25;   %5th Gear Ratio
nf = 40/14;     %Final Gear Ratio
np = 62/22;     %Primary Reduction Ratio
shift_time = 0.01;   %Time required to shift
shift_rpm = 9500;   %Shift RPM
shifts=0;           
total_shifts=0; %Keeps track of number of shifts





%Tire 1=Goodyear 13"
%Tire 2=Hoosier 13"
%Tire 3=Hoosier Small 13"
%Tire 4=Michelin 13"
tire_choice=4;



%Competition Points Breakdown - FSAE 2011
%Not Currently Used.
%IN PROGRESS

competition_acceleration=75;
competition_skidpad=50;
competition_autocross=150;
competition_fueleconomy=100;
competition_endurance=300;

competition_static_total=325;






















%Vehicle Velocity Calculations
    velocity_iteration(1)=0;
    rpm_iteration(1)=0;


    for iii=1:90
    
        rpm_iteration(iii+1)=rpm_iteration(iii)+10;
        velocity_iteration(iii+1)=(((rpm_iteration(iii)*2*pi)/60)/(nf*np*nt_1))*r;
    
    end
    
    fact_1=polyfit(rpm_iteration,velocity_iteration,1);
    veloc_fact_1=1/fact_1(1);
    
    clear velocity_iteration
    clear rpm_iteration
    velocity_iteration(1)=0;
    rpm_iteration(1)=0;
    
    
    for iii=1:90
    
        rpm_iteration(iii+1)=rpm_iteration(iii)+10;
        velocity_iteration(iii+1)=(((rpm_iteration(iii)*2*pi)/60)/(nf*np*nt_2))*r;
    
    end
    
    fact_2=polyfit(rpm_iteration,velocity_iteration,1);
    veloc_fact_2=1/fact_2(1);
    
    clear velocity_iteration
    clear rpm_iteration
    velocity_iteration(1)=0;
    rpm_iteration(1)=0;
    
    for iii=1:90
    
        rpm_iteration(iii+1)=rpm_iteration(iii)+10;
        velocity_iteration(iii+1)=(((rpm_iteration(iii)*2*pi)/60)/(nf*np*nt_3))*r;
    
    end
    
    fact_3=polyfit(rpm_iteration,velocity_iteration,1);
    veloc_fact_3=1/fact_3(1);
    
    clear velocity_iteration
    clear rpm_iteration
    velocity_iteration(1)=0;
    rpm_iteration(1)=0;

    for iii=1:90
    
        rpm_iteration(iii+1)=rpm_iteration(iii)+10;
        velocity_iteration(iii+1)=(((rpm_iteration(iii)*2*pi)/60)/(nf*np*nt_4))*r;
    
    end
    
    fact_4=polyfit(rpm_iteration,velocity_iteration,1);
    veloc_fact_4=1/fact_4(1);
    
    clear velocity_iteration
    clear rpm_iteration
    velocity_iteration(1)=0;
    rpm_iteration(1)=0;

    for iii=1:90
    
        rpm_iteration(iii+1)=rpm_iteration(iii)+10;
        velocity_iteration(iii+1)=(((rpm_iteration(iii)*2*pi)/60)/(nf*np*nt_5))*r;
    
    end
    
    fact_5=polyfit(rpm_iteration,velocity_iteration,1);
    veloc_fact_5=1/fact_5(1);
    
    clear velocity_iteration
    clear rpm_iteration
    velocity_iteration(1)=0;
    rpm_iteration(1)=0;




%For 0.5 seconds assume maximum acceleraton that the tires are able to
%produce with 50% static weight distribution.  Then you begin weight
%transfer.

%%IN PROGRESS

v_initial = 0;     %Starting Velocity
time_launching=0;    %Starting Time
t=iteration_time;
launching_time=0; %Time for maximum tractive acceleration.
launch_distance=0;
ntf=nt_1*np*nf;
t_ime=0;


%Maximum tire traction for 0.5 seconds without weight transfer.  After 0.5
%seconds, there is weight transfer until 15mph(22 ft/s).  The car launches from a stop
%and accelerates until it reaches 15mph.  At that point, the software will
%jump into the regular sequence and calculate accordingly.  The following
%while loop will record the distance traveled and time taken to get to
%15mph.

while (launching_time<0.5)
weight_on_tire=vehicle_weight/2;
friction_instant = TireFunction(weight_on_tire,tire_choice);
maximum_tractive_force=friction_instant*(weight_on_tire)*(b_length/wheelbase)/(1-(hcg/wheelbase)*friction_instant);
a_x_max=maximum_tractive_force/(vehicle_weight/g);

launch_distance = launch_distance+v_initial*t+0.5*a_x_max*t^2;
v_initial=v_initial+a_x_max*t;
launching_time=t+launching_time;

end

rpm_cur = veloc_fact_1*v_initial;

while (v_initial<22)

for i=1:126
    
    if rpm(i)>rpm_cur
        L_tor = (tor(i-1)+(rpm_cur-rpm(i-1))*(tor(i)-tor(i-1))/(rpm(i)-rpm(i-1)));
        break
        
    
    end
end    

r_x = (f_0+3.25*f_s*((v_initial)/100)^(2.5))*vehicle_weight;
d_a = 0.5*rho*((v_initial)^2)*c_d*frontal_a;
   
    
a_x=((((L_tor*ntf*drivetrain_eff)/r)-r_x-d_a))/(vehicle_weight/g);    

%Weight Distribution & Instantaneous Friction between Tires and Ground
%Calculation
weight_on_tire=vehicle_weight*(b_length/wheelbase+(a_x/g)*(hcg/wheelbase));
friction_instant = TireFunction(weight_on_tire,tire_choice);    
    
%Maximum Tractive Force Calculation
maximum_tractive_force=friction_instant*(weight_on_tire)*(b_length/wheelbase)/(1-(hcg/wheelbase)*friction_instant);
a_x_max=maximum_tractive_force/(vehicle_weight/g);

    if (a_x>a_x_max)
        
        a_x=a_x_max;
    end


launch_distance = launch_distance+v_initial*t+0.5*a_x*t^2;
v_initial=v_initial+a_x*t;
launching_time=t+launching_time;

rpm_cur = veloc_fact_1*v_initial;



end


%%IN PROGRESS
v_i = v_initial;  




%Cornering Simulation

%IN PROGRESS
%a_x_cornering is given an initial assumption and later iterrated.
%The higher the acceleration, the higher the weight transfer to the outside
%wheels.  
%Tire 1 is the left front tire.
%Tire 2 is the right front tire.
%Tire 3 is the left rear tire.
%Tire 4 is the right rear tire.

a_x_cornering=1*g;
i=0;
difference=100;
tire_1=100;
tire_3=100;


while (difference>.1) || (i<10)

cornering_force=0;
old_cornering=a_x_cornering;

    
cornering_weight_distribution=vehicle_weight*(track_a/track_width+(a_x_cornering/g)*(hcg/track_width));

    if (cornering_weight_distribution) >= (vehicle_weight)
        
        cornering_weight_distribution=vehicle_weight;
        
    end

tire_1=(vehicle_weight-cornering_weight_distribution)/2;
tire_3=tire_1;
tire_2=(cornering_weight_distribution)/2;
tire_4=tire_2;

weight_on_tire=tire_1;
friction_instant = TireFunction(weight_on_tire,tire_choice);
cornering_force=friction_instant*tire_1+cornering_force;

weight_on_tire=tire_2;
friction_instant = TireFunction(weight_on_tire,tire_choice);
cornering_force=friction_instant*tire_2+cornering_force;


weight_on_tire=tire_3;
friction_instant = TireFunction(weight_on_tire,tire_choice);
cornering_force=friction_instant*tire_3+cornering_force;

weight_on_tire=tire_4;
friction_instant = TireFunction(weight_on_tire,tire_choice);
cornering_force=friction_instant*tire_4+cornering_force;


a_x_cornering=(cornering_force/vehicle_weight)*g;

difference=(abs(old_cornering-a_x_cornering)/((old_cornering+a_x_cornering)*(1/2)))*100;
i=i+1;


end
clear i
clear difference
clear old_cornering



% weight_on_tire is used as a variable to calculate weight transfer in
% corner due to use later on in forward acceleration.  This variable is not
% to be confused with the actual fact that it is specific to the
% invidual tire.

%Cornering Parameters
cornering_accel=a_x_cornering;  
cornering_size=size(corner_total);
cornering_time=0;



for i=1:cornering_size(2)
    
    v_corner(i)=sqrt(cornering_accel*corner_total(i));
    cornering_time=cornering_dist(i)/v_corner(i)+cornering_time;
    
end
    
    
    
    
    
    
    
    
    
    
%Braking Parameters
braking_accel=(-1.5)*g;









for o=1:straight_sections
    
%Determines the brake-to velocity based on the first corner after straight section    
if o==1    
v_cornering=v_corner(1);
else
    v_cornering=v_corner(1+number_of_corners(o-1));
end

    

    
veloc_fact = veloc_fact_1;


rpm_cur = veloc_fact_1*v_i;      
ntf=nt_1*np*nf;
 

if (rpm_cur > shift_rpm )
    veloc_fact=veloc_fact_2;
    rpm_cur = veloc_fact_2*v_i;
    ntf=nt_2*np*nf;
end

if (rpm_cur > shift_rpm)
    veloc_fact = veloc_fact_3;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_3*np*nf;
end

if (rpm_cur > shift_rpm)
    veloc_fact = veloc_fact_4;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_4*np*nf;
end

if (rpm_cur > shift_rpm)
    veloc_fact = veloc_fact_5;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_5*np*nf;
end





for i=1:126
    
    if rpm(i)>rpm_cur
        L_tor = (tor(i-1)+(rpm_cur-rpm(i-1))*(tor(i)-tor(i-1))/(rpm(i)-rpm(i-1)));
        break
        
    
    end
end



x=0;
t=iteration_time;
count=0;



    while (x<dist(o))
    count=count+1;
    
        
    r_x = (f_0+3.25*f_s*((v_i)/100)^(2.5))*vehicle_weight;
    d_a = 0.5*rho*((v_i)^2)*c_d*frontal_a;
   
    
    a_x=((((L_tor*ntf*drivetrain_eff)/r)-r_x-d_a))/(vehicle_weight/g);
    
    %Weight Distribution & Instantaneous Friction between Tires and Ground
    %Calculation
    weight_on_tire=vehicle_weight*(b_length/wheelbase+(a_x/g)*(hcg/wheelbase));
    friction_instant = TireFunction(weight_on_tire,tire_choice);
    
    %Maximum Tractive Force Calculation
    maximum_tractive_force=friction_instant*(weight_on_tire)*(b_length/wheelbase)/(1-(hcg/wheelbase)*friction_instant);
    a_x_max=maximum_tractive_force/(vehicle_weight/g);
    
    
    
    if (a_x>a_x_max)
        
        a_x=a_x_max;
    end
    
    x = x+v_i*t+0.5*a_x*t^2;
    v_i=v_i+a_x*t;
    v_store(count) = v_i;   %Stores values of velocity
    x_store(count) = x;
    a_x_store(count) = a_x;
    
    rpm_cur = veloc_fact*v_i;
    
        if (rpm_cur > shift_rpm ) && (veloc_fact>veloc_fact_2)
            veloc_fact = veloc_fact_2;
            rpm_cur = veloc_fact*v_i;
            t_ime=t_ime+shift_time;
            shifts=shifts+1;
            ntf=nt_2*np*nf;
        end

        if (rpm_cur > shift_rpm ) && (veloc_fact>veloc_fact_3)
           
            veloc_fact = veloc_fact_3;
            rpm_cur = veloc_fact*v_i;
            t_ime=t_ime+shift_time;
            shifts=shifts+1;
            ntf=nt_3*np*nf;
        end

        if (rpm_cur > shift_rpm ) && (veloc_fact>veloc_fact_4)
            veloc_fact = veloc_fact_4;
            rpm_cur = veloc_fact*v_i;
            t_ime=t_ime+shift_time;
            shifts=shifts+1;
            ntf=nt_4*np*nf;
        end

        if (rpm_cur > shift_rpm ) && (veloc_fact>veloc_fact_5)
            veloc_fact = veloc_fact_5;
            rpm_cur = veloc_fact*v_i;
            t_ime=t_ime+shift_time;
            shifts=shifts+1;
            ntf=nt_5*np*nf;
        end


        for j=1:126
    
         if rpm(j)>rpm_cur
             L_tor = (tor(j-1)+(rpm_cur-rpm(j-1))*(tor(j)-tor(j-1))/(rpm(j)-rpm(j-1)));
             break
             
         end
         
         
        end
        
        
        
    end
    
r2_value = 0;

    
for i=1:7
    
    if (r2_value<0.99)   
        
        p1_accel = polyfit(x_store,v_store,i);
        f_accel = polyval(p1_accel,x_store);
        mean_veloc = mean(v_store);
        j_accel = sum((f_accel-v_store).^2);
        s_accel = sum((v_store-mean_veloc).^2);
        r2_value = 1-j_accel/s_accel;
    end

    
end

decel_accel_point_1=0;
accelerating=0;
decelerating=0;
total_accel_decel=1; 



    while (total_accel_decel>0.1)
        decel_accel_point_1=decel_accel_point_1+0.01;
        accelerating = polyval(p1_accel,decel_accel_point_1);
        decelerating = sqrt((v_cornering)^2-2*braking_accel*(dist(o)-decel_accel_point_1));
        total_accel_decel = abs(accelerating-decelerating);
    end
    
% End of the calculation for braking point of First Straight.
% End of the calculation for braking point of First Straight.


%Recalculation of Acceleration and Deceleration

v_i = v_initial;  


x=0;
t=iteration_time;
t_ime=0;
count=0;  
shifts=0;




clear v_store
clear x_store
clear a_x_store








veloc_fact = veloc_fact_1;

rpm_cur = veloc_fact*v_i;      
ntf=nt_1*np*nf;
     
if (rpm_cur > shift_rpm )
    veloc_fact = veloc_fact_2;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_2*np*nf;
end

if (rpm_cur > shift_rpm)
    veloc_fact = veloc_fact_3;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_3*np*nf;
end

if (rpm_cur > shift_rpm)
    veloc_fact = veloc_fact_4;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_4*np*nf;
end

if (rpm_cur > shift_rpm)
    veloc_fact = veloc_fact_5;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_5*np*nf;
end

for i=1:126
    
    if rpm(i)>rpm_cur
        L_tor = (tor(i-1)+(rpm_cur-rpm(i-1))*(tor(i)-tor(i-1))/(rpm(i)-rpm(i-1)));
        break
        
    
    end
end





















while (x<decel_accel_point_1)
    count=count+1;
    
        
    r_x = (f_0+3.25*f_s*((v_i)/100)^(2.5))*vehicle_weight;
    d_a = 0.5*rho*((v_i)^2)*c_d*frontal_a;
   
    
    a_x=((((L_tor*ntf*drivetrain_eff)/r)-r_x-d_a))/(vehicle_weight/g);
    
    %Weight Distribution & Instantaneous Friction between Tires and Ground
    %Calculation
    weight_on_tire=vehicle_weight*(b_length/wheelbase+(a_x/g)*(hcg/wheelbase));
    friction_instant = TireFunction(weight_on_tire,tire_choice);
    
    %Maximum Tractive Force Calculation
    maximum_tractive_force=friction_instant*(weight_on_tire)*(b_length/wheelbase)/(1-(hcg/wheelbase)*friction_instant);
    a_x_max=maximum_tractive_force/(vehicle_weight/g);
    
    if (a_x>a_x_max)
        
        a_x=a_x_max;
    end
    
    horsepower_store(count) = L_tor*rpm_cur/5252;
    x = x+v_i*t+0.5*a_x*t^2;
    v_i=v_i+a_x*t;
    a_x_store(count) = a_x;
    
    rpm_cur = veloc_fact*v_i;
    
        if (rpm_cur > shift_rpm ) && (veloc_fact>veloc_fact_2)
            veloc_fact = veloc_fact_2;
            rpm_cur = veloc_fact*v_i;
            t_ime=t_ime+shift_time;
            shifts=shifts+1;
            ntf=nt_2*np*nf;
        end

        if (rpm_cur > shift_rpm ) && (veloc_fact>veloc_fact_3)
           
            veloc_fact = veloc_fact_3;
            rpm_cur = veloc_fact*v_i;
            t_ime=t_ime+shift_time;
            shifts=shifts+1;
            ntf=nt_3*np*nf;
        end

        if (rpm_cur > shift_rpm ) && (veloc_fact>veloc_fact_4)
            veloc_fact = veloc_fact_4;
            rpm_cur = veloc_fact*v_i;
            t_ime=t_ime+shift_time;
            shifts=shifts+1;
            ntf=nt_4*np*nf;
        end

        if (rpm_cur > shift_rpm ) && (veloc_fact>veloc_fact_5)
            veloc_fact = veloc_fact_5;
            rpm_cur = veloc_fact*v_i;
            t_ime=t_ime+shift_time;
            shifts=shifts+1;
            ntf=nt_5*np*nf;
        end

  
    
    
    
    t_ime=t+t_ime;
    
    
        for j=1:126
    
         if rpm(j)>rpm_cur
             L_tor = (tor(j-1)+(rpm_cur-rpm(j-1))*(tor(j)-tor(j-1))/(rpm(j)-rpm(j-1)));
             break
             
         end
         
         
        end
        
        
        
end



avg_accel(o)= mean(a_x_store);
avg_horsepower(o)=mean(horsepower_store);

clear horsepower_store
clear a_x_store
clear v_store
clear x_store

%Time to brake before the first corner

deceleration_velocity_1 = v_i;
braking_dist(o) = dist(o)-decel_accel_point_1;
t_ime = t_ime + 2*braking_dist(o)/(deceleration_velocity_1+v_cornering);
shifts=shifts+1;
t_ime=t_ime+shift_time;


%Total Time Calculation
total_straights_time(o)=t_ime;


%Total Shifts Calculation
total_shifts(o) = shifts;




%Set-up for next Straight Calculation

v_i = v_initial;  


x=0;
t=iteration_time;
t_ime=0;
count=0;  
shifts=0;




veloc_fact = veloc_fact_1;

rpm_cur = veloc_fact*v_i;      
ntf=nt_1*np*nf;
     
if (rpm_cur > shift_rpm )
    veloc_fact = veloc_fact_2;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_2*np*nf;
end

if (rpm_cur > shift_rpm)
    veloc_fact = veloc_fact_3;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_3*np*nf;
end

if (rpm_cur > shift_rpm)
    veloc_fact = veloc_fact_4;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_4*np*nf;
end

if (rpm_cur > shift_rpm)
    veloc_fact = veloc_fact_5;
    rpm_cur = veloc_fact*v_i;
    ntf=nt_5*np*nf;
end

for i=1:126
    
    if rpm(i)>rpm_cur
        L_tor = (tor(i-1)+(rpm_cur-rpm(i-1))*(tor(i)-tor(i-1))/(rpm(i)-rpm(i-1)));
        break
        
    
    end
end






end


    
    
    
%Total Distance Traveled Calculation
total_dist=(sum(dist)+sum(cornering_dist))*number_of_laps;

    
    
%Total Time 
total_time=(sum(total_straights_time)+cornering_time)*number_of_laps;


    
%Fuel Consumption Calculation

fuel_consumption_rate=BSFC*mean(avg_horsepower);


total_fuel=fuel_consumption_rate*total_time;
%The weight-density of gasoline is 6.073lb/US Gal.
fuel_used=total_fuel/6.073;





    
    
    
disp('THE STARTING VELOCITY IN FT/S IS:')    
disp(v_i)
disp('THE STARTING VELOCITY IN MPH IS:')
v_mph=v_i/1.46667;
disp(v_mph)
disp('THE TOTAL DISTANCE IS:')
disp(total_dist)
disp('THE TOTAL TIME IS:')
disp(total_time)
disp('The number of shifts is:')
disp((sum(total_shifts))*number_of_laps);
disp('THE AVERAGE LONGITUDINAL ACCELERATION IS:')
disp(mean(avg_accel));
disp('THE MEAN HORSEPOWER IS:')
disp(mean(avg_horsepower));
disp('THE TOTAL FUEL USED IS:')
disp(fuel_used);
disp('THE CORNERING TIME IS:')
disp((cornering_time)*number_of_laps);
disp('NUMBER OF CORNERS:')
disp((number_of_corners));

disp('Shifts per run')
disp(total_shifts)





disp('-------------------------------')
disp('-------------------------------')
disp('-------------------------------')
disp('-------------------------------')
disp('-------------------------------')
disp('-------------------------------')







end

