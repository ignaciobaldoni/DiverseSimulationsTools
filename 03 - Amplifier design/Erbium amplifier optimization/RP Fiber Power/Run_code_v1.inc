{*
 Here, we set the parameters of the amplifiers which will be implemented in a secondary script where the simulation is run. 
*}

{Script needed for implementing different units inside the main script. Default from RP Fiber Power}
include "Units.inc" 

{Set the graphs you are interested to see according to the plotcode value given below}
diagram shown: 

; Diagram definitions
; 1: "Pulse Time Domain (Steady State)" 
; 2: "ASE Spectrum in Steady State" ;    (taking the previously calculated steady state)
; 2: "Wavelength Domain, Final State" 
; 3: "Powers vs. gain fiber position (1st Stage)"
; 4: "Powers vs. gain fiber position (2nd Stage)"
; 5: "Powers vs. gain fiber position (3rd Stage)"
; 6: "Output Pulse Duration per compressor Length" 
; 7: "Pulse duration vs devices"
; 11: "Dispersion parameter (5 dB/m)"
; 12: "Dispersion parameter (22 dB/m)"
; 17: "Er absortion and emission at 1500 nm"
; 18: "Er absortion and emission at 980 nm"
; 19: "Er absortion and emission at 1500 nm"
; 20: "Er absortion and emission at 980 nm"


;-------- Input pulse properties -------

Initial_comb_power := 0.105 {Average comb power in mW} 
cw_laser := 0;				{if cw:= 1, the cw laser is turned on. Not for this simulation.}

Pulse_duration  := 125;     {in femtoseconds}
Repetition_rate := 12.126   {in GHz} 
Chirp := 0 					{in GHz/ps}


;-------- Amplifier configuration parameters -------
number_of_stages := 2

; Stage 1
actFiber1_min_value := 0.2;{1.5}
actFiber1_max_value := 0.2;{1.5}
Step1 := 0.1

; Stage 2
actFiber2_min_value := 0.1
actFiber2_max_value := 0.1
Step2 := 0.1

ASE_factor := 0.5	
{0.5 A PM Isolator is implemented which filters half of the unpolarized ASE}

lMean := 1542 nm {Center wavelength}

{Outputs of interest}
Pulse_dBm := 0 		{If = 1: Output pulse is plotted in dBm}
print_outputs:= 2; 	{1: prints noise figure information, 2: prints pulse evolution}
print_outputs := 3	{3: Pulse energy, Peak power, effective gain and average output power }
No_z_steps:=100     {Number of steps along the fibers }


{Script initiation}
{*
For some applications, we need to set a length value of 0. 
However, the program does not allow a lower value than 0.1e-2. 
So, we go with that logic when we need to consider a zero value.
*}


begin
for XX:= actFiber1_min_value to actFiber1_max_value step Step1 do
begin
for XX2:= actFiber2_min_value to actFiber2_max_value step Step2 do
begin

;--------------- Stage 1 ---------------
; Pumping
Loss_P := 0.47; 			{Pump losses estimation from splice and measurements}
Pump_Power1 := 250.0 		{Pump power in mW}
l_p1:=976 nm                {Wavelength pump 1}
dir_p1 := backward;			{Direction of pumping}

{When studying pulse duration, the fiber length of the amplifier components is important. Here the measured values in our setup. Everything in meters}
isolator := 0.6 ; 
wdm := 0.10; 
fc_apc := 0.2; 
lensed_fibers := 2 *0;
initial_value := 2.35; {Microcomb output lensed fibers + 2 x DWDM for pump cancellation} 
initial_value2 := 0.1e-2
length_Test := 0.1e-2;

;fiber parameters 
L_pm:= initial_value + 0.80; {Fiber length before active fiber for Stage 1}             
L_pm2:= 0.35;  				 {Fiber length after active fiber for Stage 1}                   
L_f:= XX			{Length of the gain fiber in stage 1}       
GainFiber1 := 22;	{Erbium can be 5,22,80} 
gain_fiber_1 := 24  {If 22 dB fiber is being used, data can be rescalated for slight different absorption values}

if GainFiber1 = 5 then Loss_S := 0.45;
if GainFiber1 = 22 then Loss_S := 2;

;--------------- Stage 2 ---------------
;Parameters of the second stage
; Pumping
Loss_p2 := 0.47; 		{Pump losses estimation from splice and measurements}
Pump_power3 := 250.0 	{Pump power in mW}
l_p3:=976 nm            {Wavelength pump 3 (only pump in second stage)}
dir_p3 := backward;		{Direction of pumping of LD1}

;fiber parameters 
L_2pm:= 0.10;       {Fiber length before active fiber for Stage 2}             
L_2pm2:= 0.67;      {Fiber length after active fiber for Stage 2}                 
L_2f:= XX2			{Length of the gain fiber in stage 2}                
GainFiber2 := 22 
gain_fiber_2 := 27 	{If 22 dB fiber is being used, data can be rescalated for slight different absorption values}

if GainFiber2 = 5 then Loss_s2 := 0.5;
if GainFiber2 = 22 then Loss_s2 := 1.72;

;------------- Compressor --------------
L_comp:= 0.1e-2;
L_2comp:= 0.1e-2; 
loss_comp := 2.5  {parasitic losses on dispersion compesanting fiber [dB/m]}
 
show L_f:d4 
show L_2f:d4


;------------ Run the script -----------

if number_of_stages = 1 then include "FirstStage_Er.inc"
if number_of_stages = 2 then include "SecondStage_Er.inc"
if number_of_stages = 3 then include "ThirdStage_Er.inc"

{*
; For implementing in optimization script in python, we export the data obtained in the parameters sweep.

if number_of_stages = 2 then
write [L_f:d4, ", ",L_2f:d4, ", ",
 Dif:d4, ", ",
 Out_avg:d4, ", ",
 tau_2ndStage:d4],
 >>"5dB_22dB_v3.dat" ;fileName

;if number_of_stages = 1 then
write [L_f:d4, ", ",L_2f:d4, ", ",
 Dif_1:d4, ", ",
 Out_avg_1:d4, ", ",
 tau_1stStage:d4],
 >>"5dB_forward_low_loss.dat"  ;fileName
*}

 
end
if print_outputs = 1 then 
begin
show " "
show "-.. Amplifier parameters ..-"
show "  Length (1st):   ", L_f:d2:"m"
show "  Pumping (1st):  ", Pump_Power1:d3:"mW"
show "  Length (2nd):   ", L_2f:d2:"m"
show "  Pumping (2nd):  ", Pump_power3:d3:"mW"
end



end
end
