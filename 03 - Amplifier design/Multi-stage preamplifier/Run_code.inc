diagram shown: 1,2

; New diagram definitions
; 1: "Pulse Time Domain (Steady State)" 
; 2: "ASE Spectrum in Steady State" ;    (taking the previously calculated steady state)
; 2: "Wavelength Domain, Final State" 
; 3: "Powers vs. gain fiber position (1st Stage)"
; 4: "Powers vs. gain fiber position (2nd Stage)"
; 5: "Powers vs. gain fiber position (3rd Stage)"
; 6: "Output Pulse Duration per compressor Length" 
; 7: "Output Pulse Energy vs Gain Fiber Length variation"

; 11: Dispersion parameter (5 dB/m)
; 12: Dispersion parameter (22 dB/m)

; 17: "Yb absortion and emission at 1000 nm"
; 17: "Er absortion and emission at 1500 nm"
; 18: "Er absortion and emission at 980 nm"
; 19: "Er absortion and emission at 1500 nm"
; 20: "Er absortion and emission at 980 nm"


number_of_stages := 3
Amp := 0; {0: Erbium, 1: Ytterbium, 2: Measurements comparison}

Pulse_dBm := 0 {If Pulse_dBm is in 1, the output pulse plot will be done in dBm}
print_outputs:= 12; {1: prints noise figure information, 2: prints pulse evolution}


ASE_factor := 0.5	{A PM Isolator is implemented which filters half of the unpolarized ASE}

;If 22 dB fiber is being used, data needs to be rescalated.
gain_fiber_1 := 24
gain_fiber_2 := 27
gain_fiber_3 := 30

SplitterFactor := 7.19


if Amp = 0 
then 
begin ; {Erbium doped amplifier} 
include "Units.inc" 
;--------------- Stage 1 ---------------
;Parameters of the first stage


{mean laser wavelength}
lMean := 1542 nm
Pulse_gain_saturation := 0 {Read the code in case this wants to be changed}

Loss_S := 2.5;
Loss_P := 0.47;

Pump_Power1 := 50.0 {100}
Pump_Power2 := 0

L_pm:= 2.35                    
L_pm2:= .01                  
L_comp:= .1e-2
L_f:= 1.50			{Length of the gain fiber in stage 1}       
GainFiber1 := 5;	{22 Erbium can be 5,22,80} 
dir_p1 := forward;	{Direction of pumping of LD1}
dir_p2 := forward;	{Direction of pumping of LD2}

;--------------- Stage 2 ---------------
;Parameters of the second stage
L_2pm:= .1                   
L_2pm2:= .1                   
L_2comp:= .07
; Parameters of amplifier fiber:
L_2f:= .5			{0.5 Length of the gain fiber in stage 2}                
GainFiber2 := 22
dir_p3 := backward;	{Direction of pumping of LD1}
dir_p4 := backward;	{Direction of pumping of LD2}
Pump_power3 := 100.0
Pump_power4 := 0

;--------------- Stage 3 ---------------
;Parameters of the third stage
L_3pm:= 1.52           
L_3pm2:= .52                  
L_3comp:= .18
; Parameters of amplifier fiber:
L_3f:= .75			{Length of the gain fiber in stage 3}        
GainFiber3 := 22	{Erbium can be 5,22,80} 
dir_p5 := backward;	{Direction of pumping of LD1}
dir_p6 := backward;	{Direction of pumping of LD1}
Pump_power5 := 300.0
Pump_power6 := 0

;--------------- Input pulse --------------- 

Initial_comb_power := 0.08 {0.2[mW] Average comb power} 
cw_power := Initial_comb_power;			{Comb power - considered in cw regime for running initial simulations} 
cw_laser := 0;				{if cw:= 1, the cw laser is turned on. For a different kind of simulation, needs to be tested}

Pulse_duration  := 120;      {in femtoseconds}
Repetition_rate := 12.126   {in GHz} 
Chirp := 0 					{in GHz/ps}

if number_of_stages = 1 then include "FirstStage_Er_or_Yb_v2.inc"
if number_of_stages = 2 then include "SecondStage_Er_or_Yb_v2.inc"
if number_of_stages = 3 then include "ThirdStage_Er_or_Yb_v2.inc"


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

(*

for indice:=1 to 1 do 
begin
for indice2:=1 to 4 do 
begin
Pump_Power1 := Pump_Power_input[indice];
cw_power := Seed_input[indice2]*SplitterFactor;
show "Input pump: ", Pump_Power1:d3:"mW";
show "cw power: ", cw_power:d3:"mW";
include "FirstStage_Er_or_Yb_v2.inc";
show "-------------------------";
end;
end;
*)



if Amp = 1 
then 
begin {Ytterbium doped amplifier}
include "Units.inc" 

lMean := 1050 nm

show "Use Peter Mayers YFA"
;Parameters of the first stage

Loss_S := 0.2;
Loss_P := 0.147;

;--------------- Stage 1 ---------------
Pump_Power1 := 10.0
Pump_Power2 := 0

L_pm:= 2.35                    
L_pm2:= .1                  
L_comp:= .1
L_f:= .5      
;Ytterbium : 1 
GainFiber1 := 1 
dir_p1 := forward;
dir_p2 := forward;

;--------------- Stage 2 ---------------
;Parameters of the second stage
L_2pm:= .52                   
L_2pm2:= .52                   
L_2comp:= .1e-2
; Parameters of amplifier fiber:
L_2f:= 1.5               
GainFiber2 := 1
dir_p3 := backward;
dir_p4 := backward;
Pump_power3 := 100.0
Pump_power4 := 0

;--------------- Stage 3 ---------------
;Parameters of the third stage
L_3pm:= .52           
L_3pm2:= .52                  
L_3comp:= .01
; Parameters of amplifier fiber:
L_3f:= 1.75       
GainFiber3 := 1
dir_p5 := backward;
dir_p6 := backward;
Pump_power5 := 200.0
Pump_power6 := 0

;--------------- Input pulse --------------- 

Initial_comb_power := 0.2;	{Average comb power} 
cw_power := Initial_comb_power;			{Comb power - considered in cw regime for running initial simulations} 
cw_laser := 0;				{if cw:= 1, the cw laser is turned on. For a different kind of simulation}

Pulse_duration  := 95;      {in femtoseconds}
Repetition_rate := 12.126   {in GHz} 
Chirp := 0 					{in GHz/ps}



if number_of_stages = 1 then include "FirstStage_Er_or_Yb_v2.inc"
if number_of_stages = 2 then include "SecondStage_Er_or_Yb_v2.inc"
if number_of_stages = 3 then include "ThirdStage_Er_or_Yb_v2.inc"


if print_outputs = 1 then 
begin
show " "
show "-.. Amplifier parameters ..-"
show "  Length (1st):   ", L_f:d2:"m"
show "  Pumping (1st):  ", Pump_Power1:d3:"mW"
show "  Length (2nd):   ", L_2f:d2:"m"
show "  Pumping (2nd):  ", Pump_power3:d3:"mW"
end

show "Ytterbium needs to get the correct" 
show "dispersion parameters"
show "Ytterbium needs to get the correct" 
show "gain fiber"

end


if Amp = 2 
then 
begin ; {For testing} 
include "Units.inc"                 {General file needed for Unit conversion}
;--------------- Stage 1 ---------------
;Parameters of the stage to study

Loss_S := 2;
Loss_P := .47;

Pump_Power1 := 18
Pump_Power2 := 0

L_pm:= 2.35                    
L_pm2:= .1                  
L_comp:= .1
L_f:= 1.67			{Length of the gain fiber in stage 1}       
GainFiber1 := 5 ;	{Erbium can be 5,22,80} 
dir_p1 := forward;	{Direction of pumping of LD1}
dir_p2 := forward;	{Direction of pumping of LD2}

;--------------- Input pulse --------------- 

Initial_comb_power := 0.005*7.591040462427746 {Average comb power} 
cw_power := Initial_comb_power;			{Comb power - considered in cw regime for running initial simulations} 
cw_laser := 0;				{if cw:= 1, the cw laser is turned on. For a different kind of simulation}

Pulse_duration  := 50;      {in femtoseconds}
Repetition_rate := 0.100   {12.126 in GHz} 
Chirp := 0 					{in GHz/ps}



include "FirstStage_Er_or_Yb_v2.inc"



end
