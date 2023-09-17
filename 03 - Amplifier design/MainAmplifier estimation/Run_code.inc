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


number_of_stages := 2
Amp := 0; {0: Erbium, 1: Ytterbium, 2: Measurements comparison}

Pulse_dBm := 1 {If Pulse_dBm is in 1, the output pulse plot will be done in dBm}
print_outputs:= 12; {1: prints noise figure information, 2: prints pulse evolution}

ASE_factor := 0.5	{No PM Isolator is implemented which filters half of the unpolarized ASE}

;If 22 dB fiber is being used, data needs to be rescalated.
gain_fiber_1 := 30
gain_fiber_2 := 75


if Amp = 0 
then 
begin ; {Erbium doped amplifier} 
include "Units.inc" 
;--------------- Stage 1 ---------------
;Parameters of the first stage


{mean laser wavelength}
lMean := 1542 nm
Pulse_gain_saturation := 0 {Read the code in case this wants to be changed}

Loss_S := 1.5;
Loss_P := 0.;

Pump_Power1 := 956.0 
Pump_Power2 := 0

L_pm:= 2.35                    
L_pm2:= .01                  
L_comp:= .1e-2
L_f:= 3.70			{Length of the gain fiber in stage 1}       
GainFiber1 := 22;	{Erbium can be 5,22,80} 
dir_p1 := forward;	{Direction of pumping of LD1}
dir_p2 := forward;	{Direction of pumping of LD2}

;--------------- Stage 2 ---------------
;Parameters of the second stage
L_2pm:= .1                   
L_2pm2:= .1                   
L_2comp:= .18
; Parameters of amplifier fiber:
L_2f:= 1.9			{Length of the gain fiber in stage 2}                
GainFiber2 := 22
dir_p3 := backward;	{Direction of pumping of LD1}
dir_p4 := backward;	{Direction of pumping of LD2}
Pump_power3 := 18000.0
Pump_power4 := 0


Initial_comb_power := 10 {[mW] Average comb power} 
cw_power := Initial_comb_power;			{Comb power - considered in cw regime for running initial simulations} 
cw_laser := 0;				{if cw:= 1, the cw laser is turned on. For a different kind of simulation, needs to be tested}

Pulse_duration  := 120;      {in femtoseconds}
Repetition_rate := 12.126  {in GHz} 
Chirp := 0 					{in GHz/ps}

include "MainAmplifier.inc"

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
