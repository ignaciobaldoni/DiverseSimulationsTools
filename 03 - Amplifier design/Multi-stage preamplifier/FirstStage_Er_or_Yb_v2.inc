(*  Simulation of Er fiber amplifier for fs pulses
    Edited December 2022 by Ignacio Baldoni(based on previous files and examples)
    Uses RP Fiber power by R. Paschotta
    Additional input files required in the same folder: ...
*)

show "-------  First stage -------"

real_fiber22dB := gain_fiber_1

if GainFiber1 = 5 then include "Er5dB-RealData-iXBlue_V0.inc"
if GainFiber1 = 22 then include "Er22dB-RealData-iXBlue_V0.inc"
if GainFiber1 = 80 then include "Er80dB-RealData-iXBlue_V0.inc"

{lMean := 1542 nm} 

if GainFiber1 = 1 then include "Yb-MultipleStages.inc";
if GainFiber1 = 1 then N_Er := 10e24;
if GainFiber1 = 1 then r_core:=0.5*4.3 um;
{if GainFiber1 = 1 then lMean := 1050 nm} {mean laser wavelength}



A_eff_pm:=86 um^2               { effective mode area }
w_s_pm:=sqrt(A_eff_pm/pi)       { effective mode radius with 5um }
A_eff_comp:=20 um^2               { effective mode area }
w_s_comp:=sqrt(A_eff_comp/pi)       { effective mode radius with 5um }
I_s_pm(r):=exp(-2*(r/w_s_pm)^2) { Intensity profile in the PM fiber}
I_s_comp(r):=exp(-2*(r/w_s_comp)^2) { Intensity profile in the Dispersion Compensanting fiber}
l_s_pm:=lMean                   { Wavelength is set to mean }
loss_pm:= 0.02                  { Loss, probably in dB/m }
GVD_pm:=-(lMean^2 / (2pi*c)) * 17.5 ps/(nm km)  { GVD of the PM1550 fiber is assumed to -25509 fs^2/m}
;show "GVD pm fiber", GVD_pm:d3:"s/m"

D_compressor := -100 ps/(nm km) { Thorlabs PMDCF (and apparently standard values)}
GVD_compressor := -(lMean^2 / (2pi*c)) * D_compressor  

w_s:= 0.5 * 10.1 um; {check whether for Ytterbium is the same mode field diameter}


{GVD parameter is beta2}

; Parameters of amplifier fiber:

No_z_steps:=100              { no steps along the fiber }
(*
show "L_f:    ", L_f:d3:"m"
show "L_pm1:  ", L_pm:d3:"m"
show "L_pm2:  ", L_pm2:d3:"m"
show "L_comp: ", L_comp:d3:"m"
*)

;show "GVD_pm:  ", GVD_pm/fs^2:f0:"fs^2/m"
n2_pm:=2.2e-20                  { nonlinear refractive index n2 for this fiber, unit not given, is this used at all? Appears indicative}
;gamma_pm:=(2*pi/lMean) * n2_pm / A_eff_pm   { Gamma for this fiber, not used in any calculation - Appears indicative}
;show "gamma_pm: ", gamma_pm:d3:"rad/(W m)"


loss_comp := 0.50  { parasitic losses on dispersion compesanting fiber in dB/m - possibly overestimating}

loss_s := Loss_S    {1.2}
loss_p := Loss_P;   {0.05}
;show "loss signal: ", loss_s:d3:"dB"

; Parameters of the channels:
l_p1:=976 nm                 {Wavelength pump 1}
w_p1:=w_p                    {Mode Radius pump 1}
I_p1(r):=exp(-2*(r/w_p1)^2)  {Intensity profile pump 1}
l_p2:=976 nm                 {Wavelength Pump 2}
w_p2:=w_p                    {Mode Radius pump 2}
I_p2(r):=exp(-2*(r/w_p2)^2)  {Intensity profile pump 2}
l_p3:=976 nm                 {Wavelength Pump 3}
l_s1:=lMean                  {Signal wavelength}
w_s1:= w_s               	 {3.25 um Signal mode radius}

I_s1(r):=exp(-2*(r/w_s1)^2)  {Signal mode profile}
w_ASE:=w_s1                   {ASE mode radius}
I_ASE(r):=exp(-2*(r/w_ASE)^2) {ASE mode profile}

; Pump laser intensities
P_pump1_bw_in:= Pump_Power1 * mW       {Define here the pump power and direction of Pump 1}
P_pump2_bw_in:= Pump_Power2 * mW        {Define here the pump power and direction of Pump 3}
P_signal_fw_in:= cw_power * mW      {.08 mW Signal injection power}


; Array of ASE channels:
l1_ASE:= lMean - 42 nm               {ASE lower wavelength}
l2_ASE:= lMean + 58 nm               {ASE upper wavelength}
dl_ASE:=1 nm                   {ASE wavelength steps}
NoModes_ASE:=2                {Two ASE modes: forward and backward}
defarray c_ASE_fw[l1_ASE, l2_ASE, dl_ASE]   {ASE forward array definition}
defarray c_ASE_bw[l1_ASE, l2_ASE, dl_ASE]   {ASE backward array definition}
P_ASE_fw(x):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_ASE_fw[l],x))  {ASE forward power at pos x integrated over wavelength}
P_ASE_bw(x):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_ASE_bw[l],x))  {ASE backward power at pos x integrated over wavelength}

; Raman response function:
h_R(t):=(tau1^2+tau2^2)/(tau1*tau2^2)
 *exp(-t/tau2)*sin(t/tau1)
tau1:=12.2e-15
tau2:=32e-15
f_R:=0.18

; Define the fiber models:
def_model():=
  begin
  global allow all;
  { Fiber 1: passive PM fiber }
  set_device(1);
  set_fiber(L_pm, No_z_steps, '-');
  signal_pm:=addinputchannel(0, l_s_pm, 'I_s_pm', loss_pm, forward);
  finish_fiber();
  { Fiber 2: amplifier fiber }
  set_device(2);
  if GainFiber1 > 2 then set_fiber(L_f, No_z_steps, 'Er'); {Length of the gain fiber segment}
  if GainFiber1 <= 2 then set_fiber(L_f, No_z_steps, 'Yb'); {Length of the gain fiber segment}
  add_ring(r_core, N_Er);           {Fill core with Erbium from fiber include file}
  max_N:=N_Er;
  max_r:=r_core;
  P_in_max:=0;
  P_in_max:=maxr(P_in_max, P_pump1_bw_in);
  pump1_bw:=addinputchannel(P_pump1_bw_in, l_p1, 'I_p1', loss_p, dir_p1);
  P_in_max:=maxr(P_in_max, P_pump2_bw_in);
  pump2_bw:=addinputchannel(P_pump2_bw_in, l_p3, 'I_p2', loss_p, dir_p2);
  P_in_max:=maxr(P_in_max, P_signal_fw_in);
  signal_fw:=addinputchannel(P_signal_fw_in, l_s1, 'I_s1', loss_s, forward);
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     c_ASE_fw[l]:=addASEchannel(l, dl_ASE, NoModes_ASE, 'I_ASE', loss_s, forward);
     c_ASE_bw[l]:=addASEchannel(l, dl_ASE, NoModes_ASE, 'I_ASE', loss_s, backward);
     end;
  finish_fiber();
  { Fiber 3: passive PM fiber }
  set_device(3);
  set_fiber(L_pm2, No_z_steps, '-');
  signal_pm2:=addinputchannel(0, l_s_pm, 'I_s_pm', loss_pm, forward);
  finish_fiber();
  { Fiber 4: Dispersion compensating fiber }
  set_device(4);
  set_fiber(L_comp, No_z_steps, '-');
  signal_comp:=addinputchannel(0, l_s_pm, 'I_s_comp', loss_comp, forward);
  finish_fiber();
  end;
calc def_model()

; This derives the maximum ASE channel
P_ASE_max():=
  begin
  global l1_ASE, l2_ASE, dl_ASE, c_ASE_fw[], c_ASE_bw[];
  var P_max;
  P_max:=0;
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     P_max:=maxr(P_max, P_out(c_ASE_fw[l]));
     P_max:=maxr(P_max, P_out(c_ASE_bw[l]));
     end;
  P_max;
  end;


; These two functions are defined for the maximum ASE channel calculation
P_ASE_fw_out(l):=
  if l1_ASE <= l <= l2_ASE
  then P_out(c_ASE_fw[l])
P_ASE_bw_out(l):=
  if l1_ASE <= l <= l2_ASE
  then P_out(c_ASE_bw[l])
(* Results printing
show "----------------------------------------"
show "Initial powers:"
show "Average Power:      ", P_signal_fw_in:d3:"W"
show "Pump 1 (backward):  ", P_pump1_bw_in:d3:"W"
show "Pump 2 (backward):  ", P_pump2_bw_in:d3:"W"
show "----------------------------------------"
show "Output powers:"
*)

; {Calculate cw amplification} 
calc set_device(2)

{Results analysis}
P_pump1_bw_out:=P_out(pump1_bw)
;show "pump1 (forward):   ", P_pump1_bw_out:d3:"W"
;show "pump1 (at 1550 nm):   ", 0.2*P_pump1_bw_out:d3:"W"
P_pump2_bw_out:=P_out(pump2_bw)
;show "pump3 (backward):   ", P_pump2_bw_out:d3:"W"
P_signal_fw_out:=P_out(signal_fw)

;--------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------
diagram 3: 
"Powers vs. Position"
"First Stage"
x: 0, L_f
"position in fiber [m]", @x
y: 0, P_out(signal_fw)/uW*(1.3)
"Signal Power [uW]",@y
y2: 0, P_pump1_bw_in*(1.3)/mW
"Pump Power [mW]",@y2
frame
legpos 500, 650
hx
hy
defarray P_fiber[0,L_f,0.01]
defarray Pump_fiber[0,L_f,0.01]
defarray ASE_fw[0,L_f,0.01]
defarray ASE_bw[0,L_f,0.01]
defarray dop_conc[0,L_f,0.01]


P_signal_fiber():= 
  begin
  global allow all;
  for z:=0 to z:= L_f step 0.01
  do begin
     P_fiber[z]:=P(signal_fw,z);
	 Pump_fiber[z]:=P(pump1_bw,z);
	 ASE_fw[z]:=sum(l:=l1_ASE to l2_ASE step dl_ASE, ASE_factor*P(c_ASE_fw[l],z));
     ASE_bw[z]:=sum(l:=l1_ASE to l2_ASE step dl_ASE, ASE_factor*P(c_ASE_bw[l],z));
     dop_conc[z]:=n(z,2);
     end;
  end
calc P_signal_fiber();
P_fiber_z(x):=P_fiber~~[x]
Pump_fiber_z(x):=Pump_fiber~~[x] 
ASE_fw_z(x):=ASE_fw~~[x] 
ASE_bw_z(x):=ASE_bw~~[x] 
dop_conc_z(x):=dop_conc~~[x] 

f: P_fiber_z(x)/uW,
  color=green, width=3, maxconnect=1, style=solid, "Signal amplification (1st Stage)"
f: Pump_fiber_z(x)/mW,yscale=2,
  color=teal, width=3, maxconnect=1, style=fdashed, "Pump power (1st Stage)"
;f: ASE_fw_z(x)/uW,yscale=1,color=olive, width=3, maxconnect=1, style=solid, "Fwd ASE (1st Stage)"
;f: ASE_bw_z(x)/uW,yscale=1,color=olive, width=3, maxconnect=1, style=fdashed, "Bwd ASE (1st Stage)"
;f: dop_conc_z(x)*1,yscale=2,color=gray, width=3, maxconnect=1, style=fdashed, "n2 concentration"
  
;--------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------

show "signal (forward):    ", P_signal_fw_out:d3:"W"

G_signal1:=sp_gain(signal_fw)



P_ASE_fw:=sum(l:=l1_ASE to l2_ASE step dl_ASE, ASE_factor*P_out(c_ASE_fw[l]))
P_ASE_bw:=sum(l:=l1_ASE to l2_ASE step dl_ASE, ASE_factor*P_out(c_ASE_bw[l]))
calc (G_ASE:=0; for l:=l1_ASE to l2_ASE step dl_ASE do G_ASE:=maxr(G_ASE, sp_gain(c_ASE_fw[l])))

;show "Gain ASE:            ", G_ASE:d3:np:"dB"
;show "P_ASE_fw:            ", P_ASE_fw:d3:"W"
;show "P_ASE_bw:            ", P_ASE_bw:d3:"W"
;show "Linear NF:              ", NF(signal_fw):d3
;show "NF (ASE):            ", 10*lg(NF(c_ASE_fw[1530 nm])):d3:"dB"

if print_outputs = 1 or print_outputs = 12 then 
begin 
show "NF (1st):            ", 10*lg(NF(signal_fw)):d3:"dB"
SNR_out1 := 10*lg(P_signal_fw_out/P_ASE_fw)
show "SNR_out(1st):        ", SNR_out1:d3:"dB"
show "G_signal1:           ", G_signal1:d3:np:"dB"
show "ASE1:                ", 100*P_ASE_fw/(P_ASE_fw+P_signal_fw_out):d3:"%"
end



; ---------------------------------------------------------------------
; ---------------------------------------------------------------------
diagram 2:

"ASE spectrum and comb optical spectrum"

x: l1_ASE/nm, l2_ASE/nm
"wavelength [nm]", @x
y: -60,20 {-60, 20}
"PSD of ASE [dBm/nm]", @y
frame
hx
hy

dBm(P):=10*lg(P/1e-3)

; First store the wavelength-dependent ASE powers in arrays,
; so that we can later get interpolated values:
defarray P_ASE_fw[l1_ASE,l2_ASE,dl_ASE]
defarray P_ASE_bw[l1_ASE,l2_ASE,dl_ASE]
StoreASE():=
  begin
  global allow all;
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     P_ASE_fw[l]:=ASE_factor*P_out(c_ASE_fw[l]);
     P_ASE_bw[l]:=ASE_factor*P_out(c_ASE_bw[l]);
     end;
  end
calc StoreASE()
P_ASE_fw_l(l):=P_ASE_fw~~[l] 
;show "ASE @", lMean-12," nm: ",dBm(P_ASE_fw_l(lMean - 12 nm)*(nm/dl_ASE)):d3:"dB"
P_ASE_bw_l(l):=P_ASE_bw~~[l]

f: dBm(P_ASE_fw_l(x*nm)), {olive}
  color=green, width=3, maxconnect=1, style=fdashed, "forward ASE (1st Stage)"
;f: dBm(P_ASE_bw_l(x*nm)),  color=olive, width=3, maxconnect=1, style=fdashed, "backward ASE (1st Stage)"

FName0$:="C:\Users\Administrator.MENLOPC188\Documents\Simulations Ignacio\ASE1_seed"
FName1$:="mW_Pump"
FName2$:="mW.dat"


FName$:=(FName0$)+str(cw_power:d3)+(FName1$)+str(Pump_Power1)+(FName2$)


{ For saving the ASE in a file}
write [x:d5, ", ", (dBm(P_ASE_fw_l(x*nm))):d5], for x:= l1_ASE/nm to l2_ASE/nm step 1,
 >(FName$)
 





ASE_1530 :=  dBm(P_ASE_fw_l(lMean-12 nm));
;show "ASE @ ",lMean-12nm, ": ", ASE_1530:d3:"dB"
ASE_1000 :=  dBm(P_ASE_fw_l(lMean));



;show "---------------------------------"
;show "-- 1st Stage PULSE SIMULATIONS --"
;show "---------------------------------"
;if GainFiber1>2 then show "Gain Fiber (1st):     ", GainFiber1:d3:"dB/m"
if GainFiber1<2 then show "Gain Fiber: Ytterbium"
; ---------------------------------------------------------------------
; ---------------------------------------------------------------------

;-------------------------------------
;gridsize determination for the laser pulse
lgrid:=lMean {covers XXX nm to XXX nm if T_range and pulsegridsize are unchanged}
if Pulse_duration < 51 then T_range:=6 ps  {Enough for pulse duration description}
if Pulse_duration >= 51 then T_range:= 6 ps  {Enough for pulse duration description}
pulsegridsize:=2^10  {NOT clear what this is, probably clarified in RP Fiber Power Manual}

; Ultrashort pulses:
f_rep:= Repetition_rate * GHz       { <- set here Signal in repetition rate}
P_av:= Initial_comb_power *mW    {P_signal_fw_in <- set here Signal in power}
tau0:= Pulse_duration * fs            {<- set here Signal pulse duration at injection}
chirp0:= Chirp * GHz/ps;              {100 GHz/ps<- Assumed NO chirp, pulse is Gaussian. To adapt this needs UNIT or description?}
E0:=P_av/f_rep           {Pulse energy varies determination}
l0_p:=lMean              {Center wavelength of the pulse is set to mean wavelenth of it}

; ALTERNATIVE: recall output pulse from measured data file ("Osz_Output_Pulse.pdt")
;FName$:="Osz_Output_Pulse.pdt"
;calc load_pulse(FName$);
;calc set_pulse_size_f(pulsegridsize) {otherwise gridsize from load_pulse is taken!!!!}

defarray PSD_in_pulse[l1_ASE,l2_ASE,dl_ASE]



Initial_pulse():=
  begin
  global allow all;
  set_pulse_grid(T_range,pulsegridsize,lgrid);
  startpulse_s(E0,tau0,chirp0);
  for l:= l1_ASE to l2_ASE step dl_ASE
	do begin
	PSD_in_pulse[l]:=P_l(l)*Repetition_rate;{1e1/2}
	end;
	end
calc Initial_pulse();


;------------------------------------------------------------------------

DoPulsePropagation_PM_fiber() :=
  begin  {Pulse propagation in PM fiber}
  global allow all;
  set_pulse_grid(T_range,pulsegridsize,lgrid);
  startpulse_s(E0,tau0,chirp0);           {initial pulse as sech2 pulse }
  store_pulse(0);
  { passive PM fiber: }
  set_device(1);
  set_GVD(GVD_pm);
  set_n2(n2_pm);
  pp_fiber(1, signal_pm);{Do the pulse propagation}
  store_pulse(1);{save the pulse}
  end

calc DoPulsePropagation_PM_fiber()

; Where does it say that the input for the AMP is the output of the PM Fiber 1?

if cw_laser = 0 then calc (set_device(2); set_P_in(signal_fw, 0)) {signal_fw}

DoPulsePropagation_amplifier():=
  begin  {Pulse propagation in amplifier fiber: }
  global allow all;
  set_device(2);{Device 2 is the gain fiber!}
  set_ss(1);{Self-steppening is considered!}
  set_dnr("h_R(t)",f_R); {Raman is considered!}
  set_GVD_l("GVD_f(l)")   ;{Use instead the GVD from fiber.inc file}
  set_n2(n2_gain {XXX});  {set n2 of the Gain fiber as calculated above}
  pp_fiber(2, signal_fw); {Propagates the pulse through the fiber 2, using the optical channel "signal_fw"}
  store_pulse(2);         {save the pulse}
  if maxr(P_t(-0.95*T_range/2), P_t(+0.95*T_range/2)) > 0.01 * P_p()
  then error("Time range too small!?!");
  end

;calc DoPulsePropagation_amplifier(){no steady state!}
     { switch off the cw input - we now look at pulses }


DoPulsePropagation_postamp():=
  begin  {Pulse propagation in amplifier fiber: }
  global allow all;
  set_device(3);
  set_GVD(GVD_pm);
  set_n2(n2_pm);
  startpulse_recall(2);
  pp_fiber(3, signal_pm2);{Do the pulse propagation}
  store_pulse(3);{save the pulse}
  end

{if cw_laser = 0 then calc (set_device(3);set_P_in(signal_pm2,0))}

{GVD_compressor := -17}

DoCompressor_postamp():=
  begin  {Pulse propagation in amplifier fiber: }
  global allow all;
  set_device(4);
  set_GVD(GVD_compressor);
  set_n2(n2_pm);
  startpulse_recall(3);
  pp_fiber(4, signal_comp);{Do the pulse propagation}
  store_pulse(4);{save the pulse}
  end

{if cw_laser = 0 then calc (set_device(4);set_P_in(signal_comp,0))}

DoPumping(T):= {simulates the influence of the pump wave only for a given time}
  begin
  global T;
  calc_dyn(0, T, minr(10 us, T));
  end

f_rep_ss:= 100 kHz;   { Then, influence of the pump is only for 10 µs. Don't know where does this value come from }

 

Find_SteadyState():=
  { do multiple amplification and pumping cycles until the steady state is reached;
    return the number of iteration steps }
  begin
  global f_rep, E0, tau0, chirp0, f_rep_ss, Pulse_gain_saturation;
  var E_out, E_out_l, N;
  
  set_gain_sat(f_rep/f_rep_ss); { effectively increased gain saturation per pulse }
  {set_gain_sat(Pulse_gain_saturation)  effectively increased gain saturation per pulse }
  E_out:=0;
  N:=0;
  repeat
    recall_pulse(1);
    inc(N); {N = N+1}
    E_out_l:=E_out;
    DoPulsePropagation_amplifier();
    { store_pulse(10+N); {}
    DoPumping(1/f_rep_ss);
    E_out:=E_p();
  until abs(E_out-E_out_l) < 1e-3 * E_out;
  N; 
  DoPulsePropagation_postamp();
  DoCompressor_postamp();
  save_pulse_t("Pulse_from_first_stage.txt")
  end
(*
show "Initial pulse energy: ", E0:d3:"J"
show "Initial T pulse :     ", tau0:d3:"s"
show "Initial Peak Power:   ", 0.88*E0/tau0:d3:"W"
show "Initial bandwidth:    ", 0.314/tau0:d3:"Hz"
show "Initial chirp:        ", (chirp0*1e-12):d3:"Hz/ps"
*)

(*
show "----- Pulse evolution: -----"
*)
calc Find_SteadyState()

;show "N:       ", Find_SteadyState():d3
;show "E_sat:   ", E_sat(signal_pm2):d3:"J"
; Display some parameters of the obtained pulse:
;show "Initial Avg Power: " P_av:d3:"W"
tau_1stStage := tau_p()
if print_outputs = 2 or print_outputs = 12 then 
begin
show "Energy per pulse:    ", E_p():d3:"J"
show "Peak power:          ", P_p():d3:"W" {Peak power}
show "T pulse:             ", tau_p():d3:"s"
show "Geff:                ", 10*lg(E_p()/E0):d3:"dB"
show "Output avg power:    ", E_p()*f_rep:d3:"W"
end

(*
show "Bandwidth:           ", dl_p():d3:"(n)m"
show "TBP:                 ", TBP():d3
*)
Geff1:= 10*lg(E_p()/E0);
;show "Noise:              ", 10*lg(NF(signal_pm2)):d3:"dB"
;show "f_rep_ss:           ", f_rep_ss:d3:"Hz"



;if GainFiber1 > 2 then show "Difference with ASE: ",(dBm(P_l(lMean-12nm)*Repetition_rate)-ASE_1530):d3:"dBm";
;if GainFiber1 <= 2 then show "Difference with ASE: ",(dBm(P_l(lMean)*Repetition_rate)-ASE_1000):d3:"dBm";



;-------------------------------------------------------------------------------
;---------- From here on, the description of the different diagrams of interest
;---------------------------------------------------------------------

diagram 1:

"Pulse Time Domain (Steady State)"

time_in := -2 ps;
time_end := - time_in

x: time_in, time_end
"time [s]", @x

frame
hx
hy


; First store the wavelength-dependent powers in arrays,
; so that we can later get interpolated values:

dtime := 10 fs; 
defarray PulseTime[time_in,time_end,dtime];

Pulse1Stage():=
  begin
  global allow all;
  for t:=time_in to time_end step dtime
  do begin
     PulseTime[t]:=P_t(t);
     end;
  end
calc Pulse1Stage()
PulseTime_t(t):=PulseTime~~[t] 

(*
write [t:d5, ", ", PulseTime_t(t):d5],
 for t:=time_in to time_end step dtime,
 >"C:\Users\Administrator.MENLOPC188\Documents\Simulations Ignacio\Pulse1stStage.dat"
*)
if Pulse_dBm=0 then
begin
"Power (2nd & 3rd Stage) [W]", @y
"Power (1st Stage) [mW]", @y2
y: 0,P_t(0)*1.1; 
y2: 0,PulseTime_t(0)/mW*1.1;
f: PulseTime_t(x)/mW, yscale=2, {red}
  color=green, width=3, maxconnect=1, style=solid, "1st Stage"
end
  
if Pulse_dBm = 1 then
begin
y: -70,50;
"Power [dBm]", @y
f: dBm(PulseTime_t(x)), {red}
  color=green, width=3, maxconnect=1, style=solid, "1st Stage"
end



; ---------------------------------------------------------------------
; ---------------------------------------------------------------------

{dBm(P):=10*lg(P/1 mW)}

diagram 2:

x: l1_ASE/nm, l2_ASE/nm
;"wavelength [nm]", @x
y: dBm(P_l(l1_ASE))-20,dBm(P_l(lMean))+20;{-60,20{0,.0003 1.01*P_l(1530 nm)/mW}
;"power spectral density [pJ/nm]", @y
frame
hx
hy

; First store the wavelength-dependent powers in arrays,
; so that we can later get interpolated values:
defarray PSD[l1_ASE,l2_ASE,dl_ASE]
StorePSD():=
  begin
  global allow all;
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     PSD[l]:=P_l(l)*Repetition_rate;{1e1/2}
     end;
  end
calc StorePSD()
PSD_l(l):=PSD~~[l]
PSD_in_pulse_l(l):=PSD_in_pulse~~[l] 

f: dBm(PSD_in_pulse_l(x*nm)), 
  color=magenta, width=2, maxconnect=1, style=solid, "Initial pulse"
{ For saving the initial pulse in a file}
FName3$:="C:\Users\Administrator.MENLOPC188\Documents\Simulations Ignacio\spectrum_in_seed"
FName4$:=(FName3$)+str(cw_power:d3)+(FName1$)+str(Pump_Power1)+(FName2$)
write [x:d5, ", ", (dBm(PSD_in_pulse_l(x*nm))):d5], for x:= l1_ASE/nm to l2_ASE/nm step 1,
 >(FName4$)
 
FName5$:="C:\Users\Administrator.MENLOPC188\Documents\Simulations Ignacio\spectrum_Out_seed"
FName6$:=(FName5$)+str(cw_power:d3)+(FName1$)+str(Pump_Power1)+(FName2$)
f: dBm(PSD_l(x*nm)), 
  color=green, width=3, maxconnect=1, style=solid, "1st Stage"
write [x:d5, ", ", (dBm(PSD_l(x*nm))):d5], for x:= l1_ASE/nm to l2_ASE/nm step 1,
 >(FName6$)

(* To consider the third stage 
f: dBm(P_l(x*nm)),
  color=green, width=3, maxconnect=1, style=solid, "2nd Stage"
 *)
; ---------------------------------------------------------------------
; ---------------------------------------------------------------------

