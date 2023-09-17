include "PreAmplifier.inc"

show "------- Second stage -------"


tau_1stStage:= tau_p()
PeakPower_1stStage := P_p()
Avg_power_1st := E_p()*f_rep;

P_signal_fw_in := P_signal_fw_out;{Avg_power_1st}

real_fiber22dB := gain_fiber_2

if GainFiber2 = 22 then include "Er22dB-RealData-iXBlue_V0.inc"

{lMean := 1542 nm;}


if GainFiber2 = 1 then include "Yb-MultipleStages.inc";
if GainFiber2 = 1 then N_Er := 10e24;
if GainFiber2 = 1 then r_core:=0.5*4.3 um;
{if GainFiber2 = 1 then lMean := 1050 nm   }                 {mean laser wavelength}

;show "lMean :", lMean




loss_pm:= 0.02                  { Loss, probably in dB/m }
GVD_pm:=-(lMean^2 / (2pi*c)) * 17.5 ps/(nm km)  { GVD of the PM1550 fiber is assumed to -25509 fs^2/m}
D_compressor := -100 ps/(nm km) { Thorlabs PMDCF (and apparently standard values)}
GVD_compressor := -(lMean^2 / (2pi*c)) * D_compressor

No_z_steps:=100              { no steps along the fiber }



; Pump laser intensities
P_pump3_bw_in:= Pump_power3 * mW       {Define here the pump power and direction of Pump 1}
P_pump4_bw_in:= Pump_power4 * mW        {Define here the pump power and direction of Pump 3}


r_core := 12 um * 0.5
w_p := r_core

; Parameters of the channels:
l_p3:=976 nm                 {Wavelength pump 1}
w_p3:=w_p                    {Mode Radius pump 1}
I_p3(r):=exp(-2*(r/w_p3)^2)  {Intensity profile pump 1}
I_p3(r):=(r <= w_p3)  		 {Pump intensity profile }
l_p4:=976 nm                 {Wavelength Pump 3}
w_p4:=w_p                    {Mode Radius pump 3}
I_p4(r):=exp(-2*(r/w_p4)^2)  {Intensity profile pump 3}
I_p4(r):=(r <= w_p4)  		 {Pump intensity profile }
l_s1:=lMean                  {Signal wavelength}
w_s1:= w_s               	 {3.25 um Signal mode radius}

I_s1(r):=exp(-2*(r/w_s1)^2)  {Signal mode profile}
w_ASE:=w_s1                   {ASE mode radius}
I_ASE(r):=exp(-2*(r/w_ASE)^2) {ASE mode profile}


{Fiber 2 (a double-clad fiber)}



loss_s := Loss_S;
loss_p := Loss_P;

;show "loss signal: ", loss_s:d3:"dB"

; Array of ASE channels:
l1_ASE:= lMean - 42 nm               {ASE lower wavelength}
l2_ASE:= lMean + 58 nm               {ASE upper wavelength}
dl_ASE:=1 nm                  		 {ASE wavelength steps}
NoModes_ASE2:=2*3 					 {NoModes_ASE:=2  ; Two ASE modes: forward and backward}
defarray c_2ASE_fw[l1_ASE, l2_ASE, dl_ASE]   {ASE forward array definition}
defarray c_2ASE_bw[l1_ASE, l2_ASE, dl_ASE]   {ASE backward array definition}
P2_ASE_fw(x):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_2ASE_fw[l],x))  {ASE forward power at pos x integrated over wavelength}
P2_ASE_bw(x):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_2ASE_bw[l],x))  {ASE backward power at pos x integrated over wavelength}






; ---------


; Define the fiber models:
def_model():=
  begin
  global allow all;
  { Fiber 6: amplifier fiber }
  set_device(6);
  if GainFiber2 > 2 then set_fiber(L_2f, No_z_steps, 'Er'); {Length of the gain fiber segment}
  if GainFiber2 <= 2 then set_fiber(L_2f, No_z_steps, 'Yb'); {Length of the gain fiber segment}
  add_ring(r_core, N_Er);           {Fill core with Erbium or Yb from fiber include file}
  max_N:=N_Er;
  max_r:=r_core;
  P_in_max:=0;
  P_in_max:=maxr(P_in_max, P_pump3_bw_in);
  pump3_bw:=addinputchannel(P_pump3_bw_in, l_p3, 'I_p3', loss_p, dir_p3);
  P_in_max:=maxr(P_in_max, P_pump4_bw_in);
  pump4_bw:=addinputchannel(P_pump4_bw_in, l_p4, 'I_p4', loss_p, dir_p4);
  P_in_max:=maxr(P_in_max, P_signal_fw_in);
  signal2_fw:=addinputchannel(P_signal_fw_in, l_s1, 'I_s1', loss_s, forward);
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     c_2ASE_fw[l]:=addASEchannel(l, dl_ASE, NoModes_ASE2, 'I_ASE', loss_s, forward);
     c_2ASE_bw[l]:=addASEchannel(l, dl_ASE, NoModes_ASE2, 'I_ASE', loss_s, backward);
     end;
  finish_fiber();
  end;
calc def_model()


; These two functions are defined for the maximum ASE channel calculation
P2_ASE_fw_out(l):=
  if l1_ASE <= l <= l2_ASE
  then P_out(c_2ASE_fw[l])+P_ASE_fw_l(l)
P2_ASE_bw_out(l):=
  if l1_ASE <= l <= l2_ASE
  then P_out(c_2ASE_bw[l])+P_ASE_bw_l(l)


; Outputs from cw case initial simulation
(*
show "----------------------------------------"
show "Initial powers:"
show "Average Power:      ", P_signal_fw_in:d3:"W"
show "Pump 1 (backward):  ", P_pump3_bw_in:d3:"W"
show "Pump 2 (backward):  ", P_pump4_bw_in:d3:"W"
show "----------------------------------------"
*)
;show "Output powers:"
calc set_device(6)

;--------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------
diagram 4:
"Powers vs. Position"
"Second Stage"

x: 0, L_2f
"position in fiber [m]", @x
y: 0, P_out(signal2_fw)*1.3/mW
"Signal Power [mW]",@y
y2: 0, P_pump3_bw_in*(1.3)/mW
"Pump Power [mW]",@y2
frame
legpos 200, 250
hx
hy

defarray P2_fiber[0,L_2f,0.01]
defarray Pump2_fiber[0,L_2f,0.01]
P2_signal_fiber():= 
  begin
  global allow all;
  for z:=0 to z:= L_2f step 0.01
  do begin
     P2_fiber[z]:=P(signal2_fw,z);
	 Pump2_fiber[z]:=P(pump3_bw,z);
     end;
  end
calc P2_signal_fiber();
P2_fiber_z(x):=P2_fiber~~[x]
Pump2_fiber_z(x):=Pump2_fiber~~[x]
f: P2_fiber_z(x)/mW,{yscale=2,}
  color=red, width=3, maxconnect=1, style=solid, "Signal amplification (2nd Stage)"  
f: Pump2_fiber_z(x)/mW,yscale=2,
  color=magenta, width=3, maxconnect=1, style=fdashed, "Pump Power (2nd Stage)" 

defarray gain_ase2[l1_ASE,l2_ASE,dl_ASE]
defarray gain_amplifier2[l1_ASE,l2_ASE,dl_ASE]

GainASE2():=
	begin
	global allow all;
	for l:= l1_ASE to l2_ASE step dl_ASE
	do begin
		gain_ase2[l]:=sp_gain(c_2ASE_fw[l]);
		gain_amplifier2[l]:= sp_gain(signal2_fw);
		end;
	end
calc GainASE2();


(* Plot previous version
x: L_f,L_2f+L_f
"position in fiber [m]", @x
y: 0, 1500e-6
frame
legpos 600, 150
hx
hy
defarray P2_fiber[0,L_2f,0.01]
P2_signal_fiber():= 
  begin
  global allow all;
  for z:=0 to z:= L_2f step 0.01
  do begin
     P2_fiber[z]:=P(signal2_fw,z);
     end;
  end
calc P2_signal_fiber();
P2_fiber_z(x):=P2_fiber~~[x-L_f] 
f: P2_fiber_z(x),{yscale=2,}
  color=red, width=3, maxconnect=1, style=fdashed, "Signal amplification (2nd Stage)"  
*)
;--------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------
;show "prueba", P_ASE_fw_l(1530 nm)*10^(gain_ase2[1530 nm]/10)
P_pump3_bw_out:=P_out(pump3_bw)
;show "pump3 (backward):   ", P_pump3_bw_out:d3:"W"
P_pump4_bw_out:=P_out(pump4_bw)
;show "pump4 (backward):   ", P_pump4_bw_out:d3:"W"
P_signal_fw_out:=P_out(signal2_fw)

show "signal (forward):    ", P_signal_fw_out:d3:"W"

G_signal2:=sp_gain(signal2_fw)
(*
show "G_signal2:          ", G_signal2:d3:np:"dB"
*)
P2_ASE_fw:=sum(l:=l1_ASE to l2_ASE step dl_ASE, ASE_factor*P_out(c_2ASE_fw[l])+P_ASE_fw_l(l)*10^(gain_ase2[l]/10))
P2_ASE_bw:=sum(l:=l1_ASE to l2_ASE step dl_ASE, ASE_factor*P_out(c_2ASE_bw[l])+P_ASE_bw_l(l)*10^(gain_ase2[l]/10))
calc (G_ASE:=0; for l:=l1_ASE to l2_ASE step dl_ASE do G_ASE:=maxr(G_ASE, sp_gain(c_2ASE_fw[l])))
;show "P2_ASE_fw:         ", P2_ASE_fw:d3:"W"
;show "P2_ASE_bw:         ", P2_ASE_bw:d3:"W"
if print_outputs = 1 or print_outputs = 12 then 
begin 
show "NF (2nd):            ", 10*lg(NF(signal2_fw)):d3:"dB"
SNR_out2:= 10*lg(P_signal_fw_out/P2_ASE_fw)
show "SNR_out(2nd):        ", SNR_out2:d3:"dB"
show "ASE2:                ", 100*P2_ASE_fw/(P2_ASE_fw+P_signal_fw_out):d3:"%"
end


;--------------------------------------------------------------------------------------
diagram 2:

dBm(P):=10*lg(P/1e-3)

; First store the wavelength-dependent ASE powers in arrays,
; so that we can later get interpolated values:
defarray P2_ASE_fw[l1_ASE,l2_ASE,dl_ASE]
defarray P2_ASE_bw[l1_ASE,l2_ASE,dl_ASE]
StoreASE2():=
  begin
  global allow all;
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     P2_ASE_fw[l]:=ASE_factor*P_out(c_2ASE_fw[l])+P_ASE_fw_l(l)*10^(gain_ase2[l]/10);
     P2_ASE_bw[l]:=ASE_factor*P_out(c_2ASE_bw[l])+P_ASE_bw_l(l)*10^(gain_ase2[l]/10);
     end;
  end
calc StoreASE2()
P2_ASE_fw_l(l):=P2_ASE_fw~~[l] 
;show "Forward ASE @ CenterWavelength nm: ",dBm(P2_ASE_fw_l(lMean)*(nm/dl_ASE)):d3:"dB"
P2_ASE_bw_l(l):=P2_ASE_bw~~[l]

f: dBm(P2_ASE_fw_l(x*nm)), {teal}
  color=red, width=3, maxconnect=1, style=fdashed, "forward ASE (2nd Stage)"
;f: dBm(P2_ASE_bw_l(x*nm)),  color=teal, width=3, maxconnect=1, style=fdashed, "backward ASE (2nd Stage)"

ASE2_1530 :=  dBm(P2_ASE_fw_l(lMean-12 nm));
;show "ASE @ ",lMean-12nm, ": ", ASE2_1530:d3:"dB"
ASE2_1000 :=  dBm(P2_ASE_fw_l(lMean));


; ---------------------------------------------------------
; ------ Dynamic simulations for pulse amplification ------
; ---------------------------------------------------------

;show "---------------------------------"
;show "-- 2nd Stage PULSE SIMULATIONS --"
;show "---------------------------------"
;show "Gain Fiber (2nd):     ", GainFiber2:d3:"dB/m"

;gridsize determination for the laser pulse
lgrid:=lMean {covers XXX nm to XXX nm if T_range and pulsegridsize are unchanged}
T_range:=5 ps  {Enough for pulse duration description}
pulsegridsize:=2^10  {NOT clear what this is, probably clarified in RP Fiber Power Manual}

; Ultrashort pulses:
f_rep:= Repetition_rate * GHz       { <- set here Signal in repetition rate}
P_av:= Avg_power_1st     {P_signal_fw_in <- set here Signal in power}
tau1:= tau_1stStage;           {<- set here Signal pulse duration at injection}
chirp0:= 0 GHz/ps;              {Unnecessary given it is taken according to last result}
E0:=P_av/f_rep           {Pulse energy varies determination}
l0_p:=lMean              {Center wavelength of the pulse is set to mean wavelenth of it}

; ALTERNATIVE: recall output pulse from measured data file ("Osz_Output_Pulse.pdt")
;FName$:="Osz_Output_Pulse.pdt"
;calc load_pulse(FName$);
;calc set_pulse_size_f(pulsegridsize) {otherwise gridsize from load_pulse is taken!!!!}


DoPulsePropagation2_PM_fiber() :=
  begin  {Pulse propagation in PM fiber}
  global allow all;
  { passive PM fiber: }
  load_pulse("Pulse_from_first_stage.txt");
  {recall_pulse(4);}
  set_device(1);
  set_L(L_2pm);
  set_GVD(GVD_pm);
  set_n2(n2_pm);
  pp_fiber(1, signal_pm);{Do the pulse propagation}
  store_pulse(6);{save the pulse}
  end

calc DoPulsePropagation2_PM_fiber()

; Where does it say that the input for the AMP is the output of the PM Fiber 1?

if cw_laser = 0 then calc (set_device(6); set_P_in(signal2_fw, 0)) 

DoPulsePropagation2_amplifier():=
  begin  {Pulse propagation in amplifier fiber: }
  global allow all;
  set_device(6);{Device 2 is the gain fiber!}
  set_ss(1);{Self-steppening is considered!}
  set_dnr("h_R(t)",f_R); {Raman is considered!}
  set_GVD_l("GVD_f(l)")   ;{Use instead the GVD from fiber.inc file}
  set_n2(n2_gain {XXX});  {set n2 of the Gain fiber as calculated above}
  pp_fiber(6, signal2_fw); {Propagates the pulse through the fiber 2, using the optical channel "signal2_fw"}
  store_pulse(7);         {save the pulse}
  if maxr(P_t(-0.95*T_range/2), P_t(+0.95*T_range/2)) > 0.01 * P_p()
  then error("Time range too small!?!");
  end

;calc DoPulsePropagation2_amplifier(){no steady state!}
     { switch off the cw input - we now look at pulses }


DoPulsePropagation2_postamp():=
  begin  {Pulse propagation in amplifier fiber: }
  global allow all;
  set_device(3);
  set_L(L_2pm2);
  set_GVD(GVD_pm);
  set_n2(n2_pm);
  startpulse_recall(7);
  pp_fiber(3, signal_pm2);{Do the pulse propagation}
  store_pulse(8);{save the pulse}
  end

{if cw_laser = 0 then calc (set_device(3);set_P_in(signal_pm2,0))}

{GVD_compressor := -17}

DoCompressor2_postamp(x):=
  begin  {Pulse propagation in amplifier fiber: }
  global allow all;
  set_device(4);
  set_L(x); {set_L(L_2comp);}
  set_GVD(GVD_compressor);
  set_n2(n2_pm);
  startpulse_recall(8);
  pp_fiber(4, signal_comp);{Do the pulse propagation}
  store_pulse(9);{save the pulse}
  end

{if cw_laser = 0 then calc (set_device(4);set_P_in(signal_comp,0))}


DoPumping2(T):= {simulates the influence of the pump wave only for a given time}
  begin
  global T;
  calc_dyn(0, T, minr(10 us, T));
  end

f_rep_ss:= 100 kHz;   { Then, influence of the pump is only for 10 µs. Don't know where does this value come from }

Find_SteadyState2(x):=
  { do multiple amplification and pumping cycles until the steady state is reached;
    return the number of iteration steps }
  begin
  global f_rep, E0, tau1, chirp0, f_rep_ss,Pulse_gain_saturation;
  var E_out, E_out_l, N;
  
  set_gain_sat(f_rep/f_rep_ss); { effectively increased gain saturation per pulse }
  {set_gain_sat(Pulse_gain_saturation)  effectively increased gain saturation per pulse }
  E_out:=0;
  N:=0;
  repeat
    recall_pulse(6);
    inc(N); {N = N+1}
    E_out_l:=E_out;
    DoPulsePropagation2_amplifier();
    { store_pulse(10+N); {}
    DoPumping2(1/f_rep_ss);
    E_out:=E_p();
  until abs(E_out-E_out_l) < 1e-3 * E_out;
  N; 
  DoPulsePropagation2_postamp();
  DoCompressor2_postamp(x);
  save_pulse_t("Pulse_from_second_stage.txt")
  end
(*
show "Initial pulse energy: ", E0:d3:"J"
show "Initial T pulse :     ", tau1:d3:"s"
show "Initial Peak Power:   ", PeakPower_1stStage:d3:"W"
show "Initial bandwidth:    ", 0.314/tau1:d3:"Hz"
;show "Initial chirp:        ", (chirp0*1e-12):d3:"Hz/ps"


show "----- Pulse evolution: -----"
*)
calc Find_SteadyState2(L_2comp)

; Display some parameters of the obtained pulse:
;show "Initial Avg Power: " P_av:d3:"W"
if print_outputs = 2 or print_outputs = 12 then 
begin
show "Energy per pulse:    ", E_p():d3:"J"
show "Peak power:          ", P_p():d3:"W" {Peak power}
show "T pulse:             ", tau_p():d3:"s"
show "Geff:                ", 10*lg(E_p()/E0):d3:"dB"
show "Output avg power:    ", E_p()*f_rep:d3:"W"
end

tau_2ndStage := tau_p()
;show "Bandwidth:           ", dl_p():d3:"(n)m"
;show "TBP:                 ", TBP():d3

Geff2:= 10*lg(E_p()/E0);
;show "Noise:              ", 10*lg(NF(signal2_fw)):d3:"dB"
;show "f_rep_ss:           ", f_rep_ss:d3:"Hz"

TotalGain:= Geff2+Geff1

(*
show "-         Results          -"
*)
show "Total Gain at 2nd:   ", TotalGain:d3:"dB"

if GainFiber2 > 2 then show "Difference with ASE: ",(dBm(P_l(lMean-12nm))-ASE2_1530):d3:"dBm";
if GainFiber2 <= 2 then show "Difference with ASE: ",(dBm(P_l(lMean))-ASE2_1000):d3:"dBm";

;-------------------------------------------------------------------------------
;---------- From here on, the description of the different diagrams of interest
;-------------------------------------------------------------------------------

diagram 1:



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
defarray PulseTime2[time_in,time_end,dtime];

Pulse2Stage():=
  begin
  global allow all;
  for t:=time_in to time_end step dtime
  do begin
     PulseTime2[t]:=P_t(t);
     end;
  end
calc Pulse2Stage()
PulseTime2_t(t):=PulseTime2~~[t] 
if Pulse_dBm = 0 then 
begin
y: 0,P_t(0)*1.1;

f: PulseTime2_t(x), yscale=1, {violet}
  color=red, width=3, maxconnect=1, style=solid, "2nd Stage"
if number_of_stages = 3 then f: P_t(x),yscale=1,"3rd Stage",color=blue, width=3
end

if Pulse_dBm = 1 then
begin 
f: dBm(PulseTime2_t(x)), yscale=1, {violet}
  color=red, width=3, maxconnect=1, style=solid, "2nd Stage"
if number_of_stages = 3 then f: dBm(P_t(x)),yscale=1,"3rd Stage",color=blue, width=3
end 

; ---------------------------------------------------------------------
; ---------------------------------------------------------------------

dBm(P):=10*lg(P/1 mW)

diagram 2:

x: l1_ASE/nm, l2_ASE/nm
;"wavelength [nm]", @x
y: -100,100;{dBm(P_l(l1_ASE))-10,dBm(P_l(lMean))+20 ... -60,20{0,.0003 1.01*P_l(lMean)/mW}
;"power spectral density [pJ/nm]", @y
frame
hx
hy

; First store the wavelength-dependent powers in arrays,
; so that we can later get interpolated values:
defarray PSD2[l1_ASE,l2_ASE,dl_ASE]
StorePSD2():=
  begin
  global allow all;
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     PSD2[l]:=P_l(l);
     end;
  end
calc StorePSD2()
PSD2_l(l):=PSD2~~[l] 

f: dBm(PSD2_l(x*nm)), 
  color=red, width=3, maxconnect=1, style=solid, "2nd Stage"
if number_of_stages = 3 then f: dBm(P_l(x*nm)),color=blue, width=3, maxconnect=1, style=solid, "3rd Stage"
  
; ---------------------------------------------------------------------
; ---------------------------------------------------------------------

diagram 18:
;y: 100, 300

if number_of_stages = 2 then
begin
diagram 6:

"Pulse duration vs. compressor length"

x: 0.1, 1.5
"fiber compressor length [m]", @x
y: 10,1500 
"pulse duration [fs]", @y
frame
hx
hy

! begin {only calculated if needed}
  setcolor(blue);
  for x:=CS_x1 to CS_x2 step 0.1
  do begin     
     Find_SteadyState2(x);
     point(x + i* tau_p()/fs, "O");
     end;
  end


! set_device(4);
! set_L(L_2comp)  { restore the original fiber length }
end

(* Gain of ASE
diagram 2:
dBm(P):=10*lg(P/1 mW)
y: dBm(P_l(l1_ASE))-20,dBm(P_l(lMean))+50;{-60,20{0,.0003 1.01*P_l(lMean)/mW} 

plot_gain_ase2(x):= gain_ase2~~[x]
plot_gain_amp2(x):= gain_amplifier2~~[x]

f: plot_gain_ase2(x*nm),{yscale=2,}
  color=olive, width=3, maxconnect=1, style=fdashed, "ASE amplification (2nd Stage)" 

f: plot_gain_amp2(x*nm),{yscale=2,}
  color=olive, width=3, maxconnect=1, style=solid, "gain signal (2nd Stage)" 
*)  

(* Not working. Probably because of the definition of the the fiberst than I need to recover constanly
diagram 17:
"Blue: 1st Stage -- Red: 2nd Stage"

x: -1, 9
"Number of device (fiber)", @x
y: 0, 1800
"Pulse duration [fs]", @y
frame
hx
hy
     
! begin 
  setcolor(red);
  point(4 + i* tau_1stStage/fs, "P");
  point(8 + i* tau_2ndStage/fs, "O");
  set_device(4);
  set_L(.1e-2); 
  Find_SteadyState2();
  point(7 + i* tau_p()/fs, "O");
  set_device(3);
  set_L(.1e-2); 
  Find_SteadyState2();
  point(6 + i* tau_p()/fs, "O");
  set_device(6);
  set_L(.1e-2); 
  Find_SteadyState2();
  point(5 + i* tau_p()/fs, "O");
  end

! set_device(1)
! set_L(L_2pm)

! set_device(6);
! set_L(L_2f)  { restore the original fiber length }

! set_device(3);
! set_L(L_2pm2)  { restore the original fiber length }

! set_device(4);
! set_L(L_2comp)  { restore the original fiber length }

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
*)

