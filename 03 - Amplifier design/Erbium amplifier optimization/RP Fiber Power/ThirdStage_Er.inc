include "SecondStage_Er.inc"

show "------- Third stage: -------"

tau_2ndStage:= tau_p()
PeakPower_2ndStage := P_p()
Avg_power_2nd := E_p()*f_rep;


P_signal_fw_in := P_signal_fw_out;{Avg_power_2nd;}


real_fiber22dB := gain_fiber_3

if GainFiber3 = 5 then include "Er5dB-RealData-iXBlue_V0.inc"
if GainFiber3 = 22 then include "Er22dB-RealData-iXBlue_V0.inc"
if GainFiber3 = 80 then include "Er80dB-RealData-iXBlue_V0.inc"



; Pump laser intensities
P_pump5_bw_in:= Pump_power5 * mW       {Define here the pump power and direction of Pump 1}


loss_s := Loss_S;
loss_p := Loss_P;

; Parameters of the channels:
;l_p5:=976 nm                 {Wavelength pump 1}
w_p5:=w_p                    {Mode Radius pump 1}
I_p5(r):=exp(-2*(r/w_p5)^2)  {Intensity profile pump 1}

l_s1:=lMean                  {Signal wavelength}
w_s1:= w_s               	 {3.25 um Signal mode radius}

I_s1(r):=exp(-2*(r/w_s1)^2)  {Signal mode profile}
w_ASE:=w_s1                   {ASE mode radius}
I_ASE(r):=exp(-2*(r/w_ASE)^2) {ASE mode profile}





; Array of ASE channels:
l1_ASE:= lMean - 42nm               {ASE lower wavelength}
l2_ASE:= lMean + 58 nm               {ASE upper wavelength}
dl_ASE:=1 nm                  {ASE wavelength steps}
NoModes_ASE:=2                {Two ASE modes: forward and backward}
defarray c_3ASE_fw[l1_ASE, l2_ASE, dl_ASE]   {ASE forward array definition}
defarray c_3ASE_bw[l1_ASE, l2_ASE, dl_ASE]   {ASE backward array definition}
P3_ASE_fw(x):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_3ASE_fw[l],x))  {ASE forward power at pos x integrated over wavelength}
P3_ASE_bw(x):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_3ASE_bw[l],x))  {ASE backward power at pos x integrated over wavelength}


; ---------
{Script starts}


; Define the fiber models:
def_model():=
  begin
  global allow all;
  set_device(7);
  if GainFiber3 > 2 then set_fiber(L_3f, No_z_steps, 'Er'); {Length of the gain fiber segment}
  add_ring(r_core, N_Er);           {Fill core with Erbium from fiber include file}
  max_N:=N_Er;
  max_r:=r_core;
  P_in_max:=0;
  P_in_max:=maxr(P_in_max, P_pump5_bw_in);
  pump5_bw:=addinputchannel(P_pump5_bw_in, l_p5, 'I_p5', loss_p, dir_p5);
  P_in_max:=maxr(P_in_max, P_signal_fw_in);
  signal3_fw:=addinputchannel(P_signal_fw_in, l_s1, 'I_s1', loss_s, forward);
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     c_3ASE_fw[l]:=addASEchannel(l, dl_ASE, NoModes_ASE, 'I_ASE', loss_s, forward);
     c_3ASE_bw[l]:=addASEchannel(l, dl_ASE, NoModes_ASE, 'I_ASE', loss_s, backward);
     end;
  finish_fiber();
  end;
calc def_model()


; These two functions are defined for the maximum ASE channel calculation
P3_ASE_fw_out(l):=
  if l1_ASE <= l <= l2_ASE
  then P_out(c_3ASE_fw[l])+P2_ASE_fw_l(l)
P3_ASE_bw_out(l):=
  if l1_ASE <= l <= l2_ASE
  then P_out(c_3ASE_bw[l])+P2_ASE_bw_l(l)


;Outputs from cw regime initial simulation

(*
show "----------------------------------------"
show "Initial powers:"
show "Average Power:      ", P_signal_fw_in:d3:"W"
show "Pump 1 (backward):  ", P_pump5_bw_in:d3:"W"
show "----------------------------------------"
show "Output powers:"
*)
calc set_device(7)

;--------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------
diagram 5:
"Powers vs. Position"
"Third Stage"

x: 0, L_3f
"position in fiber [m]", @x
y: 0, P_out(signal3_fw)*1.3/mW
"Signal Power [mW]",@y
y2: 0, P_pump5_bw_in*(1.3)/mW
"Pump Power [mW]",@y2
frame
legpos 200, 250
hx
hy

defarray P3_fiber[0,L_3f,0.01]
defarray Pump3_fiber[0,L_3f,0.01]
P3_signal_fiber():= 
  begin
  global allow all;
  for z:=0 to z:= L_3f step 0.01
  do begin
     P3_fiber[z]:=P(signal3_fw,z);
	 Pump3_fiber[z]:=P(pump5_bw,z);
     end;
  end
calc P3_signal_fiber();
P3_fiber_z(x):=P3_fiber~~[x]
Pump3_fiber_z(x):=Pump3_fiber~~[x]
f: P3_fiber_z(x)/mW,{yscale=2,}
  color=blue, width=3, maxconnect=1, style=solid, "Signal amplification (3rd Stage)"  
f: Pump3_fiber_z(x)/mW,yscale=2,
  color=cyan, width=3, maxconnect=1, style=fdashed, "Pump Power (3rd Stage)"  

defarray gain_ase3[l1_ASE,l2_ASE,dl_ASE]
defarray gain_amplifier3[l1_ASE,l2_ASE,dl_ASE]
GainASE3():=
	begin
	global allow all;
	for l:= l1_ASE to l2_ASE step dl_ASE
	do begin
		gain_ase3[l]:=sp_gain(c_3ASE_fw[l]);
		gain_amplifier3[l]:=sp_gain(signal3_fw)
		end;
	end
calc GainASE3();

 
;--------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------

P_pump5_bw_out:=P_out(pump5_bw)
;show "pump1 (backward):   ", P_pump5_bw_out:d3:"W"
P_signal_fw_out:=P_out(signal3_fw)
;show "signal (forward):   ", P_signal_fw_out:d3:"W"
G_signal3:=sp_gain(signal3_fw)
;show "G_signal3:          ", G_signal3:d3:np:"dB"
P3_ASE_fw:=sum(l:=l1_ASE to l2_ASE step dl_ASE, ASE_factor*P_out(c_3ASE_fw[l])+P2_ASE_fw_l(l)*10^(gain_ase3[l]/10))
P3_ASE_bw:=sum(l:=l1_ASE to l2_ASE step dl_ASE, ASE_factor*P_out(c_3ASE_bw[l])+P2_ASE_bw_l(l)*10^(gain_ase3[l]/10))
calc (G_ASE:=0; for l:=l1_ASE to l2_ASE step dl_ASE do G_ASE:=maxr(G_ASE, sp_gain(c_3ASE_fw[l])))
;show "P3_ASE_fw:           ", P3_ASE_fw:d3:"W"
;show "P3_ASE_bw:           ", P3_ASE_bw:d3:"W"




if print_outputs = 1 or print_outputs = 12 then 
begin 
show "NF (3rd):            ", 10*lg(NF(signal3_fw)):d3:"dB"
SNR_out3 := 10*lg(P_signal_fw_out/P3_ASE_fw)
show "SNR_out(3rd):        ", SNR_out3:d3:"dB"
show "ASE3:                ", 100*P3_ASE_fw/(P3_ASE_fw+P_signal_fw_out):d3:"%"
end
;--------------------------------------------------------------------------------------
diagram 2:

dBm(P):=10*lg(P/1e-3)

; First store the wavelength-dependent ASE powers in arrays,
; so that we can later get interpolated values:
defarray P3_ASE_fw[l1_ASE,l2_ASE,dl_ASE]
defarray P3_ASE_bw[l1_ASE,l2_ASE,dl_ASE]
StoreASE3():=
  begin
  global allow all;
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     P3_ASE_fw[l]:=ASE_factor*P_out(c_3ASE_fw[l])+P2_ASE_fw_l(l)*10^(gain_ase3[l]/10);{(G_signal3/10)}
     P3_ASE_bw[l]:=ASE_factor*P_out(c_3ASE_bw[l])+P2_ASE_bw_l(l)*10^(gain_ase3[l]/10);
     end;
  end
calc StoreASE3()
P3_ASE_fw_l(l):=P3_ASE_fw~~[l] 
;show "Forward ASE @ CenterWavelength: ",dBm(P3_ASE_fw_l(lMean)*(nm/dl_ASE)):d3:"dB"
P3_ASE_bw_l(l):=P3_ASE_bw~~[l]
(*
plot_gain_ase3(x):= gain_ase3~~[x]
plot_gain3(x):=gain_amplifier3~~[x]
f: plot_gain_ase3(x*nm),{yscale=2,}
  color=cyan, width=3, maxconnect=1, style=fdashed, "ASE amplification (3rd Stage)" 
f: plot_gain3(x*nm),{yscale=2,}
  color=cyan, width=3, maxconnect=1, style=solid, "ASE amplification (3rd Stage)"
*)  
f: dBm(P3_ASE_fw_l(x*nm)), {violet}
  color=blue, width=3, maxconnect=1, style=fdashed, "forward ASE (3rd Stage)"
;f: dBm(P3_ASE_bw_l(x*nm)),  color=violet, width=3, maxconnect=1, style=fdashed, "backward ASE (3rd Stage)"

ASE3_1530 :=  dBm(P3_ASE_fw_l(lMean-12nm));
;show "ASE @ ",lMean-12nm, ": ", ASE3_1530:d3:"dB"
ASE3_1000 :=  dBm(P3_ASE_fw_l(lMean));


; ---------------------------------------------------------
; ------ Dynamic simulations for pulse amplification ------
; ---------------------------------------------------------

;show "---------------------------------"
;show "-- 3rd Stage PULSE SIMULATIONS --"
;show "---------------------------------"
;show "Gain Fiber (3rd):     ", GainFiber3:d3:"dB/m"

;gridsize determination for the laser pulse
lgrid:=lMean {covers XXX nm to XXX nm if T_range and pulsegridsize are unchanged}
T_range:=5 ps  {Enough for pulse duration description}
pulsegridsize:=2^10  {NOT clear what this is, probably clarified in RP Fiber Power Manual}

; Ultrashort pulses:
f_rep:= Repetition_rate * GHz       { <- set here Signal in repetition rate}
P_av:= Avg_power_2nd    {P_signal_fw_in <- set here Signal in power}
tau2nd:= tau_2ndStage;           {<- set here Signal pulse duration at injection}
chirp0:= 0 GHz/ps;              {Unnecessary given it is taken according to last result}
E0:=P_av/f_rep           {Pulse energy varies determination}
l0_p:=lMean              {Center wavelength of the pulse is set to mean wavelenth of it}


DoPulsePropagation3_PM_fiber() :=
  begin  {Pulse propagation in PM fiber}
  global allow all;
  { passive PM fiber: }
  load_pulse("Pulse_from_second_stage.txt");
  {recall_pulse(9);}
  set_device(1);
  set_L(L_3pm);
  set_GVD(GVD_pm);
  {if GVD0 = 0 then set_GVD(0);}
  set_n2(n2_pm);
  pp_fiber(1, signal_pm);{Do the pulse propagation}
  store_pulse(10);{save the pulse}
  end

calc DoPulsePropagation3_PM_fiber()

if cw_laser = 0 then calc (set_device(7); set_P_in(signal3_fw, 0)) 
;show "cw laser OFF"

DoPulsePropagation3_amplifier():=
  begin  {Pulse propagation in amplifier fiber: }
  global allow all;
  set_device(7);{Device 2 is the gain fiber!}
  set_ss(1);{Self-steppening is considered!}
  set_dnr("h_R(t)",f_R); {Raman is considered!}
  set_GVD_l("GVD_f(l)");
  {if GVD0 = 0 then set_GVD(0);}
  set_n2(n2_gain {XXX});  {set n2 of the Gain fiber as calculated above}
  pp_fiber(7, signal3_fw); {Propagates the pulse through the fiber 2, using the optical channel "signal3_fw"}
  store_pulse(11);         {save the pulse}
  if maxr(P_t(-0.95*T_range/2), P_t(+0.95*T_range/2)) > 0.01 * P_p()
  then error("Time range too small!?!");
  end

;calc DoPulsePropagation3_amplifier(){no steady state!}
     { switch off the cw input - we now look at pulses }


DoPulsePropagation3_postamp():=
  begin  {Pulse propagation in amplifier fiber: }
  global allow all;
  set_device(3);
  set_L(L_3pm2);
  set_GVD(GVD_pm);
  {if GVD0 = 0 then set_GVD(0);}
  set_n2(n2_pm);
  startpulse_recall(11);
  pp_fiber(3, signal_pm2);{Do the pulse propagation}
  store_pulse(12);{save the pulse}
  end

if cw_laser = 0 then calc (set_device(3);set_P_in(signal_pm2,0))

{GVD_compressor := -17}

DoCompressor3_postamp(x):=
  begin  {Pulse propagation in amplifier fiber: }
  global allow all;
  set_device(4);
  set_L(x);
  set_GVD(GVD_compressor);
  {if GVD0 = 0 then set_GVD(0);}
  set_n2(n2_pm);
  startpulse_recall(12);
  pp_fiber(4, signal_comp);{Do the pulse propagation}
  store_pulse(13);{save the pulse}
  end

if cw_laser = 0 then calc (set_device(4);set_P_in(signal_comp,0))


DoPumping3(T):= {simulates the influence of the pump wave only for a given time}
  begin
  global T;
  calc_dyn(0, T, minr(10 us, T));
  end

f_rep_ss:= 100 kHz;   { Then, influence of the pump is only for 10 µs. Don't know where does this value come from }

Find_SteadyState3(x):=
  { do multiple amplification and pumping cycles until the steady state is reached;
    return the number of iteration steps }
  begin
  global f_rep, E0, tau2nd, chirp0, f_rep_ss;
  var E_out, E_out_l, N;
  
  set_gain_sat(f_rep/f_rep_ss); { effectively increased gain saturation per pulse }
  E_out:=0;
  N:=0;
  repeat
    recall_pulse(10);
    inc(N); {N = N+1}
    E_out_l:=E_out;
    DoPulsePropagation3_amplifier();
    { store_pulse(10+N); {}
    DoPumping3(1/f_rep_ss);
    E_out:=E_p();
  until abs(E_out-E_out_l) < 1e-3 * E_out;
  N; 
  DoPulsePropagation3_postamp();
  DoCompressor3_postamp(x);
  save_pulse_t("Pulse_from_third_stage.txt")
  end
(*
show "Initial pulse energy: ", E0:d3:"J"
show "Initial T pulse :     ", tau2nd:d3:"s"
show "Initial Peak Power:   ", PeakPower_2ndStage:d3:"W"
show "Initial bandwidth:    ", 0.314/tau2nd:d3:"Hz"
;show "Initial chirp:        ", (chirp0*1e-12):d3:"Hz/ps"
*)

calc Find_SteadyState3(L_3comp)

; Display some parameters of the obtained pulse:
;show "Initial Avg Power: " P_av:d3:"W"
;show "Bandwidth:           ", dl_p():d3:"(n)m"
;show "TBP:                 ", TBP():d3

if print_outputs = 2 or print_outputs = 12 then 
begin
show "Energy per pulse:    ", E_p():d3:"J"
show "Peak power:          ", P_p():d3:"W" {Peak power}
show "T pulse:             ", tau_p():d3:"s"
tau_3rdStage := tau_p()
show "Geff:                ", 10*lg(E_p()/E0):d3:"dB"
Geff3:= 10*lg(E_p()/E0);
show "Output avg power :   ", E_p()*f_rep:d3:"W"
TotalGain:= Geff3+Geff2+Geff1
show "Overall Geff:        ", TotalGain:d4:"dB"
end

if GainFiber3 > 2 then show "Difference with ASE: ",(dBm(P_l(lMean-12nm))-ASE3_1530):d3:"dBm";
if GainFiber3 <= 2 then show "Difference with ASE: ",(dBm(P_l(lMean))-ASE3_1000):d3:"dBm";


;-------------------------------------------------------------------------------
;---------- From here on, the description of the different diagrams of interest
;-------------------------------------------------------------------------------

diagram 1:
if Pulse_dBm = 0 then
begin
y: 0,P_t(0)*1.1;
end

diagram 2:
dBm(P):=10*lg(P/1 mW)
y: dBm(P_l(l1_ASE))-20,dBm(P_l(lMean))+50;{-60,20{0,.0003 1.01*P_l(lMean)/mW} 

diagram 18:
;y: 100, 300

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------

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
     Find_SteadyState3(x);
     point(x + i* tau_p()/fs, "O");
     end;
  end

! set_device(4);
! set_L(L_3comp)  { restore the original fiber length }

