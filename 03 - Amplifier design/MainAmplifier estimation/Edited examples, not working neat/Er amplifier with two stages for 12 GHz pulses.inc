(*
Demo for program "RP Fiber Power":
two-stage ytterbium-doped single-mode fiber amplifier for nanosecond
pulses.
*)
diagram shown: 5


P_p1 := LD_1*mW
P_p2 := LD_2*mW


; 1: "Powers vs. Position in Stage 1"
; 2: "Powers vs. Position in Stage 2"
; 3: "Powers vs. Position in Stage 3"
; 4: "Amplification of a Single Pulse"
; 5: "Output Spectrum"
; 6: "Variation of Pulse Repetition Rate"

; Data for the used fibers:

if GainFiber1 = 5 then include "Er5dB-RealData-iXBlue_V0_forPulses.inc"
if GainFiber1 = 22 then include "Er22dB-RealData-iXBlue_V0.inc"
if GainFiber1 = 80 then include "Er80dB-RealData-iXBlue_V0.inc"

; Fiber 1:
N_Er1:= 1.9e25 { 2.5e24 Er concentration estimation for 5dB}
r_co1:=5 um  { core radius }
w_p1:=4.5 um  { pump mode radius }
w_s1:=5.1 um  { signal mode radius }
I_p1(r):=exp(-2*(r/w_p1)^2)  { pump intensity profile }
I_s1(r):=exp(-2*(r/w_s1)^2)  { signal intensity profile }
NoModes_ASE1:=2

; Fiber 2 (a double-clad fiber):
N_Er2:=1.9e25  { Er concentration }
r_co2:=5 um  { core radius }
w_p2:=4.50 um  { pump mode radius }
w_s2:=5.10 um  { signal mode radius }
I_p2(r):=(r <= w_p2)  { pump intensity profile }
I_s2(r):=exp(-2*(r/w_s2)^2)  { signal intensity profile }
NoModes_ASE2:=2

No_z_steps:=50  { no steps along the fiber }

; Parameters of the channels:
l_p1:=976 nm  { pump wavelength in fiber 1 }
l_p2:=976 nm  { pump wavelength in fiber 2 }
l_s:=1542 nm  { signal wavelength }


; ASE channels:
l1_ASE:=1500 nm
l2_ASE:=1584 nm
dl_ASE:=1 nm
defarray c_ASE_fw1[l1_ASE, l2_ASE, dl_ASE]
defarray c_ASE_fw2[l1_ASE, l2_ASE, dl_ASE]
defarray c_ASE_bw1[l1_ASE, l2_ASE, dl_ASE]
defarray c_ASE_bw2[l1_ASE, l2_ASE, dl_ASE]

P_ASE_fw1(z):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_ASE_fw1[l], z))
P_ASE_fw2(z):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_ASE_fw2[l], z))

P_ASE_bw1(z):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_ASE_bw1[l], z))
P_ASE_bw2(z):=sum(l:=l1_ASE to l2_ASE step dl_ASE, P(c_ASE_bw2[l], z))

; Function for defining the model
; with steady-state power levels:
def_model():=
  begin
  global allow all;
  { Fiber 1: }
  set_device(1);
  set_fiber(L_f1, No_z_steps, 'Er');
  add_ring(r_co1, N_Er1);
  pump1:=addinputchannel(P_p1, l_p1, 'I_p1', 0, dir_p1);
  signal1:=addinputchannel(0, l_s, 'I_s1', 0, forward);
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     c_ASE_fw1[l]:=addASEchannel(l, dl_ASE, NoModes_ASE1, 'I_s1', 0, forward);
     c_ASE_bw1[l]:=addASEchannel(l, dl_ASE, NoModes_ASE1, 'I_s1', 0, backward);
     end;
  finish_fiber();
  { Fiber 2: }
  set_device(2);
  set_fiber(L_f2, No_z_steps, 'Er');
  add_ring(r_co2, N_Er2);
  pump2:=addinputchannel(P_p2, l_p2, 'I_p2', 0, dir_p2);
  signal2:=addinputchannel(0, l_s, 'I_s2', 0, forward);
    { ignored the small background loss }
  for l:=l1_ASE to l2_ASE step dl_ASE
  do begin
     c_ASE_fw2[l]:=addASEchannel(l, dl_ASE, NoModes_ASE2, 'I_s2', 0, forward);
     c_ASE_bw2[l]:=addASEchannel(l, dl_ASE, NoModes_ASE2, 'I_s2', 0, backward);
     end;
  finish_fiber();
  end;

calc def_model()

; Bandpass filter between stages 1 and 2:
T_f(l):=
  { Gaussian transmission function }
  0.9*exp(-(4*ln(2)*(l-l_f)/dl_f)^2)
l_f:=1542 nm  { center wavelength }
dl_f:=100 nm  { bandwidth }

defarray P_ASE_12[l1_ASE, l2_ASE, dl_ASE]
defarray P_ASE_23[l1_ASE, l2_ASE, dl_ASE]

set_signal_in(P_in):=
  { set a certain input power for stage 1
    and transfer the output into stage 2 }
  begin
  global allow all;
  var P12,P23;
  {First Stage}
  set_device(1);
  set_P_in(signal1, P_in);
  P12:=P_out(signal1);
  for l:=l1_ASE to l2_ASE step dl_ASE
  do P_ASE_12[l]:=P_out(c_ASE_fw1[l]);
  {Second Stage}
  set_device(2);
  set_P_in(signal2, P12*T_f(l_s));
  for l:=l1_ASE to l2_ASE step dl_ASE
  do set_P_in(c_ASE_fw2[l], P_ASE_12[l]*T_f(l));
  P23:=P_out(signal2); { get it calculated }
  P_out(signal2); { get it calculated }
  end

calc set_signal_in(0 mW) { get it calculated }


; Parameters of dynamic simulation:
P_s_p:=160 mW  { peak power }
tau_s:=100 fs  { signal pulse duration }

t_p:=0.6*tau_s
{P_s_in(t):=P_s_p*(abs(t-t_p) <= tau_s/2);;P_s_in(t):=P_s_p*exp(-((t-t_p)/(tau_s/2))^30)}
P_s_in(t):= P_s_p*(2/(exp(t/tau_s)+exp(-t/tau_s)))^2 {sech2(t/t_s)}
  { time-dependent signal input power }
T_range:=20*tau_s  {1.2*tau_s time window for simulation of pulse amplification }
dt:=10 fs  { temporal step size }
E_s_in:=int(P_s_in(t), t:=0 to T_range step dt)

f_rep:=12.13 GHz  { repetition rate }
dt_p:=10 fs  { temporal step size for pumping phases }

E_input := P_s_p*tau_s/0.88
Average_Power := E_input*f_rep
;show "Comb Power (Avg):", Average_Power:d3:"mW";
;show E_input

defarray P_s_out1[0, T_range, dt]
store_P_out1():=
  begin
  global allow all;
  set_device(1);
  for t:=0 to T_range step dt
  do P_s_out1[t]:=P_out_dyn(signal1, t);
  end
  
defarray P_s_out2[0, T_range, dt]
store_P_out2():=
  begin
  global allow all;
  set_device(2);
  for t:=0 to T_range step dt
  do P_s_out2[t]:=P_out_dyn(signal2, t);
  end


defarray P_ASE_out[l1_ASE, l2_ASE, dl_ASE]
  { average ASE output powers from the last stage
    during the pumping cycle, see do_pumping() }

do_pumping(T, last_stage):=
  { simulate pumping for some time }
  begin
  global allow all;
  set_device(1);
  set_P_in_dyn(signal1, '');
  calc_dyn(0, T, dt_p);
  if last_stage >= 2
  then begin
       set_device(2);
       set_P_in_dyn(signal2, '');
       calc_dyn(0, T, dt_p);
       end;
  end

pulsegridsize:=2^5;
lgrid:=1542 nm;
tau0:= tau_s;
E0 := E_input;
chirp0 := 0;


do_pulse_amp(last_stage):=
  { simulate the amplification of one signal pulse in all stages }
  begin
  global allow all;
  set_device(1);
  {set_P_in_dyn(signal1, 'P_s_in(t)');}
  set_pulse_grid(T_range,pulsegridsize,lgrid);
  startpulse_s(E0,tau0,chirp0);           {initial pulse as sech2 pulse }
  store_pulse(0);
  calc_dyn(0, T_range, dt);
  store_P_out1();
  if last_stage >= 2
  then begin
       set_device(2);
       set_P_in_dyn(signal2, 'P_s_out1~[t]*T_f(l_s)');
       calc_dyn(0, T_range, dt);
       end; 
	   store_pulse(1){store_P_out2();}
  end


find_ss(last_stage):=
  { simulate multiple amplification and pumping cycles
    (for stages 1..last_stage)
    until the steady state is reached;
    return the number of iterations required }
  begin
  global f_rep, T_range;
  var j, E_out, E_out_l;
  j:=0;
  E_out:=0;
  repeat
    inc(j);
    do_pumping(1/f_rep, last_stage);
    E_out_l:=E_out;
    do_pulse_amp(last_stage);
    E_out:=E_s_out(last_stage);
    {showmessage("j = " + str(j) + ": E_out = " + str(E_out:d3:"J"));}
  until j >= 2 and abs(E_out - E_out_l) < maxr(0.1 uJ, 1e-3*E_out);
  j;
  end
  



E_s_out(j):=
  { signal output energy of stage j }
  begin
  global signal1, signal2, signal3, T_range, dt;
  var ch;
  if j = 1
  then ch:=signal1
  else if j = 2
  then ch:=signal2
  else ch:=signal3;
  set_device(j);
  int(P_out_dyn(ch, t), t:=0 to T_range step dt);
  end

gain(j):= {Used only in the pulse repetition rate variation plots}
  { saturated signal gain of stage j }
  begin
  global E_s_in;
  var r;
  if j = 1
  then r:=E_s_out(1)/E_s_in
  else if j = 2
  then r:=E_s_out(2)/E_s_out(1);
  10*lg(r);
  end

dBm(P):=10*lg(P/(1 mW))


LastStage := 0
if (L_f2>0.1) then LastStage := 2;
if (L_f1>0.1 and L_f2<0.1) then LastStage := 1;

show "Studying Stage ",LastStage



calc find_ss(LastStage)
print_outputs := 12

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

(*
if GainFiber2 > 2 then show "Difference with ASE: ",(dBm(P_l(lMean-12nm))-ASE2_1530):d3:"dBm";
if GainFiber2 <= 2 then show "Difference with ASE: ",(dBm(P_l(lMean))-ASE2_1000):d3:"dBm";
*)
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
y: dBm(P_l(l1_ASE))-10,dBm(P_l(lMean))+20;{-60,20{0,.0003 1.01*P_l(lMean)/mW}
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




