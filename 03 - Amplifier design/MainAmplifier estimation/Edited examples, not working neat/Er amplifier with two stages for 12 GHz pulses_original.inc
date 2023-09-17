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
T_range:=1.2*tau_s  { time window for simulation of pulse amplification }
dt:=0.5 fs  { temporal step size }
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

do_pulse_amp(last_stage):=
  { simulate the amplification of one signal pulse in all stages }
  begin
  global allow all;
  set_device(1);
  set_P_in_dyn(signal1, 'P_s_in(t)');
  calc_dyn(0, T_range, dt);
  store_P_out1();
  if last_stage >= 2
  then begin
       set_device(2);
       set_P_in_dyn(signal2, 'P_s_out1~[t]*T_f(l_s)');
       calc_dyn(0, T_range, dt);
       end;
	   store_P_out2();
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
  
;
;save_pulse_f('Pulse_F.txt');

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


; -------------
diagram 1:

"Powers vs. Position in Stage 1"

x: 0, L_f1
"position in fiber [m]", @x, size=%6
y: -80,25 {0, maxr(P_p1, 0.1)}
y2: 0, 100
frame
legpos 100, 500
hx
hy

! set_device(1)

f: dBm(P(pump1, x)),
   color=red, width=3, "pump"
;f: P(signal1, x),
   color=blue, width=3, "signal"
f: dBm(P_ASE_fw1(x)),
   color=teal, width=3, "ASE"
f: dBm(P_ASE_bw1(x)), "Backward ASE",
   color=teal, width=3, style=dashed

f: 100*n(x, 2), yscale=2,
   color=black, width=3, style=fdashed, "Er excitation (%, right scale)"

["signal gain at ", l_s:d3:"(n)m", ": ", (set_device(1); sp_gain(signal1)):f1:"dB"]


; -------------
diagram 2:

"Powers vs. Position in Stage 2"

x: 0, L_f2
"position in fiber [m]", @x, size=%6
y: -80,25{0, maxr(P_p2, 0.1)}
y2: 0, 100
frame
legpos 100, 500 {Position of legend}
hx
hy

! set_device(2)

f: dBm(P(pump2, x)),
   color=red, width=3, "pump"
;f: P(signal2, x),
   color=blue, width=3, "signal"
f: dBm(P_ASE_fw2(x)),
   color=teal, width=3, "ASE"
f: dBm(P_ASE_bw2(x)), "Backward ASE",
   color=teal, width=3, style=dashed


f: 100*n(x, 2), yscale=2,
   color=black, width=3, style=fdashed, "Er excitation (%, right scale)"

["signal gain at ", l_s:d3:"(n)m", ": ", (set_device(2); sp_gain(signal2)):f1:"dB"]
["Power: ", (set_device(2); P(signal2,L_f2)):f1:"mW"]
; -------------

; -------------
diagram 4:

"Amplification of a Single Pulse"

x: 0, T_range/fs
"time [fs]", @x
y: 10, 100
;"output power [W]", @y
y2: 20,100
frame
legpos 500, 150
hx
hy

;f: 0.9*CS_y2*P_s_in(x*fs)/P_s_p,
  width=3, style=dotted, "seed power (a. u.)"

; Use one of the two following commmands:
;! do_pulse_amp(LastStage)  { start after pumping for a long time }
! find_ss(LastStage)  { start with steady state for given repetition rate }

f: dBm(P_out_dyn(signal1, x*fs)),
  init set_device(1),
  color=red, width=3, maxconnect=100, "power from stage 1 [dBm]"

["E input = ", E_input:d3:"J"], 150, 150
! E_s1:=E_s_out(1)
["E_s1 = ", E_s1:d3:"J"], 150, 200

f: dBm(P_out_dyn(signal2, x*fs)),{/kW}
  init set_device(2),
  color=green, width=3, maxconnect=100, "power from stage 2 [dBm]"

! E_s2:=E_s_out(2)
["E_s2 = ", E_s2:d3:"J"], 150, 250


f: 100*n_av_dyn(2,x*fs),
  init set_device(1),
  color=red, width=3, yscale = 2,style=dashed{, "Er excitation in stage 1 (%)"}

f: 100*n_av_dyn(2,x*fs),
  init set_device(2),
  color=green, width=3, yscale = 2, style=dashed{, "Er excitation in stage 2 (%)"}


; -------------
diagram 5:

"Output Spectrum"

x: l1_ASE/nm, l2_ASE/nm
"wavelength (nm]", @x
y: -100,100
"average power [dBm in 1 nm]", @y
frame
hx
hy

! find_ss(LastStage)

P_ASE_out_av(l):=
  begin
  global f_rep, T_range, dt_p, LastStage, c_ASE_fw2[];{ c_ASE_fw2[];}
  set_device(LastStage);
  int(P_out_dyn(c_ASE_fw2[l], t), t:=0 to 1/f_rep step dt_p)/(1/f_rep);
  end  
(*
Pulse_output(t):=
  begin
  set_device(LastStage);
  int(P_out_dyn(signal2, t*fs), t:=0 to T_range/fs step 10 fs);
  end  
*)
  
! P_s_av:=E_s_out(LastStage)*f_rep
! do_pumping(1/f_rep, LastStage)

! (setcolor(green); point(l_s/nm + i*dBm(P_s_av), "R"));

f: dBm(P_ASE_out_av(x*nm)),
  color=olive, width=3, maxconnect=1e3


! P_ASE_tot:=sum(l:=l1_ASE to l2_ASE step dl_ASE, P_ASE_out_av(l))
! P_tot:=P_s_av+P_ASE_tot
["Signal output: ", P_s_av:d3:"W",
 "; ASE: ", P_ASE_tot:d3:"W",
 " (", 100*P_ASE_tot/P_tot:f1:"%", " of total power)"]


; -------------
diagram 6:

"Output Spectrum for the three stages"

x: l1_ASE/nm, l2_ASE/nm
"wavelength (nm]", @x
y: -50, +30
"average power [dBm in 1 nm]", @y
frame
hx
hy
legpos 500,150
! find_ss(1){find_ss(2)}

P_ASE_out_av1(l):=
  begin
  global f_rep, T_range, dt_p, c_ASE_fw1[];{ c_ASE_fw2[];}
  set_device(1);
  int(P_out_dyn(c_ASE_fw1[l], t), t:=0 to 1/f_rep step dt_p)/(1/f_rep);
  end

! P_s_av1:=E_s_out(1)*f_rep
! do_pumping(1/f_rep, 1)

! (setcolor(teal); point(l_s/nm + i*dBm(P_s_av1), "R"));



! P_ASE_tot1:=sum(l:=l1_ASE to l2_ASE step dl_ASE, P_ASE_out_av1(l))
! P_tot1:=P_s_av1+P_ASE_tot1

f: dBm(P_ASE_out_av1(x*nm)), 
  color=teal, width=3, maxconnect=1e3
  
["Signal output: ", P_s_av1:d3:"W",
 "; ASE: ", P_ASE_tot1:d3:"W",
 " (", 100*P_ASE_tot1/P_tot1:f1:"%", " of total power)"],150,150

! find_ss(2){find_ss(2)}

P_ASE_out_av2(l):=
  begin
  global f_rep, T_range, dt_p, c_ASE_fw2[];{ c_ASE_fw2[];}
  set_device(2);
  int(P_out_dyn(c_ASE_fw2[l], t), t:=0 to 1/f_rep step dt_p)/(1/f_rep);
  end

! P_s_av2:=E_s_out(2)*f_rep
! do_pumping(1/f_rep, 2)

! (setcolor(olive); point(l_s/nm + i*dBm(P_s_av2), "R"));

f: dBm(P_ASE_out_av2(x*nm)),
  color=olive, width=3, maxconnect=1e3

! P_ASE_tot2:=sum(l:=l1_ASE to l2_ASE step dl_ASE, P_ASE_out_av2(l))
! P_tot2:=P_s_av2+P_ASE_tot2
["Signal output: ", P_s_av2:d3:"W",
 "; ASE: ", P_ASE_tot2:d3:"W",
 " (", 100*P_ASE_tot2/P_tot2:f1:"%", " of total power)"],150,200
 
! find_ss(2){find_ss(2)}


; -------------
(*
diagram 7:

"Variation of Pulse Repetition Rate"

x: 5, 20, 2.5
"pulse repetition rate (kHz)", @x
y: 0, 35
;"average output power [W]", @y
frame
legpos 500, 150
hx
hy

! f_rep_old:=f_rep

! begin
  setcolor(red);
  legend("R", "", "signal average power [W]");
  setcolor(teal);
  legend("N", "", "ASE average power [W]");
  setcolor(orange);
  legend("h", "", "gain of stage 1 [dB]");
  setcolor(green);
  legend("h", "", "gain of stage 2 [dB]");
  for x:=CS_x1 to CS_x2 step 2.5
  do begin
     var P_ASE_tot;
     f_rep:=x*kHz;
     find_ss(2);
     setcolor(red);
     point(x + i*E_s_out(2)*f_rep, "R");
     do_pumping(1/f_rep, 2);
     P_ASE_tot:=sum(l:=l1_ASE to l2_ASE step dl_ASE, P_ASE_out_av(l));
     do_pulse_amp(2);
     setcolor(teal);
     point(x + i*P_ASE_tot, "N");
     setcolor(orange);
     point(x + i*gain(1), "h");
     setcolor(green);
     point(x + i*gain(2), "h");
     end;
  end

! f_rep:=f_rep_old
*)