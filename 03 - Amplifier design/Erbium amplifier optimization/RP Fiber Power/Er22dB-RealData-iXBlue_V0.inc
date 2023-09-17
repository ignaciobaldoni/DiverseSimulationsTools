; Radhard fiber by iX:  A2918M1F1C020R01
; Edited by M Lezius 10-18-2021 for speciality gain fiber provided by iXBlue

include "Er22dB-RealData-iXBlue.inc"   { load the spectroscopic data of the fiber, here A2918M1F1C020R01}

r_core:=0.5*5.0 um  { core radius (taken as MFD) }
NA_core:=0.22       { numerical aperture }
r_clad:=125 um/2    { cladding radius }

w_m(a,V):=a*(0.65+1.619/V^1.5+2.879/V^6)  { Marcuse formula for mode field radius }
w(l):=w_m(r_core,(2pi/l)*r_core*NA_core)  { wavelength-dependent mode field radius }
;w_s:=w(1550 nm)                          { mode field radius of signal, not used, use measured instead }
w_s:=0.5 * 10.1 um                        { from MFD given by IXFiber }
w_p:=0.5 * 5.0 um                         { mode field radius of pump }
I_s(r):=exp(-2*(r/w_s)^2)                 { intensity distribution of signal field }
I_p(r):=(r <= r/w_p)                      { intensity distribution of pump field }

s_a_Er(l):=6.72e-25*abs_Er(l)/abs_Er(1530 nm)  {state absorption}
s_e_Er(l):=6.55e-25*em_Er(l)/em_Er(1530 nm)    {state emission}


; XXX check the scaling: probably slightly lower values
; in order to obtain a realistic radiative lifetime

n:=1.45 {Dispersion will be calculated elsewhere, this is for radiative lifetime estimate}
tau_rad:=1/(8*pi*n^2*c*int(s_e_Er(l)/l^4, l:=l1_1500_Er*nm to l2_1500_Er*nm step 1e-9)) {Appears to be the radiative lifetime}
;show "tau_rad: ", tau_rad:d3:"s"
xi_p:=(r_core/r_clad)^2

xi_s:=int(r*I_s(r), r:=0 to r_core step 0.1 um)/int(r*I_s(r), r:=0 to 3*r_core step 0.1 um)  { overlap factor for signal wave }
abs_1530:=abs_Er(1530 nm)
show "overlap factor:", xi_s:d3;


;N_Er2:=abs_1530/(4.343*s_a_Er(1530 nm)*xi_p)
;show "Calculated N_Er2: ",N_Er1:d3:"/m^3"   { show the dopant concentration }
N_Er:=abs_1530/(5.77*s_a_Er(1530 nm)*xi_s)       { calculate the dopant concentration}
show "dop. actual gain:    ", N_Er:d3:"/m^3"   

N_Er1:= 1.94e25; {1.94e25 Measured value by iXblue}
show "Meas. N_Er @ 30dB:   ", N_Er1:d3:"/m^3"    
;show "overlap factor:     ", xi_s:d3
(*
loss_p:= 0.010  { parasitic losses of pump wave in dB/m }
loss_s:= 0.30  { parasitic losses of signal and ASE waves in dB/m }
*)

;active fiber parameter, 
A_eff_gain:=pi*r_core^2 {The effective surface of the active core}
n2_gain:=2.2e-20        {n2_fused_silica = 2.2e-20 @ 1570 nm}
;GVD_gain:=-(lMean^2 / (2pi*c)) * (-9.5) ps/(nm km)    {OLD: Calculates the GVD of this fiber directly. Assumption: D=-9.5 ps/(nm km)@1570nm}
;show "GVD_gain:  ", GVD_gain/fs^2:f0:"fs^2/m"         {Print the GVD of the gain fiber}
;gamma_gain:=(2*pi/lMean) * n2_gain / A_eff_gain
;show "gamma_gain: ", gamma_gain:d3:"rad/(W m)"

; Calculate dispersion in units of fs/nm/m from D_f[l] provided in fiber spectroscopic data:
GVD1_f(l) :=
  if ld1 <= l / nm <= ld2 then
    -l^2 / (2pi * c)* (D_f~[l / nm] * fs / nm)
  { use this function in main file to define the GVD of a fiber: set_GVD_l("GVD_f(l)") }

; For plotting dispersion parameter in units of fs/nm/m from D_f[l] provided in fiber spectroscopic data:
DD1_f(l) :=
  if ld1 <= l / nm <= ld2 then
    (D_f~[l / nm])

; Calculate dispersion in units of fs/nm/m from D_f[l] provided in fiber spectroscopic data:
GVD_f(l) :=
  if ld1 <= l / nm <= ld2 then
    -l^2 / (2pi * c) * 1.92* (D_f~[l / nm] * fs / nm)
  { use this function in main file to define the GVD of a fiber: set_GVD_l("GVD_f(l)") }

{1.92 in order to get 17.7 fs/nm m in our fiber}

; For plotting dispersion parameter in units of fs/nm/m from D_f[l] provided in fiber spectroscopic data:
DD_f(l) :=
  if ld1 <= l / nm <= ld2 then
    1.92* (D_f~[l / nm])

{1.92 in order to get -17.7 fs/nm m in our fiber}



; Automatically set values in input fields:
;IF_NoRings:=1
;IF_r1:=r_core
;IF_r2:=
;IF_r3:=
;IF_conc11:=N_Er
;IF_conc12:=
;IF_conc21:=
;IF_conc22:=
;IF_conc31:=
;IF_conc32:=
;IF_shape_p1:=
;IF_shape_p2:=
;IF_w_p1:=w_p
;IF_w_p2:=w_p
;IF_shape_s1:=Gaussian
;IF_shape_s2:=Gaussian
;IF_w_s1:=w_s
;IF_w_s2:=w_s
;IF_shape_ASE:=Gaussian
;IF_w_ASE:=w_s
;IF_loss_p1:=loss_p
;IF_loss_p2:=loss_p
;IF_loss_s1:=loss_s
;IF_loss_s2:=loss_s
;IF_loss_ASE:=loss_s


diagram 12: 

"GVD vs. Frequency"
"22 dB/m"

x: 1450 , 1599 
"wavelength [nm]", @x
y: -50,0
" [dBm/m]", @y
frame
hx
hy

f: DD_f(x*nm),
  color=olive, width=3, maxconnect=1, style=solid, "DD_f"
f: DD1_f(x*nm),
  color=teal, width=3, maxconnect=1, style=solid, "DD1_f"


diagram 19:

"Er absortion and emission at 1500 nm"
"22 dB/m"
x: 1450 , 1600 
"wavelength [nm]", @x
y: 0,50
" [dBm/m]", @y
frame
hx
hy

f: abs_Er(x*nm),
  color=olive, width=3, maxconnect=1, style=solid, "Er absortion"
f: em_Er(x*nm),
  color=olive, width=3, maxconnect=1, style=fdashed, "Er emission"


diagram 20:

"Er absortion and emission at 980 nm"
"22 dB/m"
x: 913 , 1100 
"wavelength [nm]", @x
y: 0,50
" [dBm/m]", @y
frame
hx
hy

f: abs_Er(x*nm),
  color=olive, width=3, maxconnect=1, style=solid, "Er absortion"
f: em_Er(x*nm),
  color=olive, width=3, maxconnect=1, style=fdashed, "Er emission"
