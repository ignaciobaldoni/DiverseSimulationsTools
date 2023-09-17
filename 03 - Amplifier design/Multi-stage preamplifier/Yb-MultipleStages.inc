include "Yb-germanosilicate.inc"   { load the spectroscopic data of the fiber, here A2918M1F1C020R01}
;include "Liekki Yb.inc"

l1_Yb:=850
l2_Yb:=1140
dl_Yb:=1


r_core:=0.5*4.3 um  { core radius (taken as MFD) }
NA_core:=0.22       { numerical aperture }
r_clad:=125 um/2    { cladding radius }

w_m(a,V):=a*(0.65+1.619/V^1.5+2.879/V^6)  { Marcuse formula for mode field radius }
w(l):=w_m(r_core,(2pi/l)*r_core*NA_core)  { wavelength-dependent mode field radius }
;w_s:=w(1050 nm)                          { mode field radius of signal, not used, use measured instead }
w_s:=0.5 * 10.1 um                        { from MFD given by Thorlans (usual values) }
w_p:=0.5 * 6.5 um                         { mode field radius of pump given by iXblue }
I_s(r):=exp(-2*(r/w_s)^2)                 { intensity distribution of signal field }
I_p(r):=(r <= r/w_p)                      { intensity distribution of pump field }


; XXX check the scaling: probably slightly lower values
; in order to obtain a realistic radiative lifetime

n:=1.45 {Dispersion will be calculated elsewhere, this is for radiative lifetime estimate}
;tau_rad:=1/(8*pi*n^2*c*int(s_e_Yb(l)/l^4, l:=l1_Yb*nm to l2_Yb*nm step 1e-9)) {Appears to be the radiative lifetime}
;show "tau_rad: ", tau_rad:d3:"s"

xi_s:=int(r*I_s(r), r:=0 to r_core step 0.1 um)/int(r*I_s(r), r:=0 to 3*r_core step 0.1 um)  { overlap factor for signal wave }                                

show "Estimated N_Yb:      ", N_Yb:d3:"/m^3"    { show the dopant concentration }

;show "overlap factor:      ", xi_s:d3

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

; Calculate dispersion in units of fs/nm/m from D_Yb_f[l] provided in fiber spectroscopic data:
GVD1_f(l) :=
  if ld1 <= l / nm <= ld2 then
    -l^2 / (2pi * c)* (D_Yb_f~[l / nm] * fs / nm)
  { use this function in main file to define the GVD of a fiber: set_GVD_l("GVD_f(l)") }

; For plotting dispersion parameter in units of fs/nm/m from D_Yb_f[l] provided in fiber spectroscopic data:
DD1_f(l) :=
  if ld1 <= l / nm <= ld2 then
    (D_Yb_f~[l / nm])

; Calculate dispersion in units of fs/nm/m from D_Yb_f[l] provided in fiber spectroscopic data:
GVD_f(l) :=
  if ld1 <= l / nm <= ld2 then
    -l^2 / (2pi * c) * .89* (D_Yb_f~[l / nm] * fs / nm)
	

  { use this function in main file to define the GVD of a fiber: set_GVD_l("GVD_f(l)") }

{1.92 in order to get 17.7 fs/nm m in our fiber}

; For plotting dispersion parameter in units of fs/nm/m from D_Yb_f[l] provided in fiber spectroscopic data:
DD_Yb_f(l) :=
  if ld1 <= l / nm <= ld2 then
    .89* (D_Yb_f~[l / nm])

{.89 in order to get -17.7 fs/nm m in our fiber}


diagram 17:

"Yb absortion and emission at 980 nm"
x: l1_Yb , l2_Yb 
"wavelength [nm]", @x
y: 0.01e-25,40e-25
" Cross section[dBm/m]", @y
frame
hx
hy

f: s_a_Yb(x*nm),
  color=olive, width=3, maxconnect=1, style=solid, "Yb absortion"
f: s_e_Yb(x*nm),
  color=olive, width=3, maxconnect=1, style=fdashed, "Yb emission"
  
diagram 11: 

"GVD vs. Frequency"
"22 dB/m"

x: ld1 , ld2 
"wavelength [nm]", @x
y: -50,0
" [dBm/m]", @y
frame
hx
hy

f: DD_Yb_f(x*nm),
  color=olive, width=3, maxconnect=1, style=solid, "DD_Yb_f"
f: DD1_f(x*nm),
  color=teal, width=3, maxconnect=1, style=solid, "DD1_f"
