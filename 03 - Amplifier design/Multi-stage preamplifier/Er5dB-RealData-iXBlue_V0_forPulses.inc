include "Er5dB-RealData-iXBlue.inc"   { load the spectroscopic data of the fiber, here A2918M1F1C020R01}

r_core:=0.5*10.1 um  { core radius (taken as MFD) }
;NA_core:=0.22       { numerical aperture }
;r_clad:=125 um/2    { cladding radius }
;w_m(a,V):=a*(0.65+1.619/V^1.5+2.879/V^6)  { Marcuse formula for mode field radius }
;w(l):=w_m(r_core,(2pi/l)*r_core*NA_core)  { wavelength-dependent mode field radius }
;w_s:=w(1550 nm)                          { mode field radius of signal, not used, use measured instead }

w_s:=0.5 * 10.1 um                        { from MFD given by Thorlabs (usual values) }
w_p:=0.5 * 6.5 um                         { mode field radius of pump given by iXblue }
I_s(r):=exp(-2*(r/w_s)^2)                 { intensity distribution of signal field }
I_p(r):=(r <= r/w_p)                      { intensity distribution of pump field }

s_a_Er(l):=6.72e-25*abs_Er(l)/abs_Er(1530 nm)  {state absorption}
s_e_Er(l):=6.55e-25*em_Er(l)/em_Er(1530 nm)    {state emission}

xi_s:=int(r*I_s(r), r:=0 to r_core step 0.1 um)/int(r*I_s(r), r:=0 to 3*r_core step 0.1 um)  { overlap factor for signal wave }
abs_1530:=abs_Er(1530 nm)
N_Er:=(abs_1530/4.34)/(s_a_Er(1530 nm)*xi_s)       { calculate the dopant concentration}
                                
show "Measured N_Er: ", N_Er:d3:"/m^3"    { show the dopant concentration }

N_Er:=(abs_1530/4.34)/(s_a_Er(1530 nm)*xi_s)       { calculate the dopant concentration}
;show "Calculated N_Er: ",N_Er:d3:"/m^3"                                { show the dopant concentration }
;show "overlap factor:  ", xi_s:d3

loss_p:= 0.010  { parasitic losses of pump wave in dB/m }
loss_s:= 0.30  { parasitic losses of signal and ASE waves in dB/m }

;active fiber parameter, 
A_eff_gain:=pi*r_core^2 {The effective surface of the active core}
n2_gain:=2.2e-20        {n2_fused_silica = 2.2e-20 @ 1570 nm}

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