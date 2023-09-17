;PM DCF from Thorlabs

ld1 := 1000
ld2 := 1650
dld := 50
defarray PM_f[ld1,ld2,dld]

readlist l, PM_f[l]:
1000.0, -44
1050.0, -34
1100.0,	-27
1150.0,	-20
1200.0,	-14
1250.0, -7.
1300.0, -1.
1350.0,	4.
1400.0,	8.
1450.0,	12.
1500.0,	15.
1550.0,	17.5
1600.0, 19.5
1650.0, 30.

; Calculate dispersion in units of fs/nm/m from D_f[l] provided in fiber spectroscopic data:
GVD_pm_f(l) :=
  if ld1 <= l / nm <= ld2 then
    -l^2 / (2pi * c)* (PM_f~[l / nm] * fs / nm)
  { use this function in main file to define the GVD of a fiber: set_GVD_l("GVD_f(l)") }

; For plotting dispersion parameter in units of fs/nm/m from D_f[l] provided in fiber spectroscopic data:
DD_pm_f(l) :=
  if ld1 <= l / nm <= ld2 then
    (PM_f~[l / nm])


diagram 16: 

"PM Fiber GVD vs. Frequency"


x: 1000 , 1600 
"wavelength [nm]", @x
y: -40,30
" [dBm/m]", @y
frame
hx
hy

f: DD_pm_f(x*nm),
  color=olive, width=3, maxconnect=1, style=solid, "DCF_f"



