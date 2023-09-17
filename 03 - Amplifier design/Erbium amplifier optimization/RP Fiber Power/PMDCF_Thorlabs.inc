;PM DCF from Thorlabs

ld1_dcf := 1290
ld2_dcf := 1650
dld_dcf := 20
defarray DCF_f[ld1,ld2,dld]

readlist l, DCF_f[l]:
1290.0, -66
1310.0, -68
1330.0, -70
1350.0, -72
1370.0, -74
1390.0, -76
1410.0,	-78
1430.0,	-80
1450.0, -82
1470.0, -84
1490.0, -85
1510.0,	-87
1530.0,	-94
1550.0,	-100
1570.0,	-107
1590.0,	-114
1610.0,	-121
1630.0, -122
1650.0,	-123

; Calculate dispersion in units of fs/nm/m from D_f[l] provided in fiber spectroscopic data:
GVD_DCF_f(l) :=
  if ld1_dcf <= l / nm <= ld2_dcf then
    -l^2 / (2pi * c)* (DCF_f~[l / nm] * fs / nm)
  { use this function in main file to define the GVD of a fiber: set_GVD_l("GVD_f(l)") }

; For plotting dispersion parameter in units of fs/nm/m from D_f[l] provided in fiber spectroscopic data:
DD_DCF_f(l) :=
  if ld1_dcf <= l / nm <= ld2_dcf then
    (DCF_f~[l / nm])


diagram 13: 

"DCF GVD vs. Frequency"


x: 1470 , 1650 
"wavelength [nm]", @x
y: -140,-70
" [dBm/m]", @y
frame
hx
hy

f: DD_DCF_f(x*nm),
  color=olive, width=3, maxconnect=1, style=solid, "DCF_f"