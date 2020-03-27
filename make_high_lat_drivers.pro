subaurora = 1.0
variability = 0.0
cc = 0.0

ByBase =   0.0
BzBase =  -2.0
;BzPert =  -5.0
BzPert =   0.0

BzChangeTime = 15.0 * 60.0

DeltaOCFLBOnset    =  -1.0
DeltaOCFLBRecovery =  10.0
;DeltaOCFLBGrowth   =   5.0
DeltaOCFLBGrowth   =   0.0

rBz = 0.0
rBy = 0.0
rH  = 0.0

filename = 'gravity_waves_v6.bin'

; The way that the AMIE files work is that you can have a 
; non-uniform distribution of times.  We want nPrepDays worth
; constant potential patterns, so we will write one pattern out
; at the very beginning, then skip nPrepDays and start the rest of
; the patterns.  I.e., there will be one pattern, then a large gap,
; then data every minute for a day.

nPrepDays = 3

nHoursToRun = 3.0
HourStart   = 11.0

nTimes = 60*nHoursToRun + 2               ; 1440+1

iTime = [2000,3,18,0,0,0]
c_a_to_r, iTime, basetime

time = dindgen(nTimes-2)/(nTimes-2) * nHoursToRun*3600.0D + $
       basetime + HourStart * 3600.0D + $
       (nPrepDays)*24.0D*3600.0D

endtime = basetime + (nPrepDays+1)*24.0D*3600.0D
time = [basetime,time,endtime]

imf = fltarr(nTimes, 4)
ae = fltarr(nTimes, 4)
dst = fltarr(nTimes, 2)
hpi = fltarr(nTimes, 2)
cpcp = fltarr(nTimes)
Emax = fltarr(nTimes)
Vmax = fltarr(nTimes)
AU = fltarr(nTimes)
AH = fltarr(nTimes)
AL = fltarr(nTimes)

; bx has no meaning, really.  Can be ignored

bx = fltarr(nTimes)+0.0 + 0.0*randomn(s,nTimes)

; by sets a circular potential at the pole.  You could set this to 
; zero if you want a dawn-dusk symmetric potential pattern.

by = fltarr(nTimes) + ByBase + rBy * randomn(s,nTimes)

; This sets the strength of the background two cell convection pattern

bz = fltarr(nTimes)+BzBase + rBz*randomu(s,nTimes)

; This sets the resolution:

dlat = 0.5
dmlt = 0.25
lowlat = 45.0
nLats = (90.0-lowlat)/dlat + 1
nMlts = 24.0/dmlt + 1

lats = findgen(nLats)*dlat + lowlat
mlts = findgen(nMlts)*dmlt 

; This sets where the potential and aurora maximizes (potential has
; max/min at the ocflb, while the aurora is equatorward of this.)

ocflbBase = 70.0 - 3.0*cos(mlts*!pi/12.0)

Ef = dblarr(nMlts, nLats)
El = dblarr(nMlts, nLats)
Em = dblarr(nMlts, nLats)

Vexb = dblarr(nMlts, nLats, 4)

; -----calculate the dipole field----- 
Br = dblarr(nMlts, nLats)
Bl = dblarr(nMlts, nLats)
Bm = dblarr(nMlts, nLats)
Mb = double(8.0e22)
mudiv = double(1.e-7)
r0 = double(6372+110)*1000.0
const = mudiv*Mb/r0^3
for iMlt = 0, nMlts-1 do begin
   for iLat = 0, nLats-1 do begin
      Br[iMlt, iLat] = -2*const*sin(lats[iLat]/180*!pi)
      Bl[iMlt, iLat] = const*cos(lats[iLat]/180.*!pi)
      Bm[iMlt, iLat] = const*sqrt(1+3.*sin(lats[iLat]/180.*!pi)^2)
   endfor
endfor

;-------------------------------------

;nVars = 3
;nVars = 4
nVars = 6

data = fltarr(nTimes, nVars, nMlts, nLats)

potential = fltarr(nMlts, nLats)
eflux     = fltarr(nMlts, nLats)
avee      = fltarr(nMlts, nLats)
; area is used to calculate the hemispheric power
area      = fltarr(nMlts, nLats)
for iLat = 0,nLats-1 do $
   area(*,iLat) = 111.0*(111.0*360.0/nMlts) * cos(lats(iLat)*!dtor)*1000.0*1000.0

; I have a substorm going off on the second day, starting at 10 UT.
; You can "turn this off" by setting the time after the run is
; finished (i.e. > 48 hours) (if the time is > 100, that is just to
; turn that feature off!)

SubstormGrowth   = (24.0 * nPrepDays + 12.0)*3600.0D + basetime
SubstormOnset    = (24.0 * nPrepDays + 13.0)*3600.0D + basetime
SubstormRecovery = (24.0 * nPrepDays + 113.5)*3600.0D + basetime
DtRecovery = 3600.0D

c_r_to_a, itime, SubstormGrowth
print, 'growth starts at : ', itime
c_r_to_a, itime, SubstormOnset
print, 'onset starts at : ', itime
c_r_to_a, itime, SubstormRecovery
print, 'Recovery starts at : ', itime

; This is stuff for the substorm

for i=0,nTimes-1 do begin

   ocflb = ocflbBase

   potential_width0 = 10.0
   potential_width = potential_width0
;   fac_pt_base = 5.0

   fac_ef_base = 1.0
   fac_ae_base = 6.0

   ; This is the amplitude of the potential in the harang
   ; discontinuity.  I would set this to zero, or you could play with
   ; it for a strong jet....

   amp_h = 7.5 + randomn(s,1) * rH
   amp_h = amp_h(0)

   harang_mlt_peak = 1.0
   harang_mlt_width = 2.0

   ; This is the amplitude of the saps potential:
   saps_efield = 50.0 ; mV/m
   saps_width = 2.0; deg * 120 km/deg
   saps_potential_max = -saps_efield*saps_width*120.0/1000.0 ; kV
   saps_peak_offset = -7.0 ; deg from the OCFLB
   saps_mlt_peak = 23.0
   saps_mlt_width = 4.0

   ; This is for the auroral energy flux for the substorm:
   ss_peak_mlt = 23.0
   ss_mlt_width = 2.5
   amp_ef_ss = 0.0

   diffuse_mlt_width = 6.0

   ss_equator_mlt = 1.0
   ss_equator_width = 2.0
   ss_equator_amp = 0.0

   ; These are for some localized features to give us the correct 
   ; FAC pattern:

   ss_upper_lat_width = 2.0
   ss_upper_pot = 0.0 ; changes through SS
   ss_upper_mlt_peak = 1.0
   ss_upper_mlt_width = 1.0

   ss_lower_lat_width = 2.0
   ss_lower_pot = 0.0 ; changes through SS
   ss_lower_mlt_peak = 23.0
   ss_lower_mlt_width = 1.0

   ; This is the amplitude of the potential in the PBIs
   ; I would set this to zero, or you could play with
   ; it for a strong jet....

   amp_pbi = 30.0 + randomn(s,1) * rH
   amp_pbi = amp_pbi(0)*1000.0

   ; This is the amplitude, multiplied by Bz for the strength of the CPCP
   amp = 7.5*1000.0
   ; Amplitude of the E-flux
   amp_ef = 5.0                 ; /1000.0
   ; Amplitude of the average energy
   amp_ae = 3.0

   pbi = 0

   saps_pot = 0.0

   ; Substorm growth - expand the potential pattern
   if (time(i) gt SubstormGrowth and $
       time(i) le SubstormOnset) then begin
      Dt = (time(i)-SubstormGrowth)/(SubstormOnset-SubstormGrowth)
      ocflb = ocflbBase + Dt*DeltaOCFLBGrowth
      ; Change Bz over 15 minute
      Dt = (time(i)-SubstormGrowth)/BzChangeTime
      if (dt gt 1.0) then Dt = 1.0
      bz(i) = bz(i) + BzPert*Dt
  
      ;Dt = (time(i)-SubstormGrowth)/3600.0
      ;amp_h = amp_h + 10.0 * dt

      ; Let's do something with shielding too...

   endif

   if (time(i) gt SubstormOnset) then begin

      ocflb = ocflbBase + DeltaOCFLBGrowth
      bz(i) = bz(i) + BzPert
      
   endif
   
   print, 'bz : ', bz(i), mm(ocflb)

   
;   ; Onset - retract the auroral oval, and broaden it greatly and
;   ; strengthen the precipitation
;
;   if (time(i) gt SubstormOnset and $
;       time(i) le SubstormRecovery) then begin
;      Dt = (time(i)-SubstormOnset)/(SubstormRecovery-SubstormOnset)
;
;      ocflb = ocflbBase + DeltaOCFLBGrowth + $
;              Dt*DeltaOCFLBOnset*cos(mlts*!pi/12.0)
;      ; release flow shear
;      ;amp_h = amp_h + 10.0* (1-dt)
;      
;      fac_ef_base = 1.0 ;+ 2.0 * dt
;      amp_ef = amp_ef ;+ (8.0-amp_ef)*dt
;      amp_ef_ss = (15.0-amp_ef)*dt
;      ss_equator_amp = 15.0*dt
;
;      amp_ae = amp_ae + (amp_ae*2-amp_ae)*dt
;;      amp_ae = amp_ae + (2.5-amp_ae)*dt
;      potential_width = potential_width0 + 5.0*dt
;;      fac_pt_base = 5.0 + 10.0*dt
;      ; Change Bz over 10 minutes
;      Dt = (time(i)-SubstormOnset)/1800.0
;      if (dt gt 1.0) then Dt = 1.0
;      bz(i) = bz(i) - 2.0*(1.0-Dt)
;
;      ; Change amp_sub over 15 minutes
;      Dt = (time(i)-SubstormOnset)/900.0
;      if (Dt gt 1.0) then Dt = 1.0
;      ; add by me
;      ss_upper_pot = 10.0*1000.0 * dt
;      ss_lower_pot = -10.0*1000.0 * dt
;
;      ; SAPS
;      saps_pot = Dt * saps_potential_max
;
;   endif
;
;   ; recovery - go back to the original state
;
;   if (time(i) gt SubstormRecovery and $
;       time(i) le SubstormRecovery+DtRecovery) then begin
;      Dt = (time(i)-SubstormRecovery)/(DtRecovery)
;      ocflb = ocflbBase + (1.0-Dt)*(DeltaOCFLBGrowth + $
;              DeltaOCFLBOnset*cos(mlts*!pi/12.0))
;
;      fac_ef_base = 1.0 ;+ 2.0 * (1.0-dt)
;      amp_ef = amp_ef ;+ (8.-amp_ef)*(1.0-dt)
;      amp_ef_ss = (15.0-amp_ef)*(1.0-dt)
;      ss_equator_amp = 15.0*(1.0-dt)
;      amp_ae = amp_ae + (2*amp_ae-amp_ae)*(1.0-dt)
;      potential_width = potential_width0 + 5.0*(1.0-dt)
;;      fac_pt_base = 5.0 + 10.0*(1.0-dt)
;
;      ss_upper_pot = 10.0*1000.0 * (1.0-dt)
;      ss_lower_pot = -10.0*1000.0 * (1.0-dt)
;
;      ; SAPS
;      saps_pot = (1.0-Dt) * saps_potential_max
;
;   endif

   amp_h = amp_h * 1000.0

   m = (90.0-mean(ocflb))/2

   for iMlt = 0,nMlts-1 do begin
      ;fac = fltarr(nLats) + fac_pt_base/2
      ;add by me increase the potential width
      width_pot = fltarr(nLats) + potential_width ; fac_pt_base/3.*2.

      fac_ef = fltarr(nLats) + fac_ef_base * (0.3+1.7*(cos(mlts(iMlt)*!pi/12.0)+1.0))/3.0
      fac_ae = fltarr(nLats) + fac_ae_base

      d = lats-ocflb(iMlt)

      ; build the potential pattern

      l = where(d ge 0)
      width_pot(l) = (90.0-ocflb(iMlt))/3.0 ; fac_pt_base

      line = exp(-abs(lats-ocflb(iMlt))/width_pot)
      dawn_dusk = 1.0-sin(mlts(iMlt)*!pi/12.0)*0.1
      potential(iMlt,*) = -bz(i) * amp * line * sin(mlts(iMlt)*!pi/12.0)*dawn_dusk

      ; Some By
      line = -amp/2 * by(i) * exp(-(90.0-lats)/m)
      potential(iMlt,*) = potential(iMlt,*) + line

      ; Harang Discontinuity 

      dMlt = acos(cos((mlts(iMlt)-harang_mlt_peak)*!pi/12.0))*12.0/!pi
      dMlt = exp(-dMlt/harang_mlt_width)

      line = -amp_h * exp(-2*(lats-ocflb(iMlt))^2/width_pot^2) * dMlt
      potential(iMlt,*) = potential(iMlt,*) + line
;      potential(iMlt,*) = line

      ; Add some SAPS:
      d = lats-(ocflb(iMlt)+saps_peak_offset)
      l = where(d ge 0)
      width_pot(l) = (90.0-(ocflb(iMlt)+saps_peak_offset))/4.0 ; fac_pt_base
      l = where(d lt 0)
      width_pot(l) = saps_width
      dMlt = acos(cos((mlts(iMlt)-saps_mlt_peak)*!pi/12.0))*12.0/!pi
      dMlt = exp(-dMlt/saps_mlt_width)

      line = saps_pot * exp(-abs(d)/width_pot) * dMlt
;      l = where(d lt 0)
;      line(l) = saps_pot * exp(-d(l)^2/width_pot^2) * dMlt
      potential(iMlt,*) = potential(iMlt,*) + line*1000.0
;      potential(iMlt,*) = line*1000.0

      if (pbi) then begin
         ; PBIs - centered at midnight
         line = -intensity * amp_pbi * exp(-(lats-ocflb_pbi(iMlt))^2/1.5^2) * $
                exp(50.0*((cos((mlts(iMlt)-loctim)*!pi/12.0)-1.0)/2.0))
         l = where(line gt 0,c)
         if (c gt 0) then line(l) = 0.0
        potential(iMlt,*) = potential(iMlt,*) + line
;         potential(iMlt,*) = line
      endif

      ; Add some variability
;      fac(l) = 5.0
      vari_line  = $
         variability*randomn(s,nLats) * exp(-abs(lats-ocflb(iMlt))/width_pot)
      vari_line2 = $
         variability*randomn(s,nLats) * exp(-abs(lats-ocflb(iMlt))/width_pot)

      potential(iMlt,*) = potential(iMlt,*) + $
                          vari_line*potential(iMlt,*)
   
      if (subaurora) then begin

;         dMlt = acos(cos((mlts(iMlt)-ss_lower_mlt_peak)*!pi/12.0))*12.0/!pi
;         dMlt = exp(-dMlt/ss_lower_mlt_width)
;         xline = ocflb(iMlt)-3  
;         sublower = ss_lower_pot*exp(-abs(lats-xline)/ss_lower_lat_width) *dMlt
;         potential(iMlt,*) = potential(iMlt,*) + sublower

         dMlt = acos(cos((mlts(iMlt)-ss_upper_mlt_peak)*!pi/12.0))*12.0/!pi
         dMlt = exp(-dMlt/ss_upper_mlt_width)
         xline = ocflb(iMlt)+2  
         subupper = ss_upper_pot*exp(-abs(lats-xline)/ss_upper_lat_width) *dMlt
         potential(iMlt,*) = potential(iMlt,*) + subupper

      endif

      ; build the energy flux pattern

;      if (mlts[iMlt] ge 12) then begin
;         offset = 2.*cos((mlts[iMlt]-18)/12*!pi)
;      endif else begin
;         offset = 2.*cos((mlts[iMlt]-6)/12*!pi)
;      endelse
      offset = 2.0
      d = lats-ocflb(iMlt)+offset
      l = where(d gt 0)
      fac_ef(l) = fac_ef(l)*1.25

      line = exp(-abs(lats-ocflb(iMlt)+offset)^2/fac_ef^2)

      dMlt = acos(cos((mlts(iMlt))*!pi/12.0))*12.0/!pi
      dMlt = exp(-dMlt/diffuse_mlt_width)

      ;;l = where(line gt 0.5,c)
      ;if (c gt 0) then line(l) = 1.0
      eflux(iMlt,*) = amp_ef * line * dMlt
;      eflux(iMlt,*) = eflux(iMlt,*) + $
;                      cc * vari_line * eflux(iMlt,*) + $
;                      sign(cc) * (1.0 - abs(cc)) * vari_line2 * eflux(iMlt,*)

      ; Add substorm enhancement:

      dMlt = acos(cos((mlts(iMlt)-ss_peak_mlt)*!pi/12.0))*12.0/!pi
      dMlt = exp(-dMlt/ss_mlt_width)

;      fac_ef = fac_ef*0.75

;      offset = 2.0
;      d = lats-ocflb(iMlt)+offset
;      l = where(d gt 0)
;      fac_ef(l) = fac_ef(l)*2.0
;      line = exp(-abs(lats-ocflb(iMlt)+offset)^2/fac_ef^2)
      eflux(iMlt,*) = eflux(iMlt,*) + amp_ef_ss * line * dMlt

      ; Equatorward edge arc
      dMlt = acos(cos((mlts(iMlt)-ss_equator_mlt)*!pi/12.0))*12.0/!pi
      dMlt = exp(-dMlt/ss_equator_width)
      offset = saps_peak_offset + 3.0
      line = exp(-abs(lats-ocflb(iMlt)-offset)^2/1^2)
      eflux(iMlt,*) = eflux(iMlt,*) + ss_equator_amp * line * dMlt

      if (pbi) then begin
         ; PBIs
         line = intensity * 10.0 * exp(-(lats-ocflb_pbi(iMlt))^2/1.5^2) * $
                exp(50.0*((cos((mlts(iMlt)-loctim)*!pi/12.0)-1.0)/2.0))
         l = where(line lt 0,c)
         if (c gt 0) then line(l) = 0.0
         eflux(iMlt,*) = eflux(iMlt,*) + line
      endif


;      if (subaurora) then begin
;         amp_left = 0.0
;         amp_right = 0.0
;         mltdecay = 1.0
;         width = 1.
;
;         ;; if (mlts(iMlt) le 1 and mlts(iMlt) ge 0) then amp_right = amp_dec                               
;         ;; if (mlts(iMlt) le 24. and mlts(iMlt) ge 19.) then amp_right = amp_dec * exp(-abs(mlts(iMlt)-24)/mltdecay)
;         ;; if (mlts(iMlt) ge 1. and mlts(iMlt) le 6) then amp_right = amp_dec * exp(-abs(mlts(iMlt)-1)/mltdecay)     
;         ;; xline = ocflb(iMlt) + 2.
;         ;; subup = amp_right*exp(-abs(lats-xline)/width)
;
;
;         if (mlts(iMlt) le 1 and mlts(iMlt) ge 0) then amp_right = amp_dec
;         if (mlts(iMlt) le 24 and mlts(iMlt) ge 23) then amp_right = amp_dec
;         if (mlts(iMlt) le 23. and mlts(iMlt) ge 18.) then amp_right = amp_dec * exp(-abs(mlts(iMlt)-23)/mltdecay)
;         if (mlts(iMlt) ge 1. and mlts(iMlt) le 6) then amp_right = amp_dec * exp(-abs(mlts(iMlt)-1)/mltdecay)     
;         xline = ocflb(iMlt)-1
;         subup = amp_right*exp(-abs(lats-xline)/width)
;         eflux(iMlt,*) = eflux(iMlt,*) + subup
;
;      endif


      ; build the averaged energy pattern

      line = exp(-(lats-ocflb(iMlt))^2/fac_ae^2)
      avee(iMlt,*) = amp_ae * line * (cos(mlts(iMlt)*!pi/12.0)+5.0)/6
      avee(iMlt,*) = avee(iMlt,*) + $
                     cc * vari_line * avee(iMlt,*) + $
                     sign(cc) * (1.0 - abs(cc)) * vari_line2 * avee(iMlt,*)

      if (pbi) then begin
         ; PBIs
         line = intensity * 3.0 * exp(-(lats-ocflb_pbi(iMlt))^2/1.5^2) * $
                exp(50.0*((cos((mlts(iMlt)-loctim)*!pi/12.0)-1.0)/2.0))
         l = where(line lt 0,c)
         if (c gt 0) then line(l) = 0.0
         avee(iMlt,*) = avee(iMlt,*) + line
      endif

   endfor

   ; set base values for aurora

   l = where(avee lt 0.5)
   avee(l) = 0.5

   l = where(eflux lt 0.01)
   eflux(l) = 0.01

   print, 'eflux : ', mm(eflux)
   
   
   Df = 2*!pi/(nMlts-1)
   ;;;;;; solve for electric field
   for iMlt = 0, nMlts-1 do begin
      for iLat = 1, nLats-2 do begin
         pre = iMlt-1
         aft = iMlt+1
         if (pre < 0) then pre = nMlts-1
         if (aft eq nMlts) then aft = 0
         Ef(iMlt,iLat) = -1./r0/(cos(lats(ilat)/180.*!pi)+0.01)*(potential(aft,iLat)-potential(pre,iLat))/2./Df
         pre = iLat-1
         aft = iLat+1
         El(iMlt,iLat) = -1./r0*(potential(iMlt,aft)-potential(iMlt,pre))/2./((lats(iLat+1)-lats(iLat-1))/180.*!pi)
         Em(iMlt,iLat) = sqrt(Ef(iMlt,iLat)^2 + El(iMlt,iLat)^2)
      endfor
   endfor


   Vexb(*,*,0) = -Ef*Bl/Bm^2
   Vexb(*,*,1) = Ef*Br/Bm^2
   Vexb(*,*,2) = -Br*El/Bm^2
   Vexb(*,*,3) = sqrt(Vexb(*,*,1)^2 + Vexb(*,*,2)^2)
   
   Sig_p = 40 * avee/(16.+avee^2)*eflux^0.5
   Sig_h = 0.45*avee^0.85*Sig_p
   
   
   Jp = dblarr(nMlts,nLats,3)
   Jp(*,*,0) = 0
   Jp(*,*,1) = Sig_p*El
   Jp(*,*,2) = Sig_p*Ef

   Jh = dblarr(nMlts,nLats,3)
   Jh(*,*,0) = 0 
   Jh(*,*,1) = -Sig_h*Ef
   Jh(*,*,2) = Sig_h*El

   
   Bh = dblarr(nMlts,nLats,4)
   for iMlt = 0, nMlts-1 do begin
      for iLat = 1, nLats-2 do begin
         Bh(iMlt,iLat,0) = 0
         Bh(iMlt,iLat,1) = Jh(iMlt,iLat,2)*4.*!pi*1.e-7/2.
         Bh(iMlt,iLat,2) = -Jh(iMlt,iLat,1)*4*!pi*1.e-7/2.
         Bh(iMlt,iLat,3) = sqrt(Bh(iMlt,iLat,1)^2 + Bh(iMlt,iLat,2)^2)
      endfor
   endfor

   AU(i) = max(Bh(*,*,1))*1.e9
   AL(i) = min(Bh(*,*,1))*1.e9
   AE(i) = AU(i) - AL(i)

   tmp = (reform(Jp(*,*,0))*Br(*,*) + Jp(*,*,1)*Bl(*,*))/Bm^2
   Jp_perp = dblarr(nMlts, nLats, 3)
   Jp_perp(*,*,0) = Jp(*,*,0) - tmp * Br
   Jp_perp(*,*,1) = Jp(*,*,1) - tmp * Bl
   Jp_perp(*,*,2) = Jp(*,*,2)

   Jp_para = dblarr(nMlts, nLats)

   for iMlt = 0, nMlts-1 do begin
      for iLat = 1, nLats-2 do begin
         theta = lats(iLat)/180.*!pi
         pre = lats(iLat-1)/180.*!pi
         aft = lats(iLat+1)/180.*!pi
         Jp_para(iMlt,iLat) = 1./(cos(theta)+0.01)/r0 * (Jp_perp(iMlt,iLat+1)*(cos(aft)+0.01) - Jp_perp(iMlt,iLat-1)*(cos(pre)+0.01))/(aft-pre) 
         pre = iMlt-1
         if (pre < 0) then pre = nMlts-1
         aft = iMlt+1
         if (aft eq nMlts) then aft = 0
         Jp_para(iMlt,iLat) += 1./(cos(theta)+0.01)/r0 * (Jp_perp(aft,iLat)-Jp_perp(pre,iLat))/2./Df 
         Jp_para(iMlt,iLat) = -float(Jp_para(iMlt,iLat))
      endfor
   endfor
   
   ;print, Jp_para(3,4)
   ; fill the data array
   data(i,0,*,*) = potential
   data(i,1,*,*) = eflux
   data(i,2,*,*) = avee
   data(i,3,*,*) = Jp_para
   data(i,4,*,*) = Ef*1000.0
   data(i,5,*,*) = El*1000.0


   ;calculate the max electric field
   Emax(i) = float(max(Em))

   ;calculate the max ExB drift
   Vmax(i) = float(max(Vexb(*,*,3)))

   ; calculate hemispheric power

   hpi(i,0) = total(area*eflux*0.001)/1.0e9

   ; calculate cpcp

   cpcp(i) = max(potential)-min(potential)

   print, 'hp: ', i, hpi(i,0), cpcp(i), mm(ef)*1000.0

endfor

imf(*,0) = 400.0
imf(*,1) = bx
imf(*,2) = by
imf(*,3) = bz

Version = 1.0

;vars = ['Potential (kV)', 'Total Energy Flux (ergs/cm2/s)', 'Mean Energy (ergs)', 'J Parallel (J/m3/s)']
vars = ['Potential (kV)', 'Total Energy Flux (ergs/cm2/s)', 'Mean Energy (ergs)', 'J Parallel (J/m3/s)',$
       'Eeast (mV/m)','Enorth (mV/m)']
;vars = ['Potential (kV)', 'Total Energy Flux (ergs/cm2/s)', 'Mean Energy (ergs)']
; write amie file

amie_write_binary, filename, Vars, lats, mlts, time, data, $
                   imf = imf, ae = ae, dst = dst, hpi = hpi, cpcp = cpcp, $
                   Version = Version



ppp = 5

setdevice, 'hp.ps', 'p', 5, 0.9
space =0.02
pos_space,ppp, space, sizes, nx=1
get_position, ppp, space, sizes, 0, pos,/rect
plot,(time-basetime)/3600.-24*nPrepDays, hpi[*,0], yrange=[0,100], ytitle='HP (GW)', xtickname=replicate(' ' ,10), xstyle=1, ystyle=1, pos=pos, xrange=[9.5,13.5], /noerase

get_position, ppp, space, sizes, 1, pos, /rect
plot,(time-basetime)/3600.-24*nPrepDays, cpcp*1.e-3, yrange=[0,100], ytitle='CPCP (kV)', xstyle=1, ystyle=1, pos=pos, xrange=[9.5,13.5], /noerase, xtickname=replicate(' ' ,10)

get_position, ppp, space, sizes, 2, pos, /rect
plot,(time-basetime)/3600.-24*nPrepDays, Emax*1000, yrange=[0,100], ytitle='E max (m V/m)', xstyle=1, ystyle=1, pos=pos, xrange=[9.5,13.5], /noerase, xtickname=replicate(' ' ,10)

get_position, ppp, space, sizes, 3, pos, /rect
plot,(time-basetime)/3600.-24*nPrepDays, Vmax, yrange=[0,1000], ytitle='Vexb max (m/s)', xstyle=1, ystyle=1, pos=pos, xrange=[9.5,13.5], /noerase, xtickname=replicate(' ' ,10)

get_position, ppp, space, sizes, 4, pos, /rect
plot,(time-basetime)/3600.-24*nPrepDays, AE, yrange=[0,300], ytitle='AE (nT)', xtitle='Simulation Hours', xstyle=1, ystyle=1, pos=pos, xrange=[9.5,13.5], /noerase

closedevice

end

