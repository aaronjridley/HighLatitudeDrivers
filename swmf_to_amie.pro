
IonosphereFiles = findfile('./*.idl');
nFilesIono   = n_elements(IonosphereFiles)
IsIdlFiles = 1
iPot_ = 5
iEflux_ = 6
iAve_ = 7

if (nFilesIono eq 1 and strlen(IonosphereFiles(0)) lt 4) then begin
   IonosphereFiles = findfile('./*.tec') ;
   nFilesIono   = n_elements(IonosphereFiles)
   IsIdlFiles = 0
   iPot_ = 10
   iEflux_ = 7
   iAve_ = 8
endif

for iFile=0,nFilesIono-1 do begin

   IonoFile = IonosphereFiles(iFile)
   print, 'Reading file....       ',IonoFile

   if (IsIdlFiles) then begin
      iono_read_file, IonoFile, nVarsIono, nLatsIono, nLonsIono, $
                      IonoVars, IonoTime, IonoData
   endif else begin
      iono_read_tec_file, IonoFile, nVarsIono, nLatsIono, nLonsIono, $
                          IonoVars, IonoTime, IonoData
   endelse

   if (iFile eq 0) then begin

      nLats = nLatsIono
      nMlts = nLonsIono

      imf   = fltarr(nFilesIono,4)
      ae    = fltarr(nFilesIono,4)
      dst   = fltarr(nFilesIono,2)
      hp    = fltarr(nFilesIono,2)
      cpcp  = fltarr(nFilesIono)
      cpcpS = fltarr(nFilesIono)
      t     = dblarr(nFilesIono)

      Vars = ['Potential (V)', $
              'Energy Flux (ergs/cm2)','Mean Energy (eV)']

      if (IsIdlFiles) then begin
         lats  = 90.0 - reform(IonoData(0,0,0,*))
         latsS = 90.0 - reverse(reform(IonoData(1,0,0,*)) - 90.0)
         mlts = (reform(IonoData(0,1,*,0))/15.0 + 12.0) mod 24.0
      endif else begin
         lats  = 90.0 - reform(IonoData(0,3,0,*))
         latsS = 90.0 - reverse(reform(IonoData(1,3,0,*)) - 90.0)
         locN = where(lats  gt 40.0,nLats)
         locS = where(latsS gt 40.0,nLatsS)
         lats = lats(locN)
         latsS = latsS(locS)
         mlts = (reform(IonoData(0,4,*,0))/15.0 + 12.0) mod 24.0
      endelse

      nVars = n_elements(Vars)
      AmieData  = fltarr(nFilesIono, nVars, nmlts, nlats)
      AmieDataS = fltarr(nFilesIono, nVars, nmlts, nlats)

      c_r_to_a, itime, ionotime
      c_a_to_ymd, itime, ymd

   endif

   t(iFile) = ionotime

   c_r_to_a, itime, t(iFile)

   for iMLT = 0, nMlts-1 do begin
      i = iMLT 

      AmieData(iFile,0,iMLT,*) = IonoData(0,iPot_,i,locN)*1000.0   ; kV -> V
      AmieData(iFile,1,iMLT,*) = IonoData(0,iEflux_,i,locN)*1000.0 ; W/m2 -> ergs/cm2/s
      AmieData(iFile,2,iMLT,*) = IonoData(0,iAve_,i,locN)

      d = reverse(IonoData(1,iPot_,i,*))*1000.0
      AmieDataS(iFile,0,iMLT,*) = d(locS)
      d = reverse(IonoData(1,iEflux_,i,*))*1000.0
      AmieDataS(iFile,1,iMLT,*) = d(locS)
      d = reverse(IonoData(1,iAve_,i,*))
      AmieDataS(iFile,2,iMLT,*) = d(locS)

   endfor

   cpcp(iFile)  = max(AmieData (iFile,0,*,*)) - min(AmieData (iFile,0,*,*))
   cpcpS(iFile) = max(AmieDataS(iFile,0,*,*)) - min(AmieDataS(iFile,0,*,*))

endfor

Version = 1.0

FileOut = 'b'+ymd+'n.swmf'
amie_write_binary, fileout, Vars, lats, mlts, t, AmieData, $
                   imf = imf, ae = ae, dst = dst, hpi = hp, cpcp = cpcp, $
                   Version = Version

FileOut = 'b'+ymd+'s.swmf'
amie_write_binary, fileout, Vars, latsS, mlts, t, AmieDataS, $
                   imf = imf, ae = ae, dst = dst, hpi = hp, cpcp = cpcpS, $
                   Version = Version


end
