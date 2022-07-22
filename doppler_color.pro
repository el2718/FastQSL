pro doppler_color, gamma=gamma, redvector=redvector, greenvector=greenvector, bluevector=bluevector, blackFirst=blackFirst
;+
; NAME:
;   DOPPLER_COLOR
;
; PURPOSE:
;   Generate a colortable from blue (0,0,255) via white (255,255,255) to red (255,0,0) 
;     for displying Doppler maps.
;
; KEYWORDS:
;   Gamma:  The value of gamma correction.  A value of 1.0 indicates a
;     linear ramp, i.e., no gamma correction.  Higher values of 
;     gamma give more contrast.  Values less than 1.0 yield lower contrast.
;
; MODIFICATION HISTORY:
;   R. Liu, October, 2006. 
;   Added gamma correction.  R. Liu, August 31, 2010.
;-

if keyword_set(gamma) then gamma=gamma else gamma=1.0
if keyword_set(blackFirst) then blackFirst=1B else blackFirst=0B
steps = 128
scaleFactor = (FINDGEN(steps) / (steps - 1))^gamma

; Do first 128 colors (blue to white)

  ; Red vector: 0 -> 255
redVector = 0 + (255 - 0) * scaleFactor; 

  ; Green vector: 0 -> 255
greenVector = 0 + (255 - 0) * scaleFactor

  ; Blue vector: 255 -> 255
if blackFirst eq 1B then blueVector = [0,REPLICATE(255, steps-1)] else blueVector = REPLICATE(255, steps)

; Do second 128 colors (white to red).

   ; Red vector: 255 -> 255
scaleFactor = (FINDGEN(steps) / (steps - 1))^(1./gamma)

redVector = [redVector, REPLICATE(255, steps)]

   ; Green vector: 255 -> 0
greenVector = [greenVector, 255 + (0 - 255) * scaleFactor]

   ; Blue vector: 255 -> 0
blueVector = [blueVector, 255 + (0 - 255) * scaleFactor]

TVLCT, redVector, greenVector, blueVector

;image=fltarr(255,50)
;data=findgen(255)
;for i=0,49 do image[*,i]=data
;tv,image
end
