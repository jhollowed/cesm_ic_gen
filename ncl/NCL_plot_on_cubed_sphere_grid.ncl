
ad "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin
f  = addfile("/scratch/cjablono_root/cjablono1/cjablono/cesm2.1.3/se16L30/run/se16L30.cam.h1.0001-01-01-00000.nc","r")       ; Open netCDF file

field = "se_fvcubed_native"

output_format = "X11"
;output_format = "png"
;output_format = "eps"

time  = f->time
lat1d = f->lat
lon1d = f->lon
; lon1d@units = "degrees_east"
; lat1d@units = "degrees_north"
print("lat1d min/max = "+min(lat1d)+" "+max(lat1d))

level = 29
n=10  ; day 10
n=0  ; day 0
var = f->Q(:,:,:)
;var = f->PS(:,:)
;var = var / 100.
;var@units="hPa"

wks = gsn_open_wks(output_format,field)
print ("Writing file "+field+"."+output_format)

gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

res                     = True
res@gsnMaximize         = True         ; Maximize size of plot in frame
res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart = 24           ; at color index 24.

res@cnFillOn            = True         ; Turn on contour fill
;res@cnFillMode          = "AreaFill"   ; Style of fill. You can also use "CellFill" and "RasterFill"
res@cnLinesOn           = True        ; Turn on contour lines
res@cnLineLabelsOn      = False        ; Turn off contour line labels

res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.

res@sfXArray            = lon1d        ; Required to tell NCL where to
res@sfYArray            = lat1d        ; overlay data on globe.
; res@trGridType          = "TriangularMesh"
; res@gsnPolar            = "SH"         ; Southern Hemisphere for polar projection
; res@gsnPolar            = "NH"         ; Northern Hemisphere for polar projection
res@mpFillOn            = False      ; no map fill
res@mpOutlineOn         = False         ; no map
res@mpCenterLonF        = 180          ; plots from 0 - 360 longitudes

; currently automatic contour level spacing, manual spacing can be invoked
; res@cnLevelSelectionMode = "ManualLevels" ;-- set contour levels manually
; res@cnMinLevelValF = 940. ;-- minimum contour level
; res@cnMaxLevelValF = 1030. ;-- maximum contour level
; res@cnLevelSpacingF = 5 ;-- contour level spacing

  res@gsnDraw   = False                             ; don't draw yet
  res@gsnFrame  = False                             ; don't advance yet

print("plotting min/max " + min(var(:,:,:)) + " " + max(var(:,:,:)) )
;print("plotting min/max " + min(var(n,:)) + " " + max(var(n,:)) )
res@tiMainString        = field+" time="+time(n)+" days"    ;
plot = gsn_csm_contour_map_ce(wks,var(n,level,:),res)     ; Draw a contour plot on a cylindrical equidistant map
;plot = gsn_csm_contour_map_ce(wks,var(n,:),res)     ; Draw a contour plot on a cylindrical equidistant map
;plot = gsn_csm_contour_map_polar(wks,var(n,:),res)     ; Draw a contour plot polar stereographic map

draw(plot)
frame(wks)

end

