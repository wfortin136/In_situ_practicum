; Create time plot of Temperature Average per Rank

begin
wks = gsn_open_wks("x11","gsun01n")

num_time_steps = 240
num_ranks = 1792
in = new(num_time_steps, file)
data = new((/num_ranks+1, num_time_steps/), double)
data2 = new(num_time_steps, double)
xvals = new(num_time_steps, integer)
ranks = new(num_ranks+1, integer)
day = new(num_time_steps, integer)
minute = new(num_time_steps, integer)
xticks = new((/2,num_time_steps/), string)
thickness = new(num_ranks+1, double)

do step=0, num_time_steps-1
  val = step+1
  name = "T"+val+".nc"
  in(step) = addfile(name, "r")
  ;print(name)
  vNames = getfilevarnames(in(step)) 
  curfile = in(step)
  in(step) = addfile(name, "r")
  xvals(step) =step
  data2(step) = curfile->$vNames(3)$
  ;do step2 = 0, num_ranks-1  
    data(1:,step) = curfile->$vNames(0)$
  ;end do
end do

do step=0, num_time_steps-1
  ;print(in(step))
  
  ;print(data(step, 4))
  if((step*15)/60 .ne. 0) then
    xticks(0,step) = step
  end if
  xticks(1,step) = mod(step*15,60)
  ;print(day(step))
end do
;print(xticks)
do i=0, num_ranks
  ranks(i) = i+1
  thickness(i) = .1 
end do

thickness(num_ranks)=8.
data(0,:) = data2

;print(data2)
;print(data(89,0:239))

resources = True
resources@tiMainString = "Sub-Domain Average Temperature (K)"
resources@tiXAxisString          = "Simulation Time (Days)"
resources@tiYAxisString          = "Temperature (K)"
resources@trYMinF                = 210
resources@trYMaxF                = 260
resources@xyLineColors = ranks
resources@xyLineThicknesses   = thickness ;(/8., .1/)
resources@tmXBMode = "Explicit"
resources@tmXBValues = (/0, 48, 96, 144, 192, 240/)
resources@tmXBLabels = (/"0", "1", "2", "3", "4", "5"/)
resources@tmXBMinorValues = (/12, 24, 36, 60, 72, 84, 108, 120, 132, 144, 156, 168, 180, 204, 216, 228/)
resources@lgItemCount = 2
resources@lgLegendOn = True
resources@lgLabelStrings = (/"Global", "Ranks"/)
;restick = True
;restick@ttmFormat = "%d - %M"
;time_axis_labels(xticks, resources, restick)
plot = gsn_csm_xy(wks, xvals, data ,resources)

txid = new(2, graphic)
tx = True
tx@txFont        = "helvetica-bold" 
tx@txFontHeightF = 0.018
tx@txJust        = "CenterLeft"

txid(0) = gsn_add_text(wks, plot, "LEGEND", 120, 220.5, tx)
txid(1) = gsn_add_text(wks, plot, "Global Average", 138, 218, tx)
txid(1) = gsn_add_text(wks, plot, "Local Rank Averages", 138, 216, tx)

draw(plot)
frame (wks)
end
