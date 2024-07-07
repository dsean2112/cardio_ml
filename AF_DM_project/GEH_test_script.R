# Load samples12, annotations12 from cardio_ml > data
sample <- 166

input12 <- samples12[sample,,]
ann12 <- annotations12[sample,,]

Rpeaks <- peak_isolation(input12[,1], ann12[,1], wave_value = 2)

XYZ <- kors(input12)
XYZ_M <- find_median_beat(XYZ, Rpeaks)

origin_point <- find_origin(XYZ_M, Rpeaks)

leads <- c(1:6,10:12)

GEH_QRS <- find_QRS_intervals(input12, ann12[,leads])
RToff <- find_RToff_interval(leads12 = input12, ann12 = ann12[,leads])
RTpeak <- find_RTpeak_interval(input12, ann12[,leads])
GEH_Rpeak <- find_geh_Rpeak(XYZ_M, origin)


geh(
  XYZ_M = XYZ_M,
  origin_point = origin_point,
  GEH_Ronset = GEH_Rpeak - GEH_QRS$QR,
  GEH_Rpeak = GEH_Rpeak,
  GEH_Roffset = GEH_Rpeak + GEH_QRS$RS,
  GEH_Tpeak = GEH_Rpeak + RTpeak,
  GEH_Toffset = GEH_Rpeak + RToff
)