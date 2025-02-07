library(tidyverse)
library(readxl)

ps_adult <- read_xlsx(
  "data/adult_chinook_pcbs.xlsx",
  sheet = "Data",
  # range = "A1:AY156",
  na = "NA",
  col_types = c(
    "numeric", # Record
    "text", # Matrix
    "skip", # Chemistry (Y/N)
    "skip", # Used in Paper (Y/N)
    "text", # SizeClass
    "text", # SampleID
    "numeric", # FishID
    "numeric", # MarineArea
    "text", # CollectionEvent
    "skip", # Year
    "skip", # CollectionMonth
    "date", # CollectionDate
    "skip", # CollectionSeasonYear
    "numeric", # MonthsAtSea
    "skip", # CollectionSeason
    "skip", # ForkLength_mm
    "numeric", # ForkLength_cm
    "skip", # CalcTL_cm
    "text", # ScaleAge_GilbertRich
    "numeric", # SWAge (yrs)
    "numeric", # FWAge
    "numeric", # Sex
    "text", # OutMigrationLH
    "text", # Origin
    "numeric", # PrcntSolids
    "numeric", # PrcntLipids
    "numeric", # TotalPCBs ng/g ww
    "numeric", # Sum11PBDEng/g ww
    "skip", # PBDE47 ng/g ww
    "skip", # PBDE99 ng/g ww
    "skip", # PBDE99rND ng/g ww
    "skip", # Sum47&49zND ng/g ww
    "skip", # Sum47&49rND ng/g ww
    "skip", # Stock_AggLevel5_BestEstimate
    "skip", # Prob....35
    "skip", # Second Best Estimate
    "skip", # Prob....37
    "skip", # Third Best Estimate
    "skip", # Probability
    "skip", # 1- Best
    "skip", # Odds Ratio
    "skip", # CWT releasedate
    "skip", # AdClip_yn
    "skip", # CW_yn
    "skip", # CWT_TagCode
    "skip", # CWT_BroodYear
    "skip", # CWT_FirstReleaseDate
    "skip", # CWT_LastReleaseDate
    "skip", # CWT_ReleaseLocation
    "skip", # CWT_HatcheryLocation
    "skip" # CWT_StockLocation
  )
)
