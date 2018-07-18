module HReg.Types
( LinRegValues(..)
) where
import           Data.Time

data LinRegValues = LinRegValues
  { c        :: Double                     -- concentration of this series
  , kappa0   :: Double
  , kappaInf :: Double
  , kappaT   :: [(Maybe DiffTime, Double)] -- the optional time and corresping
                                           --   kappa_t
  } deriving (Show)
