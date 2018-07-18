module HReg.Types
( LinRegValues(..)
, RegressionResult(..)
) where
import           Data.Time

-- | what we get from the measurements
data LinRegValues = LinRegValues
  { c        :: Double                     -- concentration of this series
  , kappa0   :: Double
  , kappaInf :: Double
  , kappaT   :: [(Maybe DiffTime, Double)] -- the optional time and corresping
                                           --   kappa_t
  }

-- | what we calculate by regression from the measurements
data RegressionResult = RegressionResult
  { intercept      :: Double
  , slope          :: Double
  , interceptError :: Double
  , slopeError     :: Double
  }

-- | make regression resuls nicely printable
instance Show RegressionResult where
  show a =
    "intercept: " ++ show (intercept a) ++ "±" ++ show (interceptError a) ++ "\n" ++
    "slope: " ++ show (slope a) ++ "±" ++ show (slopeError a)
