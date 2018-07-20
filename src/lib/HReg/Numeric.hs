module HReg.Numeric
( regression
) where
import           Data.List
import           Data.Maybe
import           Data.Time
import           HReg.Types

--------------------------------------------------------------------------------
-- helper functions
--------------------------------------------------------------------------------
-- | check if all elements of a list are identical
allTheSame :: (Eq a) => [a] -> Bool
allTheSame xs = all (== head xs) (tail xs)

-- | check if all elements of a list are "True"
allTrue :: [Bool] -> Bool
allTrue xs = all (== True) xs

map2Sum :: (Num a) => (a -> a) -> (a -> a -> a) -> [a] -> [a] -> a
map2Sum f g a b = sum [ (f i) `g` (f j) | i <- a, j <- b ]

-- | Default time generator for ith data point (0 based)
genT :: (Integral a, Floating b) => a -> b
genT i
  | i <= 9    = fromIntegral $ (i + 1) * 30
  | otherwise = fromIntegral $ (i - 4) * 60


--------------------------------------------------------------------------------
-- regression logic
--------------------------------------------------------------------------------
-- | calculate c0/c = (kappa0 - kappaInf) / (kappaT - kappaInf)
c0pc :: LinRegValues -> [(Maybe DiffTime, Double)]
c0pc vals =
  [ (fst i, (kappa0' - kappaInf') / (snd i - kappaInf'))
  | i <- kappaT'
  ]
  where
    kappa0' = kappa0 vals
    kappaInf' = kappaInf vals
    kappaT' = kappaT vals

-- | sum of all y values
-- |   sum_i^n sum_j^m y_ij
sumY :: (Num a) => [[a]] -> a
sumY y = sum . concat $ y

-- | Calculates the sum over the c*x products
-- |   sum_i^n sum_j^m c_j x_i
sumCX :: (Num a) => [a] -> [a] -> a
sumCX c x = sum $ (*) <$> c <*> x

-- | Sum of all squared c and x values
-- |   sum_i^n sum_j^m c_j^2 x_i^2
sumC2X2 :: (Num a) => [a] -> [a] -> a
sumC2X2 c x = sum $ (*) <$> map (^2) c <*> map (^2) x

-- | Zipping products
-- |   sum_i^n sum_j^m c_j x_i y_ij
sumCXY :: (Num a) => [a] -> [a] -> [[a]] -> a
sumCXY c x y = sum $ zipWith (*) ((*) <$> c <*> x) (concat y)

-- | Mean of a list
-- |   sum_i^n x_i / n
mean :: (Fractional a) => [a] -> a
mean x = (sum x) / (fromIntegral (length x))

-- | Do the regression
regression :: [LinRegValues] -> Maybe RegressionResult
regression vals
  | timesAllTheSame = Just RegressionResult
    { intercept = intercept'
    , slope = slope'
    , interceptError = interceptError'
    , slopeError = slopeError'
    , rSquared = rSquared'
    }
  | otherwise = Nothing
  where
    intercept' =
      (sCX * sCXY - sC2X2 * sY) /
      (sCX^2 - m * n * sC2X2)
    slope' =
      (sCX * sY - m * n * sCXY) /
      (sCX^2 - m * n * sC2X2)
    sumResiduals2' =
      sum . map (^2) $ z
    restS' =
      sqrt $ sR2 / (m * (n - 2))
    interceptError' =
      rS * sqrt (
      sC2X2 /
      (n * m * sC2X2 - sCX^2)
      )
    slopeError' =
      rS * sqrt (
      n * m /
      (n * m * sC2X2 - sCX^2)
      )
    rSquared' =
     1 - ssReg / ssTot
    --
    c' = map c vals
    mVals = map c0pc vals
    t' = nub . map (map $ fmap ((/10^12) . fromIntegral . diffTimeToPicoseconds) . fst) $ mVals
    timesAllTheSame = length t' == 1
    t'' = head t'
    t''' =
      [ fromMaybe (genT i) (t'' !! i)
      | i <- [0 .. length t'' - 1]
      ]
    x' = t'''
    y' = map (map snd) mVals
    sCX = sumCX c' x'
    sCXY = sumCXY c' x' y'
    sC2X2 = sumC2X2 c' x'
    sY = sumY y'
    m = fromIntegral . length $ c'
    n = fromIntegral . length $ x'
    n' = length $ x'
    a = intercept'
    b = slope'
    z = zipWith (-) (concat y') [ a + b * i * j | j <- c', i <- x' ]
    sR2 = sumResiduals2'
    rS = restS'
    ssReg = sumResiduals2'
    yBar = map mean y'
    ssTot = sum . map (^2) $ zipWith (-) (concat y') (concat [ replicate n' i | i <- yBar ])
