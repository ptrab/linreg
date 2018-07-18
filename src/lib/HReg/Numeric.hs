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

-- | Do the regression
regression :: [LinRegValues] -> Maybe RegressionResult
regression vals
  | timesAllTheSame = Just RegressionResult
    { intercept = intercept'
    , slope = slope'
    , interceptError = interceptError'
    , slopeError = slopeError'
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
      sqrt $ sR2 / (n - 2 * m)
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
    a = intercept'
    b = slope'
    z = zipWith (-) (concat y') [ a + b * i * j | j <- c', i <- x' ]
    sR2 = sumResiduals2'
    rS = restS'

{-
-- | calculate the intercept of a set of fixed c0 measurements
intercept :: [LinRegValues] -> Maybe Double
intercept vals
  | timesAllTheSame = Just $
    (sCX * sCXY - sC2X2 * sY) /
    (sCX^2 - m * n * sC2X2)
  | otherwise = Nothing
  where
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

slope :: [LinRegValues] -> Maybe Double
slope vals
  | timesAllTheSame = Just $
    (sCX * sY - m * n * sCXY) /
    (sCX^2 - m * n * sC2X2)

  | otherwise = Nothing
  where
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

sumResiduals2 :: [LinRegValues] -> Maybe Double
sumResiduals2 vals
  | timesAllTheSame && isJust a && isJust b = Just $
    sum . map (^2) $ z
  | otherwise = Nothing
  where
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
    --
    a = intercept vals
    b = slope vals
    -- this is safe because of laziness and the check at top level
    a' = fromJust a
    b' = fromJust b
    z = zipWith (-) (concat y') [ a' + b' * i * j | j <- c', i <- x' ]

restS :: [LinRegValues] -> Maybe Double
restS vals
  | timesAllTheSame && isJust sR2 = Just $
    sqrt $
    sR2' /
    (n - 2 * m)
  | otherwise = Nothing
  where
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
    --
    a = intercept vals
    b = slope vals
      -- this is safe because of laziness and the check at top level
    a' = fromJust a
    b' = fromJust b
    z = zipWith (-) (concat y') [ a' + b' * i * j | j <- c', i <- x' ]
    --
    sR2 = sumResiduals2 vals
      -- this is safe because of laziness and the check at top level
    sR2' = fromJust sR2

interceptError :: [LinRegValues] -> Maybe Double
interceptError vals
  | timesAllTheSame && isJust rS = Just $
    rS' * sqrt (
    sC2X2 /
    (n * m * sC2X2 - sCX^2)
    )
  | otherwise = Nothing
  where
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
    --
    a = intercept vals
    b = slope vals
      -- this is safe because of laziness and the check at top level
    a' = fromJust a
    b' = fromJust b
    z = zipWith (-) (concat y') [ a' + b' * i * j | j <- c', i <- x' ]
    --
    sR2 = sumResiduals2 vals
      -- this is safe because of laziness and the check at top level
    sR2' = fromJust sR2
    --
    rS = restS vals
      -- this is safe because of laziness and the check at top level
    rS' = fromJust rS

slopeError :: [LinRegValues] -> Maybe Double
slopeError vals
  | timesAllTheSame && isJust rS = Just $
    rS' * sqrt (
    n * m /
    (n * m * sC2X2 - sCX^2)
    )
  | otherwise = Nothing
  where
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
    --
    a = intercept vals
    b = slope vals
      -- this is safe because of laziness and the check at top level
    a' = fromJust a
    b' = fromJust b
    z = zipWith (-) (concat y') [ a' + b' * i * j | j <- c', i <- x' ]
    --
    sR2 = sumResiduals2 vals
      -- this is safe because of laziness and the check at top level
    sR2' = fromJust sR2
    --
    rS = restS vals
      -- this is safe because of laziness and the check at top level
    rS' = fromJust rS
-}
