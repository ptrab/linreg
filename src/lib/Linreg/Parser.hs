{-# LANGUAGE OverloadedStrings #-}
module Linreg.Parser
( linRegValues
) where
import           Data.Attoparsec.Text.Lazy
import           Data.Time
import           Linreg.Types
import Prelude hiding (takeWhile)

-- | Parsing a single block with same "t" as shown in this repository (see "input" file).
linRegValues :: Parser [LinRegValues]
linRegValues = do
  -- Parse first line of a block, starting with a "t", followed by many times in seconds
  timeSeries <- do
    skipSpace
    _ <- string "t"
    _ <- takeWhile isHorizontalSpace
    times <- many1 $ do
      time <- double
      _ <- takeWhile isHorizontalSpace
      return time
    endOfLine
    return times
  -- Parse multiple lines of concentration series. First number is c0, then followed by
  -- (length times) concentrations corresponding to the times and then kappa0 and kappaInf
  linregSeries <- many1 $ do
    _ <- takeWhile isHorizontalSpace
    c' <- double
    _ <- takeWhile isHorizontalSpace
    kappaT' <- count (length timeSeries) $ do
      cAtTime <- double
      _ <- takeWhile isHorizontalSpace
      return cAtTime
    _ <- takeWhile isHorizontalSpace
    _ <- string "kappa0"
    _ <- takeWhile isHorizontalSpace
    _ <- string "="
    _ <- takeWhile isHorizontalSpace
    kappa0' <- double
    _ <- takeWhile isHorizontalSpace
    _ <- string "kappaInf"
    _ <- takeWhile isHorizontalSpace
    _ <- string "="
    _ <- takeWhile isHorizontalSpace
    kappaInf' <- double
    _ <- takeWhile isHorizontalSpace
    endOfLine
    return LinRegValues
      { c        = c'
      , kappa0   = kappa0'
      , kappaInf = kappaInf'
      , kappaT   = zip (map (picosecondsToDiffTime . (* 10^12) . round) $ timeSeries) kappaT'
      }
  return linregSeries
