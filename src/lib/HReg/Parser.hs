module HReg.Parser
( parse_LinRegValues
) where
import           Data.Attoparsec.Text.Lazy
import qualified Data.Text                 as T
import qualified Data.Text.IO              as T
import           HReg.Types
import Control.Applicative
import Data.Time

-- | Parse the input file. It contains in a very simple format:
-- |   c (first value in row)
-- |   kappa0 (second value in row)
-- |   kappaInf (third value in row)
-- |   kappaT (many values, possibly associated with a time. If no time is
-- |     there, than it shall be a whitespace seperated double.
-- |     If time is there it shall follow the format "(Seconds, kappaT)"
parse_LinRegValues :: Parser LinRegValues
parse_LinRegValues = do
  c' <- double
  _ <- many1 $ char ' ' <|> char '\t'

  kappa0' <- double
  _ <- many1 $ char ' ' <|> char '\t'

  kappaInf' <- double
  _ <- many1 $ char ' ' <|> char '\t'

  kappaT' <- many1 parse_kappaT

  endOfLine

  return LinRegValues
    { c        = c'
    , kappa0   = kappa0'
    , kappaInf = kappaInf'
    , kappaT   = kappaT'
    }

  where
    parse_kappaT :: Parser (Maybe DiffTime, Double)
    parse_kappaT = do
      hasTime <- option Nothing (Just <$> char '(')
      case hasTime of
        Nothing -> do
          kappaT' <- double
          _ <- many' $ char ' ' <|> char '\t'
          return (Nothing, kappaT')
        Just x -> do
          _ <- many' $ char ' ' <|> char '\t'
          time' <- secondsToDiffTime <$> decimal
          _ <- many' $ char ' ' <|> char '\t'
          _ <- many' $ char ','
          _ <- many' $ char ' ' <|> char '\t'
          kappaT' <- double
          _ <- many' $ char ' ' <|> char '\t'
          _ <- char ')'
          _ <- many' $ char ' ' <|> char '\t'
          return (Just time', kappaT')
