import           Data.Attoparsec.Text.Lazy
import qualified Data.Text.IO              as T
import           Linreg.Numeric
import           Linreg.Parser
import           System.Environment

main :: IO ()
main = do
  inputFilePath <- head <$> getArgs
  rawVals <- T.readFile inputFilePath
  let vals = parseOnly (many1 parse_LinRegValues) rawVals
      results = regression <$> vals
  case results of
    Right (Just x) -> putStr $ show x
    Right Nothing -> putStr "There are problems with your data. Maybe the times differ from measurement to measurement?"
    Left _ -> putStr "Could not parse your input"
