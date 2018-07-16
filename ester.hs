-- meiner Meinung nach ist alles irgendwie hochgradig redundant

-- Berechnet c_0/c = (kappa_0 - kappa_inf) / (kappa_t - kappa_inf)
c0pc :: [Float] -> [Float]
c0pc c = let k0 = head (tail c)
             ki = head (tail (tail c))
             kt = tail (tail (tail c))
         in [ (k0 - ki) / (x - ki) | x <- kt ]

-- Liest die Messwerte ein und gibt direkt c_0/c aus
getNumbers :: String -> [[Float]]
getNumbers c = map c0pc (map (map readFloat) (map words (lines c)))
    where
        readFloat = read :: String -> Float

-- Liest die Konzentrationen ein, c
getConcentrations :: String -> [Float]
getConcentrations c = map head (map (map readFloat) (map words (lines c)))
    where
        readFloat = read :: String -> Float

-- sum_i^n sum_j^m y_ij
sumY :: [[Float]] -> Float
sumY y = sum (concat y)

-- sum_i^n sum_j^m c_j x_i
sumCX :: [Float] -> [Float] -> Float
sumCX c x = sum [ i * j | i <- x, j <- c ]

-- sum_i^n sum_j^m c_j^2 x_i^2
sumC2X2 :: [Float] -> [Float] -> Float
sumC2X2 c x = sum [ i^2 * j^2 | i <- x, j <- c ]

-- sum_i^n sum_j^m c_j x_i y_ij
sumCXY :: [Float] -> [Float] -> [[Float]] -> Float
sumCXY c x y = sum (zipWith(*) [ i * j | j <- c, i <- x ] (concat y))

-- Berechnet den y-Achsenabschnitt, 1
intercept :: [Float] -> [Float] -> [[Float]] -> Float
intercept c x y = ( sumCX c x * sumCXY c x y  - sumC2X2 c x * sumY y ) / ( (sumCX c x)^2 - m * n * sumC2X2 c x )
    where
        m = fromIntegral (length c)
        n = fromIntegral (length x)

-- Berechnet die Steigung, k
slope :: [Float] -> [Float] -> [[Float]] -> Float
slope c x y = ( sumCX c x * sumY y - m * n * sumCXY c x y) / ( (sumCX c x)^2 - m * n * sumC2X2 c x)
    where
        m = fromIntegral (length c)
        n = fromIntegral (length x)

-- Summe der Residuenquadrate, sum_i^n sum_j^m r_ij^2
sumResiduals2 :: [Float] -> [Float] -> [[Float]] -> Float
sumResiduals2 c x y = sum [ z^2 | z <- (zipWith(-) (concat y) [ a + b * i * j | j <- c, i <- x ]) ]
    where
        a = intercept c x y
        b = slope c x y

-- Reststreuung der Messwerte y, s
restS :: [Float] -> [Float] -> [[Float]] -> Float
restS c x y = sqrt ( (sumResiduals2 c x y) / ( n - 2 * m ) )
    where
        m = fromIntegral (length c)
        n = fromIntegral (length x)

-- Streuung des y-Achsenabschnitts, s(A)
interceptError :: [Float] -> [Float] -> [[Float]] -> Float
interceptError c x y = (restS c x y) * sqrt (sumC2X2 c x / ( n * m * sumC2X2 c x - (sumCX c x)^2 ))
    where
        m = fromIntegral (length c)
        n = fromIntegral (length x)

-- Streuung der Steigung, s(B)
slopeError :: [Float] -> [Float] -> [[Float]] -> Float
slopeError c x y = (restS c x y) * sqrt (n * m / ( n * m * sumC2X2 c x - (sumCX c x)^2 ))
      where
          m = fromIntegral (length c)
          n = fromIntegral (length x)

-- Messzeiten 30s bis 300s alle 30s, 360s bis 600s alle 60s
-- ... wuerde ich eigentlich ganz gerne noch mit einlesen
times :: [Float]
times = map (60*) ([ x / 2 | x <- [1 .. 10] ] ++ [6.0 .. 10.0])

-- Fuehrt die Rechnungen dann mal aus
calc :: FilePath -> IO()
calc fname = do
    contents <- readFile fname
    putStr "1: "
    print $ intercept (getConcentrations contents) times (getNumbers contents)
    putStr "s(1): "
    print $ interceptError (getConcentrations contents) times (getNumbers contents)
    putStrLn ""
    putStr "k [l mol^-1 s^-1]: "
    print $ slope (getConcentrations contents) times (getNumbers contents)
    putStr "s(k) [l mol^-1 s^-1]: "
    print $ slopeError (getConcentrations contents) times (getNumbers contents)

main :: IO()
main = do
    putStrLn "Which input file to process?"
    inputFile <- getLine
    putStrLn ""
    putStrLn "y_ij = 1 + c_j k x_i"
    putStrLn ""
    calc inputFile
