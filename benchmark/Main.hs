{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
module Main (main) where
import Criterion.Main
import PWPs.SimplePolynomials as SP

longPoly :: (Integral b, Floating a) => b -> Poly a
longPoly m = Poly (replicate (2^m) pi)

mulLongPolys :: Int -> Poly Double
mulLongPolys n = longPoly n * longPoly n

addLongPolys :: Int -> Poly Double
addLongPolys n = longPoly n + longPoly n

convPolys :: Int -> [(Double, Poly Double)]
convPolys n = convolvePolys (0, 1, Poly [1]) (0, 1, longPoly n)

main :: IO ()
main = defaultMain [
  bgroup "con" [ bench "a1"  $ nf addLongPolys 1
               , bench "a5"  $ nf addLongPolys 5
               , bench "a10" $ nf addLongPolys 10
               , bench "a15" $ nf addLongPolys 15
               , bench "a20" $ nf addLongPolys 20
{-}               , bench "m1" $ nf mulLongPolys 1
               , bench "m3" $ nf mulLongPolys 3
               , bench "m5" $ nf mulLongPolys 5
               , bench "m7" $ nf mulLongPolys 7
               , bench "m9" $ nf mulLongPolys 9
               , bench "c1" $ nf convPolys 1
               , bench "c2" $ nf convPolys 2
               , bench "c3" $ nf convPolys 3
               , bench "c4" $ nf convPolys 4-}
               ]
  ]
