{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-|
Module      : Deltas
Description : Polynomials extended with delta functions
Copyright   : (c) Peter Thompson, 2023
License     : BSD-2-Clause
Maintainer  : peter.thompson@pnsol.com
Stability   : experimental

We extend the polynomials with delta functions so that we can have 'pieces' of either.
We can add and multiply polynomials and add deltas; multiplying deltas only makes sense
in the context of multiplying CDFs, when it will be OK.
We will only be combining a delta and a poly over a zero interval, where the delta dominates.
We can define the operators in a way that keeps the deltas intact and makes the
piecewise integration/differentiation work.
-}
module PWPs.PolyDeltas 
(
      PolyDelta (..)
    , plus
    , times
    , minus
 --   , fromInteger
    , scalePD
    , differentiate
    , integrate
    , evaluatePD
    , convolvePolyDeltas
    , comparePDToZero
)
where
import PWPs.ConvolutionClasses
import PWPs.SimplePolynomials as SP

{- |
Represents either a polynomial or a shifted dirac delta. A dirac delta contains the "size"; this is
for probabilistic choice (and as such, it should be constrained between 0 and 1). Its position is
stored as its basepoint when doing piecewise operations.
-} 
data PolyDelta a = P (Poly a) | D a
    deriving (Eq, Show)

instance Functor PolyDelta where
    fmap f (D x) = D (f x)
    fmap f (P x) = P (fmap f x)

plusPD :: (Eq a, Fractional a) => PolyDelta a -> PolyDelta a -> PolyDelta a
plusPD (P x) (P y) = P (SP.plus x y)
plusPD (P _) (D y) = D y
plusPD (D x) (P _) = D x
plusPD (D x) (D y) = D (x + y)

timesPD :: (Eq a, Fractional a) => PolyDelta a -> PolyDelta a -> PolyDelta a
timesPD (P x) (P y) = P (SP.times x y)
timesPD (P _) (D y) = D y
timesPD (D x) (P _) = D x
timesPD (D x) (D y) = D (x * y)

minusPD :: Num a => PolyDelta a -> PolyDelta a
minusPD = fmap negate

scalePD :: Num a => a -> PolyDelta a -> PolyDelta a
scalePD x (P a) = P (SP.scalePoly x a)
scalePD x (D y) = D (x * y)

evaluatePD :: Num a => a -> PolyDelta a -> a
evaluatePD point (P x) = SP.evaluatePoly point x
evaluatePD _ (D x) = x -- this will make the piecewise integration work

integratePD :: (Eq a, Fractional a) => PolyDelta a -> PolyDelta a
integratePD (P x) = P (SP.integrate x)
integratePD (D x) = D x

differentiatePD :: (Eq a, Num a, Fractional a) => PolyDelta a -> PolyDelta a
differentiatePD (P x) = P (SP.differentiate x)
differentiatePD (D x) = D x

instance (Eq a, Num a, Fractional a) => Num (PolyDelta a) where
    (+)           = plusPD
    (*)           = timesPD
    negate        = minusPD
    abs           = undefined
    signum        = undefined
    fromInteger n = D $ Prelude.fromInteger n

instance (Eq a, Num a, Fractional a) => Calculable (PolyDelta a)
    where
        plus             = plusPD
        times            = timesPD
        minus            = minusPD
        zero             = D 0 -- is this the best choice?
        fromInteger n    = D (Prelude.fromInteger n) -- is this the best choice?
        differentiate    = differentiatePD
        integrate        = integratePD 

instance (Eq a, Num a, Fractional a) => Evaluable a (PolyDelta a) 
    where
        evaluate  = evaluatePD
        boost x y = plusPD y (P (makePoly x))
        scale     = scalePD

-- | Removes excess basepoints if the objects on either side are the same
aggregate :: Eq a => [(a, PolyDelta a)] -> [(a, PolyDelta a)]
aggregate ((bx, x):(by, y):xs)
    | x == y    = aggregate $ (bx, x):xs
    | otherwise = (bx, x) : aggregate ((by, y):xs)
aggregate xs = xs

convolvePolyDeltas :: (Num a, Fractional a, Ord a)
                   => (a, a, PolyDelta a) -> (a, a, PolyDelta a) -> [(a, PolyDelta a)]
{- |
When both arguments are polynomials, we use convolvePolys and just map the type.
For a delta, lower == upper (invariant to be checked), and the effect of the delta is to translate the other
argument (whichever it is) along by this amount. Need to ensure there is still an initial interval based at zero.
-} 
convolvePolyDeltas (lf, uf, P f) (lg, ug, P g) = 
    if (uf < lf) || (ug < lg) then error "Negative interval width"
                              else aggregate $ map (\(x, p) -> (x, P p)) (convolvePolys (lf, uf, f) (lg, ug, g))
convolvePolyDeltas (lf, uf, D f) (lg, ug, g) 
    | lf /= uf     = error "Non-zero delta interval"
    | ug < lg      = error "Negative interval width"
    | f == 0       = [(0, P zero)] -- convolving with a zero-sized delta gives nothing
    | lg + lf == 0 = [(0, scalePD f g), (ug, P zero)] -- degenerate case
    | otherwise    = aggregate [(0, P zero), (lg + lf, scalePD f g), (ug + lf, P zero)]
convolvePolyDeltas (lf, uf, f) (lg, ug, D g) = convolvePolyDeltas (lg, ug, D g) (lf, uf, f)  -- commutative

instance (Num a, Fractional a, Ord a) => CompactConvolvable a (PolyDelta a)
    where
        convolveIntervals = convolvePolyDeltas

{-|
    We measure whether or not a polydelta is consistently above or below zero, or equals zero
-}
comparePDToZero :: (Fractional a, Eq a, Ord a) => (a, a, PolyDelta a) -> Maybe Ordering
comparePDToZero (lf, uf, P f) = SP.compareToZero (lf, uf, f) -- simple polynomial case
comparePDToZero (lf, uf, D f) 
    | lf /= uf      = error "Non-zero delta interval"
    | f == 0        = Just EQ
    | f > 0         = Just GT
    | otherwise     = Just LT

instance (Fractional a, Eq a, Ord a) => Comparable a (PolyDelta a)
    where
        compareZero = comparePDToZero