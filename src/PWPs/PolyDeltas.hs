{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

{-|
Module      : Deltas
Description : Polynomials extended with delta functions
Copyright   : (c) Peter Thompson, 2024
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
    , polyDeltaRoot
)
where
import PWPs.ConvolutionClasses
import PWPs.SimplePolynomials as SP

{- |
A PolyDelta either a polynomial, a (shifted, scaled) Delta or a (shifted, scaled) Heaviside. 
A delta has a mass, and a Heaviside has a starting value and a rise; 
for probabilities all should be constrained between 0 and 1. 
The position of both Ds and Hs is stored as its basepoint when doing piecewise operations.
-} 
data PolyDelta a = P (Poly a) | D a | H a a 
    deriving (Eq, Show)

instance Functor PolyDelta where
    fmap f (D x) = D (f x)
    fmap f (H x y) = H (f x) (f y)
    fmap f (P x) = P (fmap f x)

plusPD :: (Eq a, Fractional a) => PolyDelta a -> PolyDelta a -> PolyDelta a
-- Polynomials have zero mass at a single point, so they are dominated by Ds and Hs
plusPD (P x) (P y) = P (SP.plus x y)
plusPD (P _) (H x y) = H x y
plusPD (P _) (D x) = D x
plusPD (D x) (P _) = D x
plusPD (D x) (D x') = D (x + x') 
plusPD (H x y) (P _) = H x y
plusPD (H x y) (H x' y') = H (x + x') (y + y')
plusPD _ _ = error "Cannot add Delta to Heaviside"

timesPD :: (Eq a, Fractional a) => PolyDelta a -> PolyDelta a -> PolyDelta a
timesPD (P x) (P y) = P (SP.times x y)
timesPD (P _) (D x) = D x
timesPD (D x) (P _) = D x
timesPD (D x) (D x') = D (x * x')
timesPD (P _) (H x y) = H x y
timesPD (H x y) (P _) = H x y
timesPD (H x y) (H x' y') = H (x * x') (y * y')
timesPD _ _ = error "Cannot multiply Delta by Heaviside"

minusPD :: Num a => PolyDelta a -> PolyDelta a
minusPD = fmap negate

scalePD :: Num a => a -> PolyDelta a -> PolyDelta a
scalePD x (P a) = P (SP.scalePoly x a)
scalePD x (D y) = D (x * y)
scalePD x (H y z) = H (x * y) (x * z)

evaluatePD :: Num a => a -> PolyDelta a -> [a]
evaluatePD point (P x) = [SP.evaluatePoly point x]
evaluatePD _ (H x y) = [x, y] 
evaluatePD _ (D x) = [x]

integratePD :: (Eq a, Fractional a) => PolyDelta a -> PolyDelta a
integratePD (P x) = P (SP.integrate x)
integratePD (D x) = H 0 x 
integratePD (H _ _) = error "Integration of a Heaviside disallowed" -- would require more sophisticated joining of pieces

differentiatePD :: (Eq a, Num a, Fractional a) => PolyDelta a -> PolyDelta a
differentiatePD (P x) = P (SP.differentiate x)
differentiatePD (H x y) = D (y - x)
differentiatePD (D _) = error "Differentiation of Delta is illegal"

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
        zero             = P 0 -- is this the best choice?
        fromInteger n    = D (Prelude.fromInteger n) -- is this the best choice?
        differentiate    = differentiatePD
        integrate        = integratePD 

boostPD :: (Eq a, Num a, Fractional a) => a -> PolyDelta a -> PolyDelta a
boostPD x (P y) = plusPD (P y) (P (makePoly x))
boostPD _ (D y) = D y
boostPD x (H y z) = H (x + y) (x + z)

instance (Eq a, Num a, Fractional a) => Evaluable a (PolyDelta a) 
    where
        evaluate  = evaluatePD
        boost     = boostPD
        scale     = scalePD

-- | Removes excess basepoints if the objects on either side are the same
aggregate :: Eq a => [(a, PolyDelta a)] -> [(a, PolyDelta a)]
aggregate []    = error "Empty list of polydeltas"
aggregate [x]   = [x] -- need at least two elements to do anything
aggregate ((bx, x):(by, y):xs)
    | x == y    = aggregate $ (bx, x):xs -- throw away the second basepoint
    | otherwise = (bx, x) : aggregate ((by, y):xs)

convolvePolyDeltas :: (Num a, Fractional a, Ord a)
                   => (a, a, PolyDelta a) -> (a, a, PolyDelta a) -> [(a, PolyDelta a)]
{- |
When both arguments are polynomials, we check the intervals are non-zero then use convolvePolys and just map the type.
For a delta, lower == upper (invariant to be checked), and the effect of the delta is to translate the other
argument (whichever it is) along by this amount. Need to ensure there is still an initial interval based at zero.
-} 
convolvePolyDeltas (lf, uf, P f) (lg, ug, P g) = 
    if (uf <= lf) || (ug <= lg) then error "Invalid polynomial interval width"
                                else aggregate $ map (\(x, p) -> (x, P p)) (convolvePolys (lf, uf, f) (lg, ug, g))
convolvePolyDeltas (lf, uf, D f) (lg, ug, P g) 
    | lf /= uf     = error "Non-zero delta interval"
    | ug < lg      = error "Negative interval width"
    | f == 0       = [(0, P zero)] -- convolving with a zero-sized delta gives nothing
    | lg + lf == 0 = [(0, scalePD f (P g)), (ug, P zero)] -- degenerate case of delta at zero
    -- Shift the poly by the basepoint of the delta and insert a new initial zero interval
    | otherwise    = aggregate [(0, P zero), (lg + lf, scalePD f (P (shiftPoly lf g))), (ug + lf, P zero)]
convolvePolyDeltas (lf, uf, P f) (lg, ug, D g) = convolvePolyDeltas (lg, ug, D g) (lf, uf, P f)  -- commutative
convolvePolyDeltas (lf, uf, D f) (lg, ug, D g) -- both deltas
    | lf /= uf || lg /= ug  = error "Non-zero delta interval"
    | f * g == 0            = [(0, P zero)] -- convolving with a zero-sized delta gives nothing
    | lg + lf == 0          = [(0, D (f * g)), (0, P zero)] -- degenerate case of deltas at zero
-- convolving with Heavisides is forbidden
convolvePolyDeltas _ _ = error "Unexpected convolution case"

instance (Num a, Fractional a, Ord a) => CompactConvolvable a (PolyDelta a)
    where
        convolveIntervals = convolvePolyDeltas

{-|
    We measure whether or not a polydelta is consistently above or below zero, or equals zero
-}
comparePDToZero :: (Fractional a, Eq a, Ord a) => (a, a, PolyDelta a) -> Maybe Ordering
comparePDToZero (lf, uf, P f) = SP.compareToZero (lf, uf, f) -- simple polynomial case
comparePDToZero (lf, uf, D f) 
    | lf /= uf      = error "Non-zero Delta interval"
    | f == 0        = Just EQ
    | f > 0         = Just GT
    | otherwise     = Just LT
comparePDToZero (lf, uf, H x y)
    | lf /= uf      = error "Non-zero Heaviside interval"
    | (x + y) == 0  = Just EQ
    | (x + y) > 0   = Just GT
    | otherwise     = Just LT

instance (Fractional a, Eq a, Ord a) => Comparable a (PolyDelta a)
    where
        compareZero = comparePDToZero

{-|
    We merge polynomials if they are equal. We merge Ds/Hs by adding them (not that we expect this case).
    We merge a zero D/H with a polynomial by discarding it. Other cases do not merge.
-}
instance (Num a, Eq a, Fractional a) => Mergeable (PolyDelta a)
    where
        mergeObject a b = case (a, b) of
            (P x, P y) -> if x == y then Just (P y) else Nothing
            (D x, P y) -> if x == 0 then Just (P y) else Nothing
            (H _ x, P y) -> if x == 0 then Just (P y) else Nothing
            (D x, D y) -> Just (D (x + y))
            (H x y, H x' y') -> Just (H (x + x') (y + y'))
            (_, _) -> Nothing
        zeroObject = zero

{-|
    Given an interval containing a given value of a PolyDelta, find its location
-}
polyDeltaRoot :: (Ord a, Num a, Eq a, Fractional a) => a -> a -> (a, a) -> PolyDelta a -> Maybe a
-- If we have a step, the interval is zero width so this is the root
polyDeltaRoot _ _ (l, u) (H _ _) = if l /= u then error "Non-zero Heaviside interval" else Just l 
-- otherwise we have a polynomial: subtract the value we are looking for so that we seek a zero crossing
polyDeltaRoot e x (l, u) (P p) = findPolyRoot e (l, u) (p `plus` makePoly (-x))
polyDeltaRoot _ _ _ (D _) = error "Can't take the root of a delta"

displayPolyDelta :: (Ord a, Num a, Eq a, Fractional a) => a -> (a, a, PolyDelta a) -> Either (a,a) [(a, a)]
displayPolyDelta _ (l, u, D x)   = if l /= u then error "Non-zero delta interval"
                                    else Left (l, x)
displayPolyDelta s (l, u, P p)   = if l >= u then error "Invalid polynomial interval"
                                    else Right (displayPoly p (l, u) s)
displayPolyDelta _ (l, u, H x y) = if l /= u then error "Non-zero heaviside interval"
                                    else Left (l, y - x)                                

instance (Ord a, Num a, Eq a, Fractional a) => Displayable a (PolyDelta a)
    where
        displayObject = displayPolyDelta