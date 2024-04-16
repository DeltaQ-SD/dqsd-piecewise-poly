{-# LANGUAGE TypeSynonymInstances #-}
{-|
Module      : IRVs
Description : Operations on improper random variables
Copyright   : (c) Peter Thompson, 2023
License     : BSD-2-Clause
Maintainer  : peter.thompson@pnsol.com
Stability   : experimental

We build IRVs as ab abstract datatype on top of piecewise polynomials and deltas. 
IRVs may be represented as PDFs or CDFs. They are polymorphic in a constrained numeric type.

We provide functions to construct IRVs as uniform distributions on an interval; a delta at a point;
and a CDF from a list of point-pairs.

We provide functions to turn an IRV into a set of point-pairs representing either the PDF or CDF.

We define operators for convolution, first-to-finish, all-to-finish and weighted choice.
We also provide an operation to extract the probability mass (<= 1).

We provide an implementation of the partial order on IRVs.

Note that we can consistently represent 'bottom' using polynomials and 'top' by including deltas.
-}
module PWPs.IRVs
(
    IRV
  , makePDF
  , makeCDF
  , constructUniform
  , constructDelta
  , zeroPDF
  , constructCDF
  , constructLinearCDF
  , firstToFinish
  , multiFtF
  , allToFinish
  , multiAtF
  , probChoice
  , multiWeightedChoice
  , (PWPs.IRVs.<+>)
  , probMass
  , invertCDF
  , asDiscreteCDF
  , asDiscretePDF
  , compareIRVs
  , support
  , top
  , bottom
  , centiles
  , cumulativeMass
  , shiftIRV
)
where

import PWPs.Piecewise
import PWPs.PolyDeltas
import PWPs.SimplePolynomials (Poly (..), makePoly)

type Distribution a = Pieces a (PolyDelta a)

data IRV a = PDF (Distribution a) | CDF (Distribution a)
    deriving (Eq, Show)

top :: (Ord a, Num a) => IRV a
top = PDF (makePieces [(0, D 1), (0, P (makePoly 0))])

bottom :: (Ord a, Num a) => IRV a
bottom = CDF (makePieces [(0, P (makePoly 0))])

cumulativeMass :: (Ord a, Enum a, Num a, Fractional a) => IRV a -> a -> a
cumulativeMass x p = last $ evaluateAtApoint p (makeCDF x)

invert :: (Eq a, Ord a, Fractional a) => Distribution a -> Distribution a
-- | Construct the inverse CDF by subtracting the CDF from 1
invert = applyObject plus (P $ makePoly 1) . minus

invertCDF :: (Eq a, Ord a, Fractional a, Enum a) => IRV a -> IRV a
invertCDF x = CDF (invert $ makeCDF x)

shiftIRV :: (Ord a, Enum a, Num a, Fractional a, Num a) => a -> IRV a -> IRV a
-- | Make a delta and convolve with it
shiftIRV s x = constructDelta s PWPs.IRVs.<+> PDF (makePDF x)

makePDF :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a -> Distribution a
-- | Force an IRV into a PDF by differentiating if necessary
makePDF (PDF x) = x
makePDF (CDF x) = differentiate x

makeCDF :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a -> Distribution a
-- | Force an IRV into a CDF by integrating if necessary
makeCDF (CDF x) = x
-- assume PDFs are 0 at 0
makeCDF (PDF x) = integrate x

constructUniform :: (Ord a, Enum a, Eq a, Fractional a, Num a) => a -> IRV a
-- | Construct a PDF with uniform probability from 0 to the given value
constructUniform x = if x <= 0 then error "Invalid interval"
                     else PDF (makePieces [(0, P (makePoly (1/x))), (x, P zero)])

constructDelta :: (Ord a, Enum a, Eq a, Fractional a, Num a) => a -> IRV a
-- | Construct a PDF that is a delta function at the given value
constructDelta x
    | x < 0 = error "Invalid value"
    | x == 0 = PDF (makePieces [(0, D 1), (0, P 0)])
    | otherwise = PDF (makePieces [(0, P 0), (x, D 1), (x, P 0)])

-- | PDF with zero value everywhere
zeroPDF :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a
zeroPDF = PDF (makePieces [(0, P $ makePoly 0)])

monotonicFromZero :: (Ord a, Num a) => [a] -> Bool
monotonicFromZero xs = if null xs then error "Empty list" else head xs == 0 && monotonic xs

constructCDF :: (Ord a, Enum a, Eq a, Fractional a, Num a) => [(a, a)] -> IRV a
-- | Construct a CDF from a list of values, treating each new value as a step up from the one before, assuming we start at 0
-- | First interval is a zero polynomial: subsequent intervals start with a delta and then have a constant polynomial.
constructCDF xs
    | length xs < 2                         = error "Insufficient points"
    | not (monotonicFromZero basepoints)    = error "Basepoints not monotonic"
    | not (monotonicFromZero probabilities) = error "Probabilities not monotonic"
    | last probabilities > 1                = error "Probability exceeds one"
    | otherwise = (CDF . makePieces) (interleave treads risers)
        where
            basepoints = map fst xs
            probabilities = map snd xs
            -- Each step up corresponds to a delta of the difference with the previous value
            risers = zip (tail basepoints) (map D (zipWith (-) (tail probabilities) probabilities))
            treads = zip basepoints (map (P . makePoly) probabilities) -- always have one more tread than riser
            interleave :: [b] -> [b] -> [b]
            interleave [] _ = []
            interleave [x] [] = [x] -- keep the last tread
            interleave _ [] = []
            interleave (x':xs') (y':ys') = x':y':interleave xs' ys'

constructLinearCDF:: (Ord a, Enum a, Eq a, Fractional a, Num a) => [(a, a)] -> IRV a
-- | Construct a CDF from a list of values, interpolating linearly between each pair of points
constructLinearCDF xs
    | length xs < 2                         = error "Insufficient points"
    | not (monotonicFromZero basepoints)    = error "Basepoints not monotonic"
    | not (monotonicFromZero probabilities) = error "Probabilities not monotonic"
    | 0 `elem` steps                        = error "Zero-width interval"
    | last probabilities > 1                = error "Probability exceeds one"
    | otherwise = (CDF . makePieces) (zip basepoints slopes)
        where
            basepoints = map fst xs
            probabilities = map snd xs
            steps   = zipWith (-) (tail basepoints) basepoints -- width of each interval
            stepUps = zipWith (-) (tail probabilities) probabilities -- increments in probabilities
            slopes  = map P (zipWith3 slope probabilities stepUps steps ++ [makePoly (last probabilities)])
            -- each polynomial starts at the current probability and linearly slopes up to the next one
            slope x y z = Poly [x, y/z]

asDiscreteCDF :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a -> Int -> [(a, a)]
-- | Turn an IRV into a list of point pairs corresponding to the CDF, with a given minimum number of points
asDiscreteCDF x n = if n <= 0 then error "Invalid number of points" else decomposeIRV n (makeCDF x)

asDiscretePDF :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a -> Int -> [Either (a,a) [(a, a)]]
-- | Return a sequence of (Left) Impulse Probablity mass (equivalent to the
--   integral of the Heaviside function at that point) or (Right) a sequence
--   of Time and Probability Density. The sequence is monotonically increasing in Time.
asDiscretePDF x n = if n <= 0 then error "Invalid number of points"
                              else displayPolyDeltaIntervals (makePDF x) spacing
    where
        spacing = (snd (support x) - fst (support x)) / Prelude.fromIntegral n

decomposeIRV :: (Ord a, Num a, Fractional a) => Int -> Distribution a -> [(a, a)]
decomposeIRV numPoints ys = zip basepoints values
    where
        pointsList = getPieces ys
        originalBasepoints = map basepoint pointsList
        spacing = last originalBasepoints / Prelude.fromIntegral numPoints
        basepoints = makePoints originalBasepoints spacing (head originalBasepoints)
        makePoints :: (Ord b, Num b) => [b] -> b -> b -> [b]
        -- make points by repeatedly adding spacing to last interval boundary until it exceeds the next interval boundary
        makePoints [] _ _ = [] -- already got the last interval boundary included, no need to go further
        makePoints (y':ys') sp prev =
            if new >= y'
                then y':makePoints ys' sp y' -- gone past the next interval, so take that instead and carry on
                else new:makePoints (y':ys') sp new -- just add spacing and repeat until we pass the boundary
            where
                new = prev + sp
        -- evaluate at the the basepoints - will get a list at each point, take the first
        values = map (head . (`evaluateAtApoint` ys)) basepoints

firstToFinish :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a -> IRV a -> IRV a
firstToFinish x y = CDF (cdfOfx `plus` cdfOfy `plus` minus (cdfOfx `times` cdfOfy))
    where
        cdfOfx = makeCDF x
        cdfOfy = makeCDF y

multiFtF :: (Ord a, Enum a, Eq a, Fractional a, Num a) => [IRV a] -> IRV a
-- | Compute the first-to-finsh of a list of IRVs by multiplying inverse CDFs then inverting
-- If there is nothing finish the result is bottom
multiFtF [] = bottom
multiFtF [x] = x
-- now know we have at least two, so head and tail are safe
multiFtF xs = CDF $ invert $ foldr times (head icdfs) (tail icdfs)
    where
        icdfs = map (invert . makeCDF) xs

allToFinish :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a -> IRV a -> IRV a
allToFinish x y = CDF (makeCDF x `times` makeCDF y)

multiAtF :: (Ord a, Enum a, Eq a, Fractional a, Num a) => [IRV a] -> IRV a
-- | Compute the last-to-finsh of a list of IRVs: if we have nothing to wait for the result is top
multiAtF [] = top
multiAtF [x] = x
-- now know we have at least two, so head and tail are safe
multiAtF xs = CDF $ foldr times (head cdfs) (tail cdfs)
    where
        cdfs = map makeCDF xs

probChoice :: (Ord a, Enum a, Eq a, Fractional a, Num a) => a -> IRV a -> IRV a -> IRV a
{- | 
The probability is for choosing the left branch.
We can do this on either PDFs or CDFs; if we have CDFs deliver a CDF, if we have both PDFs or one of each, deliver a PDF.
-}
probChoice p x y =
    case (x,y) of
        (CDF a, CDF b) -> CDF ((p >< a) `plus` ((1 - p) >< b))
        _              -> PDF ((p >< makePDF x) `plus` ((1 - p) >< makePDF y))

multiWeightedChoice :: (Ord a, Enum a, Eq a, Fractional a, Num a) => [(a, IRV a)] -> IRV a
-- If we have nothing to choose from we get nothing
multiWeightedChoice [] = bottom
-- we'll force everything into PDFs and deliver a PDF
multiWeightedChoice xs = PDF (foldr plus zero (zipWith adjust weights pdfs))
    where
        weights = map fst xs                -- :: [a]
        pdfs    = map (makePDF . snd) xs    -- :: [Distribution a]
        adjust x y = (x / sum weights) >< y

-- | To convolve, force into PDFs and then invoke piecewise convolution
infix 7 <+> -- same as *
(<+>) :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a -> IRV a -> IRV a
x <+> y = PDF (makePDF x PWPs.Piecewise.<+> makePDF y)

probMass :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a -> a
-- | Just extract the final value of the CDF
probMass = piecesFinalValue . makeCDF

compareIRVs :: (Ord a, Enum a, Eq a, Fractional a, Num a) => IRV a -> IRV a -> Maybe Ordering
{- | 
    If the two IRVs are partially ordered, return an ordering, otherwise return Nothing.
    Ordering is preserved through integration and differentiation, so go use either PDF
    or CDF - CDF -> PDF is cheaper so use PDFs.
-}
compareIRVs x y = comparePW (makePDF x) (makePDF y)

support :: (Eq a, Fractional a) => IRV a -> (a, a)
-- | return the first and last basepoints for which the value is significant
support (PDF x) = piecewiseSupport x
support (CDF x) = piecewiseSupport x

centiles :: (Ord a, Enum a, Eq a, Fractional a, Num a) => [a] -> IRV a -> [Maybe a]
-- | Given a list of probabiity values, return the times at which that value is reached; 
-- if it is never reached, return Nothing
centiles probabilities dQ
    | null probabilities            = error "Empty probability list"
    | not (monotonic probabilities) = error "Probabilities not monotonic"
    | last probabilities > 1        = error "Probability exceeds one"
    -- deal with degenerate case where the support has zero width: all centiles the same value
    | otherwise = if eps == 0 then replicate (length probabilities) (Just $ fst (support dQ))
                              else reverse $ findCentiles probabilities
        where
            intervals = getPieces $ makeCDF dQ
            eps = fst (support dQ) - snd (support dQ) / 1000 -- arbitrary parameter!
            -- findCentiles :: [a] -> [Maybe a]
            -- consume the list of probabilities and build the list of centiles
            findCentiles [] = []
            findCentiles (x:xs) = if x > probMass dQ then Nothing : findCentiles xs
                                   else findRoot x intervals : findCentiles xs
            -- Work along cumulative distribution until we find the interval containing the value,  
            -- then get the root of the polydelta
            -- findRoot :: (Ord a, Enum a, Eq a, Fractional a, Num a) => a -> [Piece a (PolyDelta a)] -> a
            findRoot _ [] = error "Empty distribution"
            -- if we have only one piece its object must be constant, so report its basepoint as the root
            findRoot _ [final] = Just (basepoint final)
            -- hereon we must have at least two pieces: if the value is in the range of the current interval,
            -- find the root, otherwise move on to the next piece
            findRoot p (next:rest) = if inInterval p next (head rest)
                                     then polyDeltaRoot eps p (basepoint next, basepoint $ head rest) (object next)
                                     else findRoot p rest
            -- inInterval :: (Ord a, Enum a, Eq a, Fractional a, Num a) => a -> Piece a (PolyDelta a) -> Piece a (PolyDelta a) -> Bool
            inInterval y first second = firstValue <= y && y <= secondValue
                where
                    firstValue  = head $ evaluatePD (basepoint first) (object first)
                    secondValue = last $ evaluatePD (basepoint second) (object first)