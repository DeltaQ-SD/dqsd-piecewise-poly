{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MonoLocalBinds #-}

{-|
Module      : IRVs
Description : Operations on improper random variables
Copyright   : (c) Peter Thompson, 2024
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

We enable extraction of centiles.

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
import PWPs.SimplePolynomials (Poly (..), makePoly, zeroPoly)
import PWPs.ConvolutionClasses 

type MyConstraints a = (Fractional a, Ord a, Num a, Enum a, Eq a)
type DistributionD a = Pieces a (PolyDelta a)
type DistributionH a = Pieces a (PolyHeaviside a)

data IRV a = PDF (DistributionD a) | CDF (DistributionH a)
    deriving (Eq, Show)

top :: (Ord a, Num a) => IRV a
top = PDF (makePieces [(0, D 1), (0, Pd (makePoly 0))])

bottom :: (Ord a, Num a) => IRV a
bottom = CDF (makePieces [(0, Ph (makePoly 0))])

cumulativeMass :: (Ord a, Enum a, Num a, Fractional a, Evaluable a (DistributionH a), Integrable (DistributionD a) (DistributionH a)) => IRV a -> a -> a
cumulativeMass x p = last $ evaluate p (makeCDF x)

invert :: (Eq a, Ord a, Fractional a) => DistributionH a -> DistributionH a
-- | Construct the inverse CDF by subtracting the CDF from 1
invert = applyObject (-) (Ph $ makePoly 1)

shiftIRV :: (Fractional a, Ord a, Num a, Enum a, Eq a, Differentiable (DistributionH a) (DistributionD a)) => a -> IRV a -> IRV a
-- | Make a delta and convolve with it
shiftIRV s x = constructDelta s PWPs.IRVs.<+> PDF (makePDF x)

makePDF :: (Fractional a, Ord a, Num a, Enum a, Eq a, Differentiable (DistributionH a) (DistributionD a)) => IRV a -> DistributionD a
-- | Force an IRV into a PDF by differentiating if necessary
makePDF (PDF x) = x
makePDF (CDF x) = differentiate x

makeCDF :: (Fractional a, Ord a, Num a, Enum a, Eq a, Integrable (DistributionD a) (DistributionH a)) => IRV a -> DistributionH a
-- | Force an IRV into a CDF by integrating if necessary
makeCDF (CDF x) = x
-- assume PDFs are 0 at 0
makeCDF (PDF x) = integrate x

constructUniform :: MyConstraints a => a -> IRV a
-- | Construct a PDF with uniform probability from 0 to the given value
constructUniform x = if x <= 0 then error "Invalid interval"
                     else PDF (makePieces [(0, Pd (makePoly (1/x))), (x, Pd zeroPoly)])

constructDelta :: MyConstraints a => a -> IRV a
-- | Construct a PDF that is a delta function at the given value
constructDelta x
    | x < 0 = error "Invalid value"
    | x == 0 = PDF (makePieces [(0, D 1), (0, Pd 0)])
    | otherwise = PDF (makePieces [(0, Pd 0), (x, D 1), (x, Pd 0)])

-- | PDF with zero value everywhere
zeroPDF :: MyConstraints a => IRV a
zeroPDF = PDF (makePieces [(0, Pd $ zeroPoly)])

monotonicFromZero :: (Ord a, Num a) => [a] -> Bool
monotonicFromZero xs = if null xs then error "Empty list" else head xs == 0 && monotonic xs

repeatedPoint :: Eq a => [a] -> Bool
repeatedPoint as = and $ zipWith (==) as (tail as)

constructCDF :: MyConstraints a => [(a, a)] -> IRV a
-- | Construct a CDF from a list of values, treating each new value as a step up from the one before, assuming we start at 0
-- | First interval is a zero polynomial: subsequent intervals start with a delta and then have a constant polynomial.
constructCDF xs
    | length xs < 2                         = error "Insufficient points"
    | not (monotonicFromZero basepoints)    = error "Basepoints not monotonic"
    | repeatedPoint basepoints              = error "Repeated basepoint"
    | not (monotonicFromZero probabilities) = error "Probabilities not monotonic"
    | last probabilities > 1                = error "Probability exceeds one"
    | otherwise = (CDF . makePieces) (interleave treads risers)
        where
            basepoints = map fst xs
            probabilities = map snd xs
            -- Each step up corresponds to a delta of the difference with the previous value
            makeStep (x, y) = H x y
            risers = zip (tail basepoints) (zipWith (curry makeStep) probabilities (tail probabilities) )
            treads = zip basepoints (map (Ph . makePoly) probabilities) -- always have one more tread than riser
            interleave :: [b] -> [b] -> [b]
            interleave [] _ = []
            interleave [x] [] = [x] -- keep the last tread
            interleave _ [] = []
            interleave (x':xs') (y':ys') = x':y':interleave xs' ys'

constructLinearCDF:: MyConstraints a => [(a, a)] -> IRV a
-- | Construct a CDF from a list of values, interpolating linearly between each pair of points
constructLinearCDF xs
    | length xs < 2                         = error "Insufficient points"
    | not (monotonicFromZero basepoints)    = error "Basepoints not monotonic"
    | repeatedPoint basepoints              = error "Repeated basepoint"
    | not (monotonicFromZero probabilities) = error "Probabilities not monotonic"
    | 0 `elem` steps                        = error "Zero-width interval"
    | last probabilities > 1                = error "Probability exceeds one"
    | otherwise = (CDF . makePieces) (zip basepoints segments)
        where
            {- 
                Each linear polynomial has the form y = sx + c, where s is given by the difference
                in successive probabilities divided by the difference in succesive basepoints.
                The constant c is fixed by the constraint that we need to pass through the point (xn,yn),
                so cn = yn - xnsn
            -}
            basepoints = map fst xs
            probabilities = map snd xs
            steps   = zipWith (-) (tail basepoints) basepoints          -- width of each interval
            stepUps = zipWith (-) (tail probabilities) probabilities    -- increments in probabilities
            slopes  = zipWith (/) stepUps steps                         -- slope of each segment
            segments  = map Ph (zipWith3 makeSegment probabilities basepoints slopes ++ [makePoly (last probabilities)])
            makeSegment y x s = Poly [y - x*s, s]

asDiscreteCDF :: MyConstraints a => IRV a -> Int -> [(a, a)]
-- | Turn an IRV into a list of point pairs corresponding to the CDF, with a given minimum number of points
asDiscreteCDF x n = if n <= 0 then error "Invalid number of points" else decomposeIRV n (makeCDF x)

asDiscretePDF :: MyConstraints a => IRV a -> Int -> [Either (a,a) [(a, a)]]
-- | Return a sequence of (Left) Impulse Probablity mass (equivalent to the
--   integral of the Heaviside function at that point) or (Right) a sequence
--   of Time and Probability Density. The sequence is monotonically increasing in Time.
asDiscretePDF x n = if n <= 0 then error "Invalid number of points"
                              else displayPolyDeltaIntervals (makePDF x) spacing
    where
        spacing = (snd (support x) - fst (support x)) / Prelude.fromIntegral n

decomposeIRV :: (Ord a, Num a, Fractional a) => Int -> DistributionH a -> [(a, a)]
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
        values = map (head . (`evaluate` ys)) basepoints

firstToFinish :: MyConstraints a => IRV a -> IRV a -> IRV a
firstToFinish x y = CDF (cdfOfx + cdfOfy - (cdfOfx * cdfOfy))
    where
        cdfOfx = makeCDF x
        cdfOfy = makeCDF y

multiFtF :: MyConstraints a => [IRV a] -> IRV a
-- | Compute the first-to-finsh of a list of IRVs by multiplying inverse CDFs then inverting
-- If there is nothing finish the result is bottom
multiFtF [] = bottom
multiFtF [x] = x
-- now know we have at least two, so head and tail are safe
multiFtF xs = CDF $ invert $ foldr (*) (head icdfs) (tail icdfs)
    where
        icdfs = map (invert . makeCDF) xs

allToFinish :: MyConstraints a => IRV a -> IRV a -> IRV a
allToFinish x y = CDF (makeCDF x * makeCDF y)

multiAtF :: MyConstraints a => [IRV a] -> IRV a
-- | Compute the last-to-finsh of a list of IRVs: if we have nothing to wait for the result is top
multiAtF [] = top
multiAtF [x] = x
-- now know we have at least two, so head and tail are safe
multiAtF xs = CDF $ foldr (*) (head cdfs) (tail cdfs)
    where
        cdfs = map makeCDF xs

probChoice :: MyConstraints a => a -> IRV a -> IRV a -> IRV a
{- | 
The probability is for choosing the left branch.
We can do this on either PDFs or CDFs; if we have CDFs deliver a CDF, if we have both PDFs or one of each, deliver a PDF.
-}
probChoice p x y =
    case (x,y) of
        (CDF a, CDF b) -> CDF ((p >< a) + ((1 - p) >< b))
        _              -> PDF ((p >< makePDF x) + ((1 - p) >< makePDF y))

multiWeightedChoice :: MyConstraints a => [(a, IRV a)] -> IRV a
-- If we have nothing to choose from we get nothing
multiWeightedChoice [] = bottom
-- we'll force everything into PDFs and deliver a PDF
multiWeightedChoice xs = PDF (sum (zipWith adjust weights pdfs))
    where
        weights = map fst xs                -- :: [a]
        pdfs    = map (makePDF . snd) xs    -- :: [DistributionD a]
        adjust x y = (x / sum weights) >< y

-- | To convolve, force into PDFs and then invoke piecewise convolution
infix 7 <+> -- same as *
(<+>) :: MyConstraints a => IRV a -> IRV a -> IRV a
x <+> y = PDF (makePDF x PWPs.Piecewise.<+> makePDF y)

probMass :: MyConstraints a => IRV a -> a
-- | Just extract the final value of the CDF
probMass = piecesFinalValue . makeCDF

compareIRVs :: MyConstraints a => IRV a -> IRV a -> Maybe Ordering
{- | 
    If the two IRVs are partially ordered, return an ordering, otherwise return Nothing.
-}
compareIRVs x y = comparePW (makeCDF x) (makeCDF y)

support :: (Eq a, Fractional a) => IRV a -> (a, a)
-- | return the first and last basepoints for which the value is significant
support (PDF x) = piecewiseSupport x
support (CDF x) = piecewiseSupport x

centiles :: (Fractional a, Ord a, Num a, Enum a, Eq a, Evaluable a (PolyHeaviside a)) => [a] -> IRV a -> [Maybe a]
-- | Given a list of probabiity values, return the times at which each value is reached; 
-- if it is never reached, return Nothing in that position
centiles probabilities dQ
    | null probabilities            = error "Empty probability list"
    | not (monotonic probabilities) = error "Probabilities not monotonic"
    | last probabilities > 1        = error "Probability exceeds one"
    -- deal with degenerate case where the support has zero width: all centiles the same value
    | otherwise = if eps == 0 then replicate (length probabilities) (Just $ fst (support dQ))
                              else findCentiles probabilities
        where
            intervals = getPieces $ makeCDF dQ
            eps = snd (support dQ) - fst (support dQ) / 1000 -- arbitrary parameter!
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
            findRoot p (first:rest@(second:_)) = 
                -- if the value is in the range of the current interval,
                if firstValue <= p && p <= secondValue
                -- find the root,
                then polyHeavisideRoot eps p (basepoint first, basepoint second) (object first)
                -- otherwise move on to the next piece
                else findRoot p rest
                where
                    firstValue  = head $ evaluate (basepoint first) (object first)
                    secondValue = last $ evaluate (basepoint second) (object first)
