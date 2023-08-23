{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}

{-|
Module      : Piecewise
Description : Sequences of objects defined over contiguous intervals
Copyright   : (c) Peter Thompson, 2023
License     : BSD-2-Clause
Maintainer  : peter.thompson@pnsol.com
Stability   : experimental

We consider 'objects' (in this case polynomials or deltas)
defined over numerical intervals, where we know how to combine objects
when the intervals are identical. When the intervals are not identical we construct the
overlaps so as to produce a combined object over a more complex shared set of intervals.

We assume the sets of intervals have a common initial point, otherwise we could not combine
them over the left-hanging piece. Note this means a piecewise list is never empty.

Each set of intervals is assumed to be in increasing order.

We allow for zero-width intervals in order to accomodate delta functions. A zero-width
interval may *only* contain a delta, and *must* be followed by a non-zero interval.
A non-zero-width interval may *not* contain a delta.

When combining piecewise objects we allow that one set of intervals may finish before the other, 
and assume that whatever the operation is will simply reproduce the tail of the longer sequence.
-}
module PWPs.Piecewise
(
      Pieces (..)
    , Piece (..)
    , makePieces
    , combinePieces
    , plus 
    , times
    , minus
    , zero
    , (><)
    , differentiate
    , integrate
    , (<+>)
    , piecesFinalValue
    , evaluateAtApoint
    , convolvePieces
    , disaggregate
    , monotonic
) where

import PWPs.SimplePolynomials hiding (scale, plus, times, minus, integrate, differentiate, evaluatePoly)
import PWPs.ConvolutionClasses 
import PWPs.PolyDeltas

data Piece a o = Piece
    {
        basepoint :: a
      , object :: o
    }
    deriving (Eq, Show)

instance Functor (Piece a) where
    --fmap :: (o -> b) -> Piece a o -> Piece a b
    fmap f x = x {object = f (object x)}

makePiece :: (a, o) -> Piece a o -- useful to have a tuple that we can have a list of to map over
makePiece (x, y) = Piece {basepoint = x, object = y}

newtype Pieces a o = Pieces { getPieces :: [Piece a o] }
    deriving (Eq, Show)
instance Functor (Pieces a) where
    fmap f = Pieces . map (fmap f) . getPieces

instance (Num a, Eq a, Ord a) => Applicative (Pieces a) where
    --pure :: Num a => a1 -> Pieces a a1
    pure f = Pieces [makePiece (0, f)]
    --(<*>) :: (Num a, Eq a, Ord a) => Pieces a (a1 -> b) -> Pieces a a1 -> Pieces a b
    f <*> g = combinePieces f g 

combinePieces :: (Num a, Eq a, Ord a) => Pieces a (a1 -> b) -> Pieces a a1 -> Pieces a b
{- |
We combine two piecewise objects by splitting intervals to obtain a consistent set.
Where one list has just a single element we just process the tail of the other list against the first list's terminal object.
-}
combinePieces f g = Pieces (doCombine (getPieces f) (getPieces g))
    where
        doCombine :: (Num a, Eq a, Ord a) => [Piece a (a1 -> b)] -> [Piece a a1] -> [Piece a b]
        doCombine [] _ = error "Empty piece list" -- lists should never be empty
        doCombine _ [] = error "Empty piece list"
        doCombine xs [y] = if basepoint (head xs) == basepoint y then map (\x_pc -> fmap (\_ -> object x_pc $ object y) x_pc) xs
                            else error "Initial points not coincident"
        doCombine [x] ys
            | by0 == bx = map ((fmap . object) x) ys
            | by0 <  bx = error "Initial points not coincident"
            | by0 >  bx = map ((fmap . object) x) (makePiece (bx, object $ head ys):ys)
            where bx = basepoint x
                  by0 = basepoint $ head ys
            
            
            --if basepoint (head ys) == basepoint x then map ((fmap . object) x) ys 
              --              else error "Initial points not coincident"
        -- now we know that both lists have at least two elements
        doCombine x@(x0:x1s@(x1:_)) y@(y0:y1s@(y1:_))
            | bx0 /= by0 = error "Initial points not coincident" -- invariant check
            -- If the second points are not identical, split the longer piece so that they now are and try again
            | bx1 > by1 = doCombine (x0:makePiece (by1, ox0):x1s) y
            | bx1 < by1 = doCombine x (y0:makePiece (bx1, oy0):y1s)
            -- bx1 == by1 - when the initial intervals are identical, we combine their objects over that interval and move on
            | otherwise = makePiece (bx0, ox0 oy0):doCombine x1s y1s
                where
                    bx0 = basepoint x0  -- invariant: always == basepoint y0
                    by0 = basepoint y0
                    bx1 = basepoint x1  -- might be = bx0 if we have a delta
                    by1 = basepoint y1  -- might be = by0 if we have a delta
                    ox0 = object x0
                    oy0 = object y0

monotonic :: Ord a => [a] -> Bool 
-- | Check that a list of values is monotonic (not strict to allow deltas)
monotonic ys = and (zipWith (<=) ys (tail ys))

makePieces :: (Num a, Eq a, Ord a) => [(a, o)] -> Pieces a o
-- check the basepoints are in order and start at 0
makePieces xs 
    | null xs                            = error "Provided list was empty"
    | fst (head xs) /= 0                 = error "List did not start at zero"
    | length xs == 1                     = Pieces [makePiece (head xs)] -- deal with the single element list case
    | not (monotonic (map fst xs))       = error "Basepoints were not in order" 
    | otherwise                          = Pieces (map makePiece xs)

-- | Use the applicative instance to construct the calculable instance.
instance (Num a, Enum a, Ord a, Fractional a) => Calculable (Pieces a (PolyDelta a))
    where
        plus x y        = plus <$> x <*> y
        times x y       = times <$> x <*> y
        minus x         = minus <$> x
        zero            = Pieces [Piece {basepoint = 0 :: a, object = P PWPs.SimplePolynomials.zero :: PolyDelta a}]
        fromInteger n   = Pieces [Piece {basepoint = 0 :: a, object = PWPs.ConvolutionClasses.fromInteger n}]
        differentiate   = differentiatePieces
        integrate       = integratePieces
 
instance (Num a, Enum a, Ord a, Fractional a) => Convolvable (Pieces a (PolyDelta a))
    where
        (<+>)           = convolvePieces

-- | Piecewise differentiation is easy: just differentiate all the objects
differentiatePieces :: (Enum a, Eq a, Fractional a, Num a) => Pieces a (PolyDelta a) -> Pieces a (PolyDelta a)
differentiatePieces = fmap differentiate

{- |
Piecewise integration is harder - need to evaluate at the boundary points to make the pieces join up.
We need to pass the integrated object to the next interation so that it can be evaluated on the basepoint
and to recognise deltas and pass them through as well as evaluating them
-}
integratePieces :: (Enum a, Num a, Eq a, Fractional a) => Pieces a (PolyDelta a) -> Pieces a (PolyDelta a)
integratePieces ps = Pieces (goInt 0 (zero :: (Enum a, Num a, Eq a, Fractional a) => PolyDelta a) (getPieces ps))
    where
        goInt :: (Enum a, Num a, Eq a, Fractional a) => a -> PolyDelta a -> [Piece a (PolyDelta a)] -> [Piece a (PolyDelta a)]
        goInt _ _ [] = [] -- stop when the list of pieces is empty 
        {- we receive the previous object and evaluate it at the basepoint, and then add a constant term so that our
           final integrated piece evaluates to the same value at the basepoint.
           When we have a delta we have to add this to the final value of the previous object, and also pass
           through the delta itself.
        -}
        goInt oldOffset previousObject (x:xs) = case object x of
            P px -> integratedPiece : goInt 0 finalObject xs -- integrating a polynomial
                where
                    bp               = basepoint x
                    basepointValue   = evaluatePD bp previousObject -- the integral at the end of the previous interval
                    integratedObject = integrate (P px)    
                    newValue         = evaluatePD bp integratedObject
                    makeObject y     = P (makePoly y)
                    finalObject      = integratedObject `plus` makeObject (oldOffset + basepointValue - newValue) -- correct the constant
                    integratedPiece  = makePiece (bp, finalObject)
            D dx -> x : goInt newOffset (D dx) xs  -- deltas integrate to themselves
                where
                    bp               = basepoint x
                    basepointValue   = evaluatePD bp previousObject -- the integral at the end of the previous interval
                    newOffset = case previousObject of
                        D _ -> error "Successive delta functions"
                        P _ -> basepointValue -- addition of the delta value will occur in next iteration

evaluateAtApoint :: (Num a, Ord a) => a -> Pieces a (PolyDelta a) -> a
{- |
To evaluate at a point we need to find the interval the point is in and then evaluate the corresponding object
i.e. find i s.t. basepoint(i) <= p < basepoint(i+1).
Our point may be beyond the last basepoint, in which case we take the final value
-}
evaluateAtApoint point as = if point < basepoint (head (getPieces as)) 
                                then error "Invalid point" 
                                else goEval point (getPieces as)
    where
        goEval _ []         = error "Empty piece list"
        goEval _ [_]        = piecesFinalValue as
        goEval p (x1:x2:xs) = if (basepoint x1 <= p) && (p < basepoint x2) 
                                then evaluatePD p (object x1) 
                                else goEval p (x2:xs)

-- | Find the value of the ultimate object at the last basepoint
piecesFinalValue :: Num a => Pieces a (PolyDelta a) -> a
piecesFinalValue (Pieces []) = error "Empty piece list"
piecesFinalValue (Pieces xs) = evaluatePD (basepoint (last xs)) (object (last xs))

{- |
Piecwise convolution requires convolving the pieces pairwise and then summing the results,
i.e. convolve every piece with every other piece and combine the results.
-}
convolvePieces :: (Ord a, Enum a, Fractional a) => Pieces a (PolyDelta a) -> Pieces a (PolyDelta a) -> Pieces a (PolyDelta a)
convolvePieces (Pieces []) _ = error "Empty piece list"
convolvePieces _ (Pieces []) = error "Empty piece list"
convolvePieces as bs = foldr plus zero [Pieces (map makePiece (convolvePolyDeltas a b)) | a <- das, b <- dbs]
    where das = disaggregate (getPieces as)
          dbs = disaggregate (getPieces bs)

-- | disaggregate takes a Pieces list and produces a list of separate bounded intervals
disaggregate :: [Piece a o] -> [(a, a, o)]
disaggregate [] = error "Empty piece list"
disaggregate [_] = [] -- ignore the last piece
disaggregate (x:xs@(x':_)) = (basepoint x, basepoint x', object x) : disaggregate xs

(><) :: (Eq a, Num a) => a -> Pieces a (PolyDelta a) -> Pieces a (PolyDelta a)
-- | multiply by a constant piecewise
infix 7 ><
x >< y = fmap (scalePD x) y

