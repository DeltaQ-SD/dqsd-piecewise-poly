{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}

{-|
Module      : Piecewise
Description : Sequences of objects defined over contiguous intervals
Copyright   : (c) Peter Thompson, 2023
License     : BSD-2-Clause
Maintainer  : peter.thompson@pnsol.com
Stability   : experimental

We consider 'objects' (which may be polynomials or deltas) defined over numerical intervals, 
assumed to be in increasing order, where we know how to combine objects
when the intervals are identical. When the intervals are not identical we construct the
overlaps so as to produce a combined object over a more complex shared set of intervals.

We assume the sets of intervals have a common initial point, otherwise we could not combine
them over the left-hanging piece. Note this means a piecewise list is never empty.

We allow for zero-width intervals in order to accomodate delta functions. A zero-width
interval should *only* contain a delta, and *must* be followed by a non-zero interval.
A non-zero-width interval should *not* contain a delta.

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
    , monotonic
    , comparePW
) where

import PWPs.ConvolutionClasses

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
instance (Num a, Eq a, Ord a, Calculable b, Evaluable a b) => Calculable (Pieces a b)
    where
        plus x y        = plus <$> x <*> y
        times x y       = times <$> x <*> y
        minus x         = minus <$> x
        zero            = Pieces [Piece {basepoint = 0 :: a, object = zero :: b}]
        fromInteger n   = Pieces [Piece {basepoint = 0 :: a, object = PWPs.ConvolutionClasses.fromInteger n}]
        differentiate   = differentiatePieces
        integrate       = integratePieces

-- | Piecewise differentiation is easy: just differentiate all the objects
differentiatePieces :: (Num a, Eq a, Ord a, Calculable b) => Pieces a b -> Pieces a b
differentiatePieces = fmap differentiate

{- |
Piecewise integration is harder - need to evaluate at the boundary points to make the pieces join up.
We need to pass the integrated object to the next interation so that it can be evaluated on the basepoint
and to recognise deltas and pass them through as well as evaluating them
-}
integratePieces :: (Num a, Eq a, Ord a, Calculable b, Evaluable a b) => Pieces a b -> Pieces a b
integratePieces ps = Pieces (goInt 0 (disaggregate (getPieces ps)))
    where
        goInt :: (Num a, Eq a, Ord a, Calculable b, Evaluable a b) => a -> [(a, a, b)] -> [Piece a b]
        goInt _ [] = [] -- stop when the list of pieces is empty 
        {- 
           We evaluate each integrated object at the initial and final points of the interval.
           We receive the integrated value from the end of the previous interval, and adjust the integrated object so
           that its initial value matches this.
           We add the new object to the ouput list and pass on the received integated value plus the evaluated integral
           at the end of the interval.
        -}
        goInt previousIntegral ((fp,lp,x):xs) = integratedPiece : goInt newIntegral xs
            where
                integratedObject = integrate x    
                basepointValue   = evaluate fp integratedObject -- the integral at the start of the current interval
                finalValue       = evaluate lp integratedObject -- the integral at the end of the current interval
                offset           = previousIntegral - basepointValue
                integratedPiece  = makePiece (fp, boost offset integratedObject) -- correct the constant
                newIntegral      = finalValue + offset

evaluateAtApoint :: (Num a, Ord a, Evaluable a b) => a -> Pieces a b -> a
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
                                then evaluate p (object x1) 
                                else goEval p (x2:xs)

-- | Find the value of the ultimate object at the last basepoint
piecesFinalValue :: (Num a, Evaluable a b) => Pieces a b -> a
piecesFinalValue (Pieces []) = error "Empty piece list"
piecesFinalValue (Pieces xs) = evaluate (basepoint (last xs)) (object (last xs))

instance (Num a, Eq a, Ord a, Evaluable a b) => Evaluable a (Pieces a b) 
    where
        evaluate = evaluateAtApoint
        boost    = fmap . boost
        scale    = (><)

{- |
Piecwise convolution requires convolving the pieces pairwise and then summing the results,
i.e. convolve every piece with every other piece and combine the results.
-}
(<+>) :: (Ord a, Enum a, Fractional a, Calculable b, Evaluable a b, CompactConvolvable a b) => Pieces a b -> Pieces a b -> Pieces a b
infix 7 <+>
(<+>) (Pieces []) _ = error "Empty piece list"
(<+>) _ (Pieces []) = error "Empty piece list"
(<+>) as bs = foldr plus zero [Pieces (map makePiece (convolveIntervals a b)) | a <- das, b <- dbs]
    where das = disaggregate (getPieces as)
          dbs = disaggregate (getPieces bs)

-- | disaggregate takes a Pieces list and produces a list of separate bounded intervals
disaggregate :: [Piece a o] -> [(a, a, o)]
disaggregate [] = error "Empty piece list"
disaggregate [x] = [(basepoint x, basepoint x, object x)] -- turn the last piece into a null interval
disaggregate (x:xs@(x':_)) = (basepoint x, basepoint x', object x) : disaggregate xs

(><) :: (Eq a, Num a, Evaluable a b) => a -> Pieces a b -> Pieces a b
-- | multiply by a constant piecewise
infix 7 ><
(><) = fmap . scale

comparePW :: (Fractional a, Eq a, Ord a, Comparable a b, Calculable b, Evaluable a b) => Pieces a b -> Pieces a b -> Maybe Ordering
-- | Check whether the pieces are all comparable, and if so all compare the same way 
comparePW x y = goCompare (Just EQ) $ disaggregate $ getPieces (plus x (minus y))
    where
        goCompare :: (Fractional a, Eq a, Ord a, Comparable a b) => Maybe Ordering -> [(a, a, b)] -> Maybe Ordering
        goCompare Nothing _     = Nothing           -- stop once we get Nothing
        goCompare prev []       = prev              -- when the list is exhauseted, keep the last result
        goCompare prev (x:xs)
            | prev == Just EQ   = goCompare next xs -- Equality is neutral
            | next == Just EQ   = goCompare prev xs -- Equality is neutral
            | prev == next      = goCompare prev xs -- The two intervals compare the same way
            | otherwise         = Nothing           -- The two intervals compare differently, so stop
            where
                next = compareZero x

