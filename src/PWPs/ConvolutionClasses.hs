{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-|
Module      : Convolutionclasses
Description : Class definition for operators with testable properties
Copyright   : (c) Peter Thompson, 2023
License     : BSD-2-Clause
Maintainer  : peter.thompson@pnsol.com
Stability   : experimental

The standard algebraic classes such as semiring do not have all the operators we need.
In particular we want integration and differentiation operators that satisfy the
Fundamental Theorem of Calculus and combine appropriately with addition and multiplication.
We also want a convolution operator that behaves correctly.
-}
module PWPs.ConvolutionClasses
(
    Calculable (..)
  , Evaluable (..)
  , CompactConvolvable (..)
  , Comparable (..)
  , Mergeable (..)
)
where

class Calculable a where
    plus          :: a -> a -> a
    times         :: a -> a -> a
    minus         :: a -> a
    zero          :: a
    fromInteger   :: Integer -> a
    differentiate :: a -> a
    integrate     :: a -> a

class Evaluable a b where
    evaluate :: a -> b -> a -- evaluate b at point a
    boost    :: a -> b -> b -- increment b by a
    scale    :: a -> b -> b -- scale b by a

{- |
    Convolution in our library is over finite intervals - this is what piecewiseness needs
-}
class CompactConvolvable a b where
    convolveIntervals :: (a, a, b) -> (a, a, b) -> [(a, b)]  -- convolution

{- |
    We express a partial order by comparing an object on an interval with zero, delivering a Maybe Ordering 
-}
class Comparable a b where
    compareZero :: (a, a, b) -> Maybe Ordering

{- |
    We want to know when two objects can be merged - this is not the same as saying they are equal
-}
class Mergeable a where
    mergeObject :: a -> a -> Maybe a
    zeroObject  :: a

{- |
    Laws:
    Usual stuff with +, *, -
    differentiate . integrate = id              } Fundamental theorem
    integrate . differentiate = add constant    } of calculus
    times distributes over integration and differentiation
    addition distributes over everything 
    convolution is commutative and associative
    differentiate (f <+> g) == (differentiate f) <+> g == f <+> (differentiate g)
    integrate (f <+> g) == (integrate f) * (integrate g)
    operations maintain ordering, if present
-}
