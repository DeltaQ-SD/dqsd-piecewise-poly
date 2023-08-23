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
Fundamental Theorem of Calculus and combine appropritely with addition and multiplication.
We also want a convolution operator that behaves correctly.
-}
module PWPs.ConvolutionClasses
(
    Calculable (..)
  , Convolvable (..)
)
where

class Calculable a where
    plus :: a -> a -> a
    times :: a -> a -> a
    minus :: a -> a
    zero :: a
    fromInteger :: Integer -> a
    differentiate :: a -> a
    integrate :: a -> a

class Calculable a => Convolvable a where
    (<+>) :: a -> a -> a -- convolution

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
-}
