{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# LANGUAGE TypeFamilies, FlexibleInstances, MultiParamTypeClasses #-}

module DeltaQ.PWPs
 ( DeltaQ(..)
 , DeltaQOps(..)
--  , DeltaQ𝛩(..)
--  , DeltaQTimeout(..)
--  , DeltaQUniform(..)
 , shiftedHeaviside
 )
where

import           DeltaQ.Model.DeltaQ

import           PWPs.IRVs (IRV)
import qualified PWPs.IRVs as PWP


instance ProbabilityMass Double where
  type ProbMassModel Double = Double
  toMassModel = id
  fromMassModel = id
  complement x = 1 - x

instance DeltaQ (IRV Double) where
  type ProbMass (IRV Double) = Double
  type Time (IRV Double) = Double

  perfection = PWP.top
  bottom = PWP.bottom

  support = (\(x,y) -> (x, Just y)) . PWP.support

  tangibleMass = PWP.probMass
  


instance DeltaQOps (IRV Double) where
  choice = PWP.probChoice
  nWayChoice = PWP.multiWeightedChoice

  convolve = (PWP.<+>)

  ftf = PWP.firstToFinish
  nWayFtf = PWP.multiFtF

  ltf = PWP.allToFinish
  nWayLtf = PWP.multiAtF

-- instance DeltaQ𝛩 (IRV Double)
{-
instance (Ord a, Num a, Enum a, Fractional a) => ImproperRandomVariable (IRV a) where
    type DelayModel (IRV a) = a
    type ProbabilityModel (IRV a) = a

    diracDelta = constructDelta
    uniform0   = constructUniform

    perfection = diracDelta 0
    bottom     = zeroPDF

    tangibleMass = probMass

instance (Ord a, Num a, Enum a, Fractional a) => Convolvable (IRV a) where
    (<+>) = (PWPs.IRVs.<+>)

instance (Ord a, Num a, Enum a, Fractional a) => NonConcurrentCombination a (IRV a) where
    weightedChoice = probChoice

instance (Ord a, Num a, Enum a, Fractional a) => ConcurrentCombination (IRV a) where
    allToFinish   = PWPs.IRVs.allToFinish
    firstToFinish = PWPs.IRVs.firstToFinish
-}
