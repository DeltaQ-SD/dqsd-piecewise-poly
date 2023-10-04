{-# LANGUAGE TypeFamilies, FlexibleInstances, MultiParamTypeClasses #-}

module DeltaQ.PWPs where

import DeltaQ.Class

import PWPs.IRVs

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