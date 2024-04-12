{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# LANGUAGE TypeFamilies, FlexibleInstances, MultiParamTypeClasses #-}

module DeltaQ.PWPs
 (
   IRV
 , DeltaQ(..)
 , DeltaQOps(..)
 , DeltaQ𝛩(..)
 , DeltaQTimeout(..)
 , DeltaQUniform(..)
 , DeltaQIntrospection(..)
 , Slazard(..)
 , DeltaQVisualisation(..)
 , shiftedHeaviside
 )
where

import           DeltaQ.Model.DeltaQ
import           DeltaQ.Model.Introspection
import           DeltaQ.Model.Utilities

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

  cumulativeMass = PWP.cumulativeMass

  cumulativeMass' a b = Just $ PWP.cumulativeMass a b

  centiles = flip PWP.centiles

instance DeltaQOps (IRV Double) where
  choice = PWP.probChoice
  nWayChoice = PWP.multiWeightedChoice

  convolve = (PWP.<+>)

  ftf = PWP.firstToFinish
  nWayFtf = PWP.multiFtF

  ltf = PWP.allToFinish
  nWayLtf = PWP.multiAtF

instance DeltaQ𝛩 (IRV Double) where
 delay = flip PWP.shiftIRV

instance DeltaQUniform (IRV Double) where
  uniform0 = PWP.constructUniform

instance DeltaQVisualisation (IRV Double) where
  asDiscreteCDF = flip PWP.displayCDF

  asDiscretePDF irv n = (\x -> [Right x]) $ (PWP.displayPDF  n irv)
