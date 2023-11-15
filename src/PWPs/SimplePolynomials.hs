{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE ScopedTypeVariables #-}

{-|
Module      : SimplePolynomials
Description : Polynomials as lists of coefficients
Copyright   : (c) Peter Thompson, 2023
License     : BSD-2-Clause
Maintainer  : peter.thompson@pnsol.com
Stability   : experimental

We define a polynomial over a numeric type as a list of coefficients of that type for increasing powers of the variable.
An empty list is not allowed: a zero polynomical must have at least one (zero) element.
-}
module PWPs.SimplePolynomials
(     Poly (..)
    , makePoly
    , zero 
    , makeMonomial
    , shiftPolyUp
    , scalePoly
    , plus
    , times
    , minus
    , PWPs.ConvolutionClasses.fromInteger 
    , integrate
    , differentiate
    , evaluatePoly
    , convolvePolys
    , compareToZero
) where

import PWPs.ConvolutionClasses
newtype Poly a = Poly [a]
    deriving (Eq,Show,Functor,Foldable)

makePoly :: a -> Poly a
-- | turn a constant into a constant polynomial
makePoly x = Poly [x] 
zeroPoly :: Num a => Poly a
zeroPoly = makePoly 0

degreePoly :: (Num a, Eq a) => Poly a -> Int
-- a constant polynomial has one coefficient and has degree 0; for Euclidian division we want the 
-- degree of the zero polynomial to be negative
degreePoly x = if trimPoly x == zeroPoly then -1 else length (trimPoly x) - 1

trimPoly :: (Num a, Eq a) => Poly a -> Poly a 
-- | remove top zeroes
trimPoly (Poly as) = Poly (reverse $ goTrim $ reverse as)
    where
        goTrim []           = error "Empty polynomial"
        goTrim xss@[_]      = xss -- can't use dropWhile as it would remove the last zero
        goTrim xss@(x:xs)   = if x == 0 then goTrim xs else xss

makeMonomial :: (Eq a, Num a) => Int -> a -> Poly a
-- | put a coefficient in the nth place only
makeMonomial n x = if x == 0 then zeroPoly else Poly (reverse (x:replicate n 0)) 

shiftPolyUp :: (Eq a, Num a) => Poly a -> Poly a
-- | effectively multiply the polynomial by x by shifting all the coefficients up one place.
shiftPolyUp (Poly xs) 
    | xs == [0] = Poly xs     -- don't shift up zero
    | otherwise = Poly (0:xs) 

scalePoly :: Num a => a -> Poly a -> Poly a
-- | scale a polynomial by a constant: more efficient than multiplying by a constant polynomial
scalePoly x (Poly xs) = Poly (map (*x) xs) 

addPolys :: (Eq a, Num a) => Poly a -> Poly a -> Poly a
{- |
   Add polynomials by simplyly adding their coefficients as long as both lists continue.
   When one list runs out we take the tail of the longer list (this prevents us from just using zipWith!).
   Addtion might cancel out the highest order terms, so need to trim just in case.
-}
addPolys (Poly as) (Poly bs) = trimPoly (Poly (go as bs))
    where
        go [] ys = ys
        go xs [] = xs
        go (x:xs) (y:ys) = (x + y) : go xs ys

mulPolys :: (Eq a, Num a) => Poly a -> Poly a -> Poly a
{- |
    multiply term-wise and then add (very simple - FFTs might be faster, but not for today)
    (a0 + a1x + a2x^2 + ...) * (b0 + b1x + b2x^2 ...)
    = a0 * (b0 + b1x + b2x^2 +...) + a1x * (b0 + b1x + ...)
    = (a0*b0) + (a0*b1x) + ...
              + (a1*b0x) +
                         + ...
    (may be an optimisation to be done by getting the shortest poly in the right place)
-}
mulPolys as bs = foldr addPolys zeroPoly (intermediateSums as bs)
    where
        intermediateSums :: (Eq a, Num a) => Poly a -> Poly a -> [Poly a]
        intermediateSums _ (Poly []) = error "Second polynomial was empty"
        intermediateSums (Poly []) _ = [] -- stop when we exhaust the first list
        -- as we consume the coeffecients of the first list, we shift up the second list to increase the power under consideration
        intermediateSums (Poly (x:xs)) ys =  scalePoly x ys : intermediateSums (Poly xs) (shiftPolyUp ys)

integratePoly :: (Eq a, Fractional a) => Poly a -> Poly a
{- |
    Integrate by puting a zero constant term at the bottom and converting ax^n into ax^(n+1)/(n+1).
    0 -> 0x is the first non-constant term, so we start at 1.
    When integrating a zero polynomial with a zero constant we get [0,0] so need to trim
-}
integratePoly (Poly as) = trimPoly (Poly (0:zipWith (/) as (iterate (+1) 1)))

instance (Eq a, Num a) => Num (Poly a) where
    (+)               = addPolys
    (*)               = mulPolys
    negate (Poly a)   = Poly (map negate a)
    abs               = undefined
    signum            = undefined
    fromInteger n     = Poly [Prelude.fromInteger n]

differentiatePoly :: Num a => Poly a -> Poly a 
-- | Simply use dx^n/dx = nx^(n-1)
differentiatePoly (Poly [])     = error "Polynomial was empty"
differentiatePoly (Poly [_])    = zeroPoly -- constant differentiates to zero
differentiatePoly (Poly (_:as)) = Poly (zipWith (*) as (iterate (+1) 1)) -- discard the constant term, everything else noves down one
instance (Eq a, Num a, Fractional a) => Calculable (Poly a)
    where
        plus              = addPolys
        times             = mulPolys
        minus (Poly a)    = Poly (map negate a)
        zero              = zeroPoly
        fromInteger n     = Poly [Prelude.fromInteger n]
        differentiate     = differentiatePoly
        integrate         = integratePoly

evaluatePoly :: Num p => p -> Poly p -> p
{- |
    Evaluate a polynomial at a point.
    Minimise the number of multiplications to evaluate the polynomial by starting from the highest coefficient
    and multiply and add alternately: a0 + a1x + a2x^2 + ... + anx^n = (((anx + an-1)x + an-2)x + ... + a0
-}
evaluatePoly point (Poly as) = foldr (\x acc -> point * acc + x) 0 as

choose :: Int -> Int -> Int
{- |
    Binomial coefficients: simple definition is n `choose` k ~ factorial n `div` (factorial k * factorial (n-k))
    Faster implementation available in math.combinatorics.exact.binomial
-}
n `choose` k
    | k <= 0     = 1
    | k >= n     = 1
    | otherwise = (n-1) `choose` (k - 1) + (n-1) `choose` k -- recursive definition

convolvePolys :: (Fractional a, Eq a, Ord a) => (a, a, Poly a) -> (a, a, Poly a) -> [(a, Poly a)]
-- | Take two polynomials f and g defined on bounded intervals and produce three contiguous pieces as a result
convolvePolys (lf, uf, Poly fs) (lg, ug, Poly gs) 
    | (lf <0) || (lg < 0) = error "Interval bounds cannot be negative"
    | (lf >= uf) || (lg >= ug) = error "Invalid interval" -- upper bounds should be strictly greater than lower bounds
    | (ug - lg) > (uf - lf) = convolvePolys (lg, ug, Poly gs) (lf, uf, Poly fs) -- if g is wider than f, swap the terms
    | otherwise = -- we know g is narrower than f 
        let
            -- sum a set of terms depending on an iterator k (assumed to go down to 0), where each term is a k-dependent
            -- polynomial with a k-dependent multiplier
            sumSeries k mulFactor poly = foldr (\n acc -> acc `plus` (mulFactor n `scalePoly` poly n)) zero [0..k]

            -- the inner summation has a similar structure each time
            innerSum m n term k = sumSeries (m+k+1) innerMult (\j -> makeMonomial (m+n+1-j) (term j)) 
                where
                    innerMult j  = fromIntegral (if even j then (m+k+1) `choose` j else negate ((m+k+1) `choose` j))

            convolemakeMonomials m n innerTerm = sumSeries n (multiplier m n) (innerTerm m n) 
                where
                    multiplier p q k = fromIntegral (if even k then q `choose` k else negate (q `choose` k))/fromIntegral (p+k+1)

            {- 
                For each term, clock through the powers of each polynomial to give convolutions of monomials, which we sum.
                We extract each coefficient of each polynomial, together with an integer recording their position (i.e. power of x), 
                and multiply the coefficients together with the new polynomial generated by convolving the monomials.
            -}
            makeTerm f = foldr plus zero [(a*b) `scalePoly` convolemakeMonomials m n f | (m,a) <- zip [0..] fs, (n,b) <- zip [0..] gs]

            firstTerm  = makeTerm (\m n k -> innerSum m n (lg ^) k `plus` minus (makeMonomial (n-k) (lf^(m+k+1))))

            secondTerm = makeTerm (\m n -> innerSum m n (\k -> lg^k - ug^k))

            thirdTerm  = makeTerm (\m n k -> makeMonomial (n-k) (uf^(m+k+1)) `plus` minus (innerSum m n (ug ^) k))

            -- remove null intervals
            trimTerms []  = []
            trimTerms [x] = [x]
            trimTerms (x:y:xs) = if fst x == fst y then trimTerms (y:xs) else x:trimTerms (y:xs)

        in trimTerms [(0, zero), (lf + lg, firstTerm), (lf + ug, secondTerm), (uf + lg, thirdTerm), (uf + ug, zero)]

{- |
We use Sturm's Theorem to count the number of roots of a polynomial in a given interval.
(See https://en.wikipedia.org/wiki/Sturm%27s_theorem)
Starting from polynomial p, construct the Sturm sequence p0, p1, . . ., where:
p0 = p
p1 = p′
pi+1 = −rem(pi−1, pi) for i > 1
where p′ is the derivative of p and rem(p, q) is the remainder of the Euclidian division of p by q. 
The length of this sequence is at most the degree of p. 
We define V(x) to be the number of sign variations in the sequence of numbers p0(x), p1(x), . . ..
Sturm’s theorem states that, if p is a square-free polynomial (one without repeated roots), then
R(l,r](p) = V (l) − V (r). This extends to non-square-free polynomials provided neither l nor r is a
multiple root of p (a circumstance we shall ignore)

We start from the tuple that emerges from disagregation.
-}
countPolyRoots :: (Fractional a, Eq a, Ord a) => (a, a, Poly a) -> Int
countPolyRoots (l, r, p) = signVariations (sturmSequence l p) - signVariations (sturmSequence r p)
    where
        signVariations :: (Fractional a, Eq a, Ord a) => [a] -> Int
        {-
        When c0, c1, c2, . . . ck is a finite sequence of real numbers, then a sign variation or sign change in the sequence
        is a pair of indices i < j such that cicj < 0, and either j = i + 1 or ck = 0 for all k such that i < k < j
        -}
        signVariations xs = length (filter (< 0) pairsMultiplied)
            where
                zeroesRemoved = filter (/= 0) xs
                -- TODO: deal with all zero corner case
                pairsMultiplied = zipWith (*) zeroesRemoved (tail zeroesRemoved)
        sturmSequence :: (Fractional a, Eq a, Ord a) => a -> Poly a -> [a]
        sturmSequence x q = map (evaluatePoly x) (doSeq [differentiatePoly q, q])
            where
                doSeq :: (Fractional a, Eq a, Ord a) => [Poly a] -> [Poly a]
                {- 
                   note that this is called with a list of length 2 and grows the list, so we don't need to match all cases
                   note that we build this backwards to avoid use of append, but this doesn't affect the number of
                   sign variations so there's no need to reverse it. 
                -}
                doSeq x@(xI:xIminusOne:xs) = if polyRemainder == zero then x else doSeq (negate polyRemainder : x)
                    where
                        polyRemainder = snd (euclidianDivision (xIminusOne, xI))

euclidianDivision :: (Fractional a, Eq a, Ord a) => (Poly a, Poly a) -> (Poly a, Poly a)
{- | 
See https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclidean_division
Take a pair of polynomials a, b, and produce the quotient and remainder q and r s.t. a = bq + r
Input: a and b ≠ 0 two polynomials; Output: q, the quotient, and r, the remainder;
Pseudocode:
    Begin
        q := 0
        r := a
        d := deg(b)
        c := lc(b)
        while deg(r) >= d do
            s := lc(r)/c x^(deg(r)-d)
            q := q + s
            r := r − sb
        end do
        return (q, r)
    end
-}
euclidianDivision (pa, pb) = if pb == zeroPoly then error "Division by zero polynomial" else goDivide (zeroPoly, pa)
    where 
        degB = degreePoly pb
        leadingCoefficient :: Poly a -> a -- coefficient of the highest power term of the poly
        leadingCoefficient (Poly x) = last x
        lcB  = leadingCoefficient pb
        -- goDivide :: (Fractional a, Eq a, Ord a) => (Poly a, Poly a) -> (Poly a, Poly a)
        goDivide (q,r) = if degreePoly r < degB then (q,r) else goDivide (q + s, r - s * pb)
            where s = makeMonomial (degreePoly r - degB) (leadingCoefficient r/lcB)

compareToZero :: (Fractional a, Eq a, Ord a) => (a, a, Poly a) -> Maybe Ordering
{-|
    We measure whether or not a polynomial is consistently above or below zero, or equals zero
-}
compareToZero (l, u, p)
    | l >= u                        = error "Invalid interval"
    | p == zeroPoly                 = Just EQ 
    | lower * upper < 0             = Nothing -- different signs at the two ends
    | countPolyRoots (l, u, p) > 0  = Nothing -- polynomial crosses zero
    | lower > 0                     = Just GT -- upper must also be > 0 fronm previous case
    | otherwise                     = Just LT -- upper and lower both < 0
    where
        lower = evaluatePoly l p
        upper = evaluatePoly u p
