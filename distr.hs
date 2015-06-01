import Data.Ord (comparing)
import Data.List (sortBy)
import Data.Monoid

data Discrete a = Discr [(a,Float)] deriving (Show)

unDiscr (Discr xs) = xs

sortedDiscr :: Discrete a -> [(a,Float)]
sortedDiscr = reverse . sortBy (comparing snd) . unDiscr

sortDescendingProb :: Discrete a -> Discrete a
sortDescendingProb = Discr . sortedDiscr

uniformD :: [a] -> Discrete a
uniformD [] = Discr []
uniformD xs  = Discr $ zip xs (repeat p)
     where p = 1.0 / fromIntegral ( length xs )

normalWithStddev :: Float -> [a] -> Discrete a
normalWithStddev sigm xs = Discr $ zip xs ks
     where ks = normalDiscreteCoefficients sigm (length xs)

normal1Sigma = normalWithStddev 1.0
normal2Sigma = normalWithStddev 2.0
normal3Sigma = normalWithStddev 3.0

estimVar :: Fractional a => Discrete a -> (a,a)
estimVar d = (estimation d,variance d)

estimation :: Fractional a => Discrete a -> a
estimation (Discr xs) = (sum $ map (\(p,q) -> p*realToFrac q) xs ) * recip ( realToFrac $ sum $ map snd xs)

variance :: Fractional a => Discrete a -> a
variance d = estimation $ transform (\x -> square (e-x)) d
    where e = estimation d
          square t = t*t

transform :: (a->b) -> Discrete a -> Discrete b
transform = fmap

instance Functor Discrete where
    fmap f (Discr xs) = Discr $ map g xs
        where g (p,q) = (f p,q)

instance Applicative Discrete where
    pure x = uniformD [x]
    Discr fs <*> Discr xs = Discr $ [ (f x,p*q) | (f,p) <-fs, (x,q) <-xs ]

instance Monad Discrete where
    -- return = pure -- x = uniformD [x]
    Discr xs >>= f = Discr $ [ (y,p*q)  | (x,p) <- xs, (y,q) <- (unDiscr.f) x ]

instance Monoid a => Monoid (Discrete a) where
    mempty = uniformD [mempty]
    Discr xs `mappend` Discr ys = Discr [ (x `mappend` y,p*q) | (x,p) <- xs, (y,q) <- ys ]

compressD :: Ord a => Discrete a -> Discrete a
compressD (Discr xs) = Discr $ compr xs

compressSortDesc :: Ord a => Discrete a -> Discrete a
compressSortDesc = sortDescendingProb . compressD

-- cov (X,Y) = E[(X-E(X))*(Y-E(Y))]
covarianceByDef :: Fractional a => Discrete (a,a) -> a
-- covarianceByDef :: ( Fractional a, Fractional b) => Discrete (a,b) -> a
covarianceByDef d = estimation $ transform (\(x,y) -> (x-ex)*(y-ey) ) d
      where ex = et fst -- estimation $ transform fst d
            ey = et snd -- estimation $ transform snd d
            et f = estimation $ transform f d

-- cov (X,Y) = E(X*Y) - E(X)*E(Y)
covarianceByFormula :: Fractional a => Discrete (a,a) -> a
covarianceByFormula d = et (\(x,y) -> x*y ) - et fst * et snd
      where et f = estimation $ transform f d

covariance :: Fractional a => Discrete (a,a) -> a
covariance = covarianceByDef

-- correlation: rho(X,Y) = cov(X,Y)/sqrt(V(X)*V(Y))
--  -1 <= rho <= 1

correlation2 :: (Eq a, Fractional a) => Discrete (a,a) -> a
correlation2 d = if ((vx == 0) || (vy == 0))
                   then 0
                   else sqr(covariance d) / (vx * vy)
     where vx = variance $ transform fst d
           vy = variance $ transform snd d
           sqr = (\x -> x*x)

-- -----------------------------------

singleStep :: Num a => a  -> Discrete a
singleStep x = uniformD $ map (x+) [(-1),0,1]








-- -----------------------------------
-- helper functions below

-- not normalized!
-- normalStandardDensity x = exp( -x*x /2) -- /sqrt(2*pi)

normalDiscreteCoefficients :: Float -> Int -> [Float]
normalDiscreteCoefficients _     0 = []
normalDiscreteCoefficients _     1 = [1.0]
normalDiscreteCoefficients sigma n = normalizeTo1 $ map normalStandardDensity ps
    where ps = map (\m -> sigma*(fromIntegral m / k2 - 1) ) [0..n-1]
          k2 = fromIntegral(n - 1) / 2.0
          normalStandardDensity x = exp( -x*x /(2 {- *sigma*sigma -} ) )
-- "sigma" is "how many sigmas are in given distribution"


normalizeTo1 :: Fractional a => [a] -> [a]
normalizeTo1 xs = map (*k) xs
    where k = 1.0 / (foldr (+) 0 xs)

compr :: (Ord a, Num b) => [(a,b)] -> [(a,b)]
-- compactify equal items (and sort)
compr = compr' . sortBy (comparing fst)
   where
         eq' :: Ord a => (a,b) -> (a,b) -> Bool
         eq' (p,_) (q,_) = (p == q)
         add' (p,m1) (q,m2) = (p,m1+m2)
         compr' [] = []
         compr' [x] = [x]
         compr' (x:y:xs) = if (eq' x y) then (compr'(add' x y : xs)) else (x : compr' (y:xs))
         compr'' = filter ((/=0.0).snd) -- throw away zero components



-- examples:
-- estimation $ liftM2 (*) (normal3Sigma [1.0..5] )  (normal3Sigma [1.0..7]) == 4.00

