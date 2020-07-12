(ns heyarne.thi.ng.utils.poisson
  (:require
   [thi.ng.math.core :as m]
   [thi.ng.geom.vector :as v]
   [thi.ng.geom.core :as g]
   [thi.ng.geom.spatialtree :as tree]))

;; Poisson-Disk-Sampling
;;
;; This namespace implements Poisson-Disk-Sampling in two dimensions with the help of a quadtree.

(defn poisson
  "Given a `container` to find points within, the numbers of samples to generate
  in each step `k` and a minimum distance between to samples `r`, returns a
  vector of samples. When no `seed` is given, starts at a random point
  inside `container`."
  ([container k r]
   (poisson container k r (g/random-point-inside container)))
  ([container k r seed]
   (let [bounds (g/bounds container)
         keep-sample? (partial g/contains-point? bounds)]
     (loop [tree (tree/quadtree bounds)
            active #{seed}]
       (if (empty? active)
         ;; we're done! let's flatten the tree. note that we throw out some of
         ;; our generated samples by calling #(contains-point? container %)
         (persistent! (tree/select-with tree (constantly true) (partial g/contains-point? container)))
         ;; while the active list is not empty, choose a random index from it and
         ;; generate up to `k` within a radius of [r, 2r]; for each point in turn,
         ;; check if it is within distance r of existing samples
         (let [sample (first active)
               candidates (into []
                                (comp
                                 (map #(-> (m/*! % (+ r (m/random r)))
                                           (m/+! sample)))
                                 (filter keep-sample?))
                                (repeatedly k v/randvec2))
               fit (when (seq candidates)
                     (reduce (fn
                               ([_ candidate]
                                (when-not (tree/points-in-circle? tree candidate r)
                                  (reduced candidate))))
                             candidates))]
           (if fit
             ;;  If a point is adequately far from existing samples, emit it as the next sample and add it to the active list
             (recur (g/add-point tree fit fit) (conj active fit))
             ;; if no such point is found, remove the sample from the list of active samples
             (recur tree (disj active sample)))))))))

(comment
  ;; let's run some benchmarks (on my Thinkpad x230 from 2012)

  (require '[criterium.core :as crit])
  (require '[thi.ng.geom.circle :as c])
  (require '[thi.ng.geom.rect :as r])

  (crit/with-progress-reporting
    (crit/bench (poisson (r/rect [-100 -100] [100 100]) 20 10)))
  ;;            Execution time mean : 33.575985 ms
  ;;   Execution time std-deviation : 805.351384 µs
  ;;  Execution time lower quantile : 32.545556 ms ( 2.5%)
  ;;  Execution time upper quantile : 35.375681 ms (97.5%)
  ;;                  Overhead used : 1.869092 ns

  (crit/with-progress-reporting
    (crit/bench (poisson (c/circle [-100 -100] 200) 20 10)))
   ;;           Execution time mean : 170.229276 ms
   ;;  Execution time std-deviation : 16.657318 ms
   ;; Execution time lower quantile : 156.769577 ms ( 2.5%)
   ;; Execution time upper quantile : 225.024192 ms (97.5%)
   ;;                 Overhead used : 1.869092 ns
  )
