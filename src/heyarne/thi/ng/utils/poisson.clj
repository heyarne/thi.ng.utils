(ns heyarne.thi.ng.utils.poisson
  (:require
   [thi.ng.math.core :as m]
   [thi.ng.geom.vector :as v]
   [thi.ng.geom.core :as g]
   [kdtree :as kdtree]
   [thi.ng.geom.spatialtree :as tree]))

;; Poisson-Disk-Sampling
;;
;; This namespace implements Poisson-Disk-Sampling in two dimensions with the help of a quadtree.

(defn- find-first
  "Returns a reducing function that returns the first item in a collection
  satisifying `pred?`"
  [pred?]
  (completing (fn [_ candidate]
                (when (pred? candidate)
                  (reduced candidate)))))

(defn poisson
  "Given a `container` to find points within, the numbers of samples to generate
  in each step `k` and a minimum distance between to samples `r`, returns a
  vector of samples. When no `seed` is given, starts at a random point
  inside `container`."
  ([container k r]
   (poisson container k r (g/random-point-inside container)))
  ([container k r seed]
   (let [in-container? (partial g/contains-point? container)]
     (loop [tree (tree/quadtree (g/bounds container))
            active #{seed}]
       (if (seq active)
         ;; while the active list is not empty, choose a random index from it and
         ;; generate up to `k` within a radius of [r, 2r]; for each point in turn,
         ;; check if it is within distance r of existing samples
         (let [sample (first active)
               fit (transduce
                    ;; take k samples and keep the ones within our container
                    (comp (take k)
                          (map #(-> (m/*! % (+ r (m/random r)))
                                    (m/+! sample)))
                          (filter in-container?))
                    ;; take the first sample that is not within distance r
                    (find-first #(not (tree/points-in-circle? tree % r)))
                    nil
                    (repeatedly v/randvec2))]
           (if fit
             ;; If a point is adequately far from existing samples, emit it as
             ;; the next sample and add it to the active list
             (recur (g/add-point tree fit fit) (conj active fit))
             ;; if no such point is found, remove the sample from the list of
             ;; active samples
             (recur tree (disj active sample))))
         ;; we're done! let's flatten the tree; we don't need to select anything
         ;; because all points in the tree are restricted to our container anyways
         (persistent! (tree/select-with tree (constantly true) (constantly true))))))))

(defn reducing-first
  ([] {})
  ([acc] (:first acc))
  ([acc cur] (if (contains? acc :first) acc (assoc acc :first cur))))

(defn poisson-kd
  ([container k r]
   (poisson-kd container k r (g/random-point-inside container)))
  ([container k r seed]
   (let [in-container? (partial g/contains-point? container)
         results (transient [])
         r-squared (* r r)]
     (loop [active #{seed}
            tree (kdtree/build-tree [(first active)])]
       (if (seq active)
         ;; while the active list is not empty, choose a random index from it and
         ;; generate up to `k` within a radius of [r, 2r]; for each point in turn,
         ;; check if it is within distance r of existing samples
         (let [sample (first active)
               fit (transduce
                    ;; take k samples and keep the ones within our container
                    (comp (take k)
                          (map #(-> (m/*! % (+ r (m/random r)))
                                    (m/+! sample)))
                          (filter #(try (> (:dist-squared (kdtree/nearest-neighbor tree %)) r-squared)
                                        (catch Exception _
                                          #break (throw _))))
                          (filter in-container?))
                    ;; take the first sample that is not within distance r
                    reducing-first
                    (repeatedly v/randvec2))]
           (if fit
             ;; If a point is adequately far from existing samples, emit it as
             ;; the next sample and add it to the active list
             (do
               (conj! results fit)
               (recur (conj active fit) (kdtree/insert tree fit)))
             ;; if no such point is found, remove the sample from the list of
             ;; active samples
             (recur (disj active sample) tree)))
         ;; we're done! let's flatten the tree; we don't need to select anything
         ;; because all points in the tree are restricted to our container anyways
         (persistent! results))))))

(let [container (r/rect [-100 -100] [100 100])]
  (poisson-kd container 20 10 (g/random-point-inside container)))

(comment
  ;; let's run some benchmarks (on my Thinkpad x230 from 2012)

  (require '[criterium.core :as crit])
  (require '[thi.ng.geom.circle :as c])
  (require '[thi.ng.geom.rect :as r])

  (crit/with-progress-reporting
    (crit/bench (poisson (r/rect [-100 -100] [100 100]) 20 10)))
  ;;            Execution time mean : 25.615302 ms
  ;;   Execution time std-deviation : 1.870445 ms
  ;;  Execution time lower quantile : 24.278897 ms ( 2.5%)
  ;;  Execution time upper quantile : 29.560016 ms (97.5%)
  ;;                  Overhead used : 1.897778 ns

  (crit/with-progress-reporting
    (crit/bench (poisson-kd (r/rect [-100 -100] [100 100]) 20 10)))
  ;;           Execution time mean : 102.236655 ms
  ;;  Execution time std-deviation : 3.393733 ms
  ;; Execution time lower quantile : 96.572241 ms ( 2.5%)
  ;; Execution time upper quantile : 108.823373 ms (97.5%)
  ;;                 Overhead used : 1.854833 ns

  (crit/with-progress-reporting
    (crit/bench (poisson (c/circle [-100 -100] 200) 20 10)))
  ;;            Execution time mean : 99.726122 ms
  ;;   Execution time std-deviation : 4.367676 ms
  ;;  Execution time lower quantile : 94.578971 ms ( 2.5%)
  ;;  Execution time upper quantile : 110.680752 ms (97.5%)
  ;;                  Overhead used : 1.897778 ns

  (crit/with-progress-reporting
    (crit/bench (poisson-kd (c/circle [-100 -100] 200) 20 10)))
  ;;            Execution time mean : 418.034098 ms
  ;;   Execution time std-deviation : 23.716442 ms
  ;;  Execution time lower quantile : 375.606904 ms ( 2.5%)
  ;;  Execution time upper quantile : 463.576173 ms (97.5%)
  ;;                  Overhead used : 1.854833 ns
)
