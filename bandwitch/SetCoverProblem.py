"""

Note you can use this file in your project as long as the following licence
appears:

Copyright Edinburgh Genome Foundry, 2017
Licence: MIT
Author: Zulko

Technical explanations:

The problem of enzyme digestion can be generalized and formalized as follows:

A set of parameters Pi
A set of elements ej
A score function score(Pi, ej) => some positive score

The coverage of Pi with threshold T:
C[Pi](T) = [elements e such that score(Pi, e) > T]

The minimal set cover for T is defined as
M(T) = [Px, Py...] such that [C[Px](T), C[Py](T)...] is minimal set cover of
all elements

Two problems:
A. Find M(T) given T.
B. Find the largest T such that the lenght of M(T) is less than some number N.

Problem A is solved by a standard greedy selection:

```
greedy_selection(T):
 compute all  C[Pi](T)
 while all elements not covered:
     select a new element.
```


Problem 2 is solved by bisection on T:

```
T_minimax = min_Pi(max_ej score(Pi, e_j))
Tmin, Tmax = 0, T_minimax
while (Tmax - Tmin) > tolerance:
    Tcenter = (Tmin + Tmax) / 2
    M_center = greedy_selection(Tcenter)
    if length(M_center) > N:
        Tmax = Tcenter
    else:
        Tmin = Tcenter
return greedy_selection(Tmin)

```
"""
def maximizing_bisection(f, x_min, x_max, tolerance=0.001, score=None):
    """Find the largest x with score(f(x)) <= 0. where score(f(x)) increasing.

    Return x, f(x), score(f(x))

    """
    if score is None:
        def score(value):
            return value

    f_xmin = f(x_min)
    score_xmin = score(f_xmin)
    if score_xmin > 0:
        return None, None, None
    f_xmax = f(x_max)
    score_xmax = score(f_xmax)
    if score_xmax <= 0:
        return x_max, f_xmax, score_xmax

    while (x_max - x_min) / 2.0 > tolerance:
        center = (x_min + x_max) / 2.000
        f_center = f(center)
        score_center = score(f_center)
        if score_center > 0:
            x_max = center
        else:
            x_min = center
            f_xmin = f_center
            score_xmin = score_center

    return x_min, f_xmin, score_xmin

def greedy_minimal_set_cover(coverages, elements=None, heuristic=None):
    """Return a (hopefully) 'minimal' full-covering set, with a greedy method.

    The greedy method consists in selecting first the element with the biggest
    coverage, then the element with the biggest coverage among yet-uncovered
    targets, etc. until all targets are covered.

    Parameters
    ----------

    coverages
      A dictionary ``{set_name: [covered elements]}``

    elements
      The full set of elements to be covered. Providing this elements enables
      to quickly check that there is a solution to the problem, i.e. that
      the union of all coverages from all elements is the full set.

    heuristic
      Function ``(element) => score`` to select the element with the highest
      score at each iteration. By default, the heuristic is the length of the
      element's coverage of yet-uncovered targets.
    """
    current_selection = []
    if heuristic is None:
        def heuristic(element, coverage, current_selection):
            return len(coverage[element])

    coverages = {k: set(v) for k, v in coverages.items()}
    if elements is not None:
        full_coverage_set = set().union(*coverages.values())
        if full_coverage_set != elements:
            #print ("AAAAAAH", elements, len(full_coverage_set))
            return None

    def key(e):
        """Function with regard to which the next best element is selected."""
        return heuristic(e, coverages, current_selection)

    while len(coverages) > 0:
        selected = max(coverages, key=key)
        covered = coverages.pop(selected)
        if len(covered) == 0:
            break
        current_selection.append(selected)
        for element, coverage in coverages.items():
            coverages[element] = coverage.difference(covered)
    return current_selection


class SetCoverProblem:

    defaults = {
        'covering_algorithm': greedy_minimal_set_cover,
        'heuristic': lambda el, coverages, selection: len(coverages[el])
    }

    def __init__(self, elements, parameters, progress_logger=None):
        self.elements = elements
        self.parameters = parameters
        if progress_logger is None:
            def progress_logger(**k):
                pass
        self.progress_logger = progress_logger
        self.compute_scores()

    @staticmethod
    def default_heuristic(element, coverages, current_selection):
        return len(coverages[element])

    def parameter_element_score(self, parameter, element):
        pass

    def compute_scores(self):
        self.scores = {element: {} for element in self.elements}
        for i, element in enumerate(self.elements):
            self.progress_logger(element_index=i)
            for j, parameter in enumerate(self.parameters):
                score = self.parameter_element_score(parameter, element)
                self.scores[element][parameter] = score

        for parameters_dict in self.scores.values():
            parameters_dict["MAX_SCORE"] = max(
                (score, parameter)
                for (parameter, score) in parameters_dict.items()
            )
        ((min_max_score, parameter), element) = min(
            (parameters_dict["MAX_SCORE"], element)
            for element, parameters_dict in self.scores.items()
        )
        self.bottleneck = dict(score=min_max_score, parameter=parameter,
                               element=element)

    def compute_coverages(self, threshold):
        return {
            parameter: set(
                element
                for element in self.elements
                if self.scores[element][parameter] >= threshold
            )
            for parameter in self.parameters
        }

    def select_elements(self, threshold=None, max_elements=None,
                        covering_algorithm ='default', heuristic='default',
                        threshold_tolerance=0.001):
        if heuristic == 'default':
            heuristic = self.defaults['heuristic']
        if covering_algorithm == 'default':
            covering_algorithm = self.defaults['covering_algorithm']
        if threshold is None:
            threshold = self.bottleneck['score']

        if max_elements is not None:
            def select(thr):
                _, selection = self.select_elements(
                    threshold=thr,
                    covering_algorithm=covering_algorithm,
                    heuristic=heuristic
                )
                return selection

            def score(selection):
                if selection is None:
                    return 1e10
                else:
                    return len(selection) - max_elements
            threshold, selection, _ = maximizing_bisection(
               f=select, x_min=0, x_max=threshold, score=score,
               tolerance=threshold_tolerance)
            return threshold, selection
        else:
            return threshold, covering_algorithm(
                coverages=self.compute_coverages(threshold),
                elements=self.elements,
                heuristic=heuristic
            )
