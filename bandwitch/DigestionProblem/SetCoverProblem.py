"""Generic class to solve minimal set cover problems.

The class offers methods useful when the coverage is determined by thresholding
scores (see below).

**Use in other projects**


Note you can use this file in your project as long as the following licence
appears:

::

    Copyright Edinburgh Genome Foundry, 2017
    Licence: MIT
    Author: Zulko


**Technical explanations:**


The problem of enzyme digestion can be generalized and formalized as follows:

A set of parameters Pi
A set of elements ej
A score function score(Pi, ej) => some positive score

The coverage of Pi with threshold T:
C[Pi](T) = [elements e such that score(Pi, e) > T]

The minimal set cover for a threshold T is defined as
M(T) = [Px, Py...] such that [C[Px](T), C[Py](T)...] is minimal set cover of
all elements

Two problems:
A. Find M(T) given T.
B. Find the largest T such that the lenght of M(T) is less than some number N.

Problem A is solved by a standard greedy selection:

.. code::

  greedy_selection(T):
    compute all  C[Pi](T)
    while all elements not covered:
      select a new element.



Problem 2 is solved by bisection on T:

.. code::

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

"""


def maximizing_bisection(f, x_min, x_max, precision=0.001, score=None):
    """Find largest x such that ``score(f(x)) <= 0`` (score(f(x)) increasing).

    Parameters
    ----------

    f
      A function f(x) => anything. The returned value must be either a score,
      or something that can be scored via the ``score`` parameter.

    x_min, x_max
      Range of allowed values for variable x

    precision
      Precision on x (the algorithm stops when x is less than this far away
      from the optimum)

    score
      Scoring function. Leave to none if function f(x) is already a score.



    Returns
    -------

    x, f(x), score(f(x))
      Everything you may need !

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

    while (x_max - x_min) > precision:
        center = (x_min + x_max) / 2.0
        f_center = f(center)
        score_center = score(f_center)
        if score_center > 0:
            x_max = center
        else:
            x_min = center
            f_xmin = f_center
            score_xmin = score_center

    return x_min, f_xmin, score_xmin


def minimal_cover(elements_set, subsets, max_subsets=None, heuristic='default',
                  selected=(), depth=0):
    """
    Parameters
    ----------
    elements_set
      The set of all ements to cover

    subsets
      A list of (name, subset)

    max_subsets
      Maximal number of subsets allowed

    heuristic
      A function ``((name, subset), selected) => value`` where ``name`` is the
      name of a subset, ``subset`` is what remains of the subset at this stage,
      ``selected`` is a list of already-selected subset names.

    selected
      (Recursion parameter, do not use.) Already-selected elements

    depth
      (Recursion parameter, do not use.). Depth of the recursion

    Returns
    -------

      None if no solution was found, else a collection of [(name, subset)...]
      in the order in which the subsets
    """


    if len(elements_set) == 0:
        return []
    if max_subsets == 0:
        return None

    if depth == 0:
        full_set = set().union(*[subset for name, subset in subsets])
        if full_set != elements_set:
            return None

    subsets = [(n, s) for (n, s) in subsets if len(s)]

    def sorting_heuristic(named_subset):
        name, subset = named_subset
        if (heuristic == 'default'):
            return len(subset)
        else:
            return heuristic(named_subset, selected)

    ordered_subsets = sorted(subsets, key=sorting_heuristic)

    while len(ordered_subsets):
        if max_subsets is not None:
            critical_subset_length = len(elements_set) / max_subsets
            max_len = max(len(s) for name, s in ordered_subsets)
            if max_len < critical_subset_length:
                return None
        name, subset = ordered_subsets.pop()
        new_elements_set = elements_set.difference(subset)
        new_subsets = [
            (name_, sub.difference(subset))
            for (name_, sub) in ordered_subsets
        ]
        new_max_subsets = None if (max_subsets is None) else max_subsets - 1
        result = minimal_cover(new_elements_set, new_subsets,
                               heuristic=heuristic,
                               selected=list(selected) + [subset],
                               max_subsets=new_max_subsets,
                               depth=depth + 1)
        if result is not None:
            return result + [name]
        ordered_subsets = [
            subset_
            for (subset_, (new_name, new_subset)) in zip(ordered_subsets,
                                                         new_subsets)
            if len(new_subset) != 0
        ]
    return None


class SetCoverProblem:
    """Generic class for the resolution of score-based set cover problems.

    Parameters
    ----------

    elements
      A set of all elements to be covered

    parameters
      A set of parameters (a parameter and a threshold will determine a subset)

    process_logger
      A ProgLog process logger.

    Note
    ----

    The somwhat complicated math problem solved by the algorithms implemented
    in this class is:

    Given set of parameters (Pi), a set of elements (ej), and a score function
    score(Pi, ej) => some positive score, we call "T-coverage of Pi" the set
    ``C[Pi](T) = (elements e such that score(Pi, e) > T)``

    The minimal set cover M(T) is defined as
    M(T) = [Px, Py...] such that [C[Px](T), C[Py](T)...] is minimal subset
    cover of all elements (ej)

    The class allows to solve the two following problems::
    1. Find M(T) given T.
    2. Find the largest T such that the lenght of M(T) is under some number N.
    """

    def __init__(self, elements, parameters, progress_logger=None):
        self.elements = elements
        self.parameters = parameters
        if progress_logger is None:
            progress_logger = lambda **k: None
            progress_logger.iter_bar = lambda **k: k.popitem()[1]
        self.progress_logger = progress_logger
        self._compute_scores()

    @staticmethod
    def _default_heuristic(named_subset, current_selection):
        name, subset = named_subset
        return len(subset)

    def _parameter_element_score(self, parameter, element):
        raise NotImplementedError('Each subset cover problem must overwrite '
                                  'this function.')

    def _compute_scores(self):
        self.scores = {element: {} for element in self.elements}
        for i, element in enumerate(self.elements):
            self.progress_logger(element_index=i)
            for j, parameter in enumerate(self.parameters):
                score = self._parameter_element_score(parameter, element)
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

    def _compute_coverages(self, threshold):
        return {
            parameter: set(
                element
                for element in self.elements
                if self.scores[element][parameter] >= threshold
            )
            for parameter in self.parameters
        }

    def _select_parameters(self, threshold=None, max_set_size=None,
                        covering_algorithm='greedy', heuristic='default',
                        threshold_tolerance=0.001, bisection=True):
        """Select a subset of parameters which collectively cover all elements.

        This lets you either find the highest-score cover of less than N
        elements, or the best free-sized set of parameters covering elements
        with the highest minimal score.

        Returns
        -------
        min_score, [param1, param2, param3...]
          Where min_score is a value lower than all the score coverages, and
          the list indicates the parameters in the order in which they have
          been selected.


        Parameters
        -----------
        threshold
          If provided and max_set_size is not provided, the method will return
          a free-sized set of parameters covering all elements. If max_set_size
          is provided, the subset will not exceed max_set_size in size. This
          can result in unsolvable.

        max_set_size
          When provided instead of threshold, the method will find the highest
          threshold which allows solutions in max_set_size or less.

        covering_algorithm
          Either 'greedy' for fast approximate solving or 'full' for full
          solving

        heuristic
          Scoring function to prioritize some parameters over others.

        threshold_tolerance
          When max_set_size is provided and not threshold, the returned
          solution will have the optimal threshold, plus or minus this
          tolerance.

        bisection
          For internal recursion. Do not modify.
        """

        if heuristic == 'default':
            heuristic = self._default_heuristic
        if threshold is None:
            threshold = self.bottleneck['score']

        if max_set_size is not None:
            if not bisection:
                return threshold, minimal_cover(
                    elements_set=self.elements,
                    subsets=self._compute_coverages(threshold).items(),
                    max_subsets=max_set_size,
                    heuristic=heuristic
                )

            def select(thr):
                _, selection = self._select_parameters(
                    threshold=thr,
                    covering_algorithm=covering_algorithm,
                    heuristic=heuristic,
                    max_set_size=None if covering_algorithm == 'greedy'
                    else max_set_size,
                    bisection=False
                )
                return selection

            def score(selection):
                if selection is None:
                    return 1e10
                else:
                    return len(selection) - max_set_size
            threshold, selection, _ = maximizing_bisection(
                f=select, x_min=0, x_max=threshold, score=score,
                precision=threshold_tolerance)
            return threshold, selection
        else:
            if covering_algorithm == 'greedy':
                max_subsets = None
            else:
                max_subsets = max_set_size
            return threshold, minimal_cover(
                elements_set=self.elements,
                subsets=self._compute_coverages(threshold).items(),
                max_subsets=max_subsets,
                heuristic=heuristic,
            )

    def min_parameter_score(self, parameter):
        """Return the smallest score on any element for the given parameter."""
        return min([self.scores[e][parameter] for e in self.elements])
