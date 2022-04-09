predicate dom_disjunctive_var(array[int] of var int: s,
  array[int] of var int: d);
predicate dom_disjunctive(array[int] of var int: s,
  array[int] of int: d);

predicate fzn_disjunctive(array[int] of var int: s, array[int] of var int: d) =
  dom_disjunctive_var(s, d);
predicate fzn_disjunctive(array[int] of var int: s, array[int] of int: d) =
  dom_disjunctive(s, d);
