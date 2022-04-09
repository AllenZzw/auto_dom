predicate dom_cumulative_var(array[int] of var int: s,
  array[int] of var int: d,
  array[int] of var int: r,
  var int: b) 
  = let { set of int: tasks = index_set(s), 
    set of int: times = min([lb(s[i]) | i in tasks]) .. max([ub(s[i]) + ub(d[i])-1 | i in tasks]) } in 
    forall( t in times )( b >= sum( i in tasks ) (bool2int( s[i] <= t /\ t < s[i] + d[i] ) * r[i]) );
  
predicate dom_cumulative(array[int] of var int: s,
  array[int] of int: d,
  array[int] of int: r,
  int: b) 
  = let { set of int: tasks = index_set(s), 
    set of int: times = min([lb(s[i]) | i in tasks]) .. max([ub(s[i]) + ub(d[i])-1 | i in tasks]) } in 
    forall( t in times )( b >= sum( i in tasks ) (bool2int( s[i] <= t /\ t < s[i] + d[i] ) * r[i]) );


predicate fzn_cumulative(array[int] of var int: s,
  array[int] of var int: d,
  array[int] of var int: r,
  var int: b) =
    dom_cumulative_var(s, d, r, b);

%% TODO: once disjunctive propagator is fixed, add
%% special case for b = 1.
predicate fzn_cumulative(array[int] of var int: s,
  array[int] of int: d,
  array[int] of int: r,
  int: b) =
    dom_cumulative(s, d, r, b);
