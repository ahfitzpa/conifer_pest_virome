#!/usr/bin/env python3
"""Tree helpers shared by the host-range scripts: Fitch parsimony and edges.

Loaded by path with importlib, because the callers run under `python3 -I`, which
drops the script's own directory from sys.path and so breaks a plain import.

Why not the draft's `felsenstein_anc`: despite the name it is neither Felsenstein
nor correct. It multiplies the children's state dictionaries, so a state survives
only if every child carries it; when children disagree, every product is 0, and
`max()` then returns whichever zero-valued key the dict happens to yield first.
That is arbitrary and order-dependent, precisely at the nodes where a host
changes -- the nodes the whole analysis is about. It also treats an unknown host
as a state, so one unknown child can zero out an otherwise clear node.

Fitch parsimony is used instead: intersection where children agree, union where
they do not, with unknowns ignored rather than counted as a state. Ties are
reported as ambiguous instead of silently resolved.
"""
from __future__ import annotations


def children(node):
    return list(getattr(node, "clades", None) or [])


def is_terminal(node) -> bool:
    return not children(node)


def edges(root):
    """(parent, child) for every edge, parents before children."""
    out = []
    stack = [root]
    while stack:
        node = stack.pop()
        for child in children(node):
            out.append((node, child))
            stack.append(child)
    return out


def terminals(node):
    if is_terminal(node):
        return [node]
    out = []
    for child in children(node):
        out.extend(terminals(child))
    return out


def node_label(node, index: dict) -> str:
    """A stable name for a node: its own, else nodeN by traversal order."""
    name = getattr(node, "name", None)
    if name:
        return str(name)
    return index.setdefault(id(node), f"node{len(index) + 1}")


def fitch(root, states: dict[str, set[str]], unknown_as_missing: bool = True):
    """Fitch parsimony over discrete states.

    states maps a terminal's name to its possible states; an empty set means
    unknown. Returns (assigned, ambiguous, changes) where
      assigned  {node id: frozenset of states at that node}
      ambiguous {node id: True} where more than one state remains
      changes   number of state changes counted on the downward pass
    """
    up: dict[int, frozenset] = {}
    changes = 0

    def down(node):
        nonlocal changes
        if is_terminal(node):
            name = getattr(node, "name", "") or ""
            s = states.get(name, set())
            up[id(node)] = frozenset(s)
            return up[id(node)]
        sets = []
        for child in children(node):
            cs = down(child)
            # 1: an unknown child contributes nothing rather than wiping the node
            if cs or not unknown_as_missing:
                sets.append(cs)
        if not sets:
            up[id(node)] = frozenset()
            return up[id(node)]
        inter = set(sets[0])
        for s in sets[1:]:
            inter &= s
        if inter:
            up[id(node)] = frozenset(inter)
        else:
            # 2: children disagree -> union, and that is one change
            union = set()
            for s in sets:
                union |= s
            up[id(node)] = frozenset(union)
            changes += 1
        return up[id(node)]

    down(root)

    # 3: resolve downward, preferring the parent's state where it is allowed,
    # so a node is only ambiguous when the data genuinely does not decide it
    assigned: dict[int, frozenset] = {}
    ambiguous: dict[int, bool] = {}

    def resolve(node, parent_state):
        s = up[id(node)]
        if parent_state and (s & parent_state):
            s = frozenset(s & parent_state)
        assigned[id(node)] = s
        ambiguous[id(node)] = len(s) > 1
        for child in children(node):
            resolve(child, s)

    resolve(root, frozenset())
    return assigned, ambiguous, changes


def single_state(states: frozenset) -> str:
    """The one state, or '' when unknown or ambiguous."""
    return next(iter(states)) if len(states) == 1 else ""
