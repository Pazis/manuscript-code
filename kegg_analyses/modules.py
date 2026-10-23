"""KEGG module definition parser and evaluator.

Grammar, loosest to tightest binding:
  seq  := alt (' ' alt)*            steps, all required
  alt  := comp (',' comp)*          alternatives, one required
  comp := ['-'] unit (('+'|'-') unit)*   complex; '-' marks an optional component
  unit := Kxxxxx | Mxxxxx | '(' seq ')'
A module's top-level seq elements are its blocks (KEGG Mapper convention).
"""
import re

TOKEN = re.compile(r"K\d{5}|M\d{5}|[()+\-, ]")


def parse(definition: str):
    tokens = TOKEN.findall(re.sub(r"\s+", " ", definition.strip()))
    if "".join(tokens) != re.sub(r"\s+", " ", definition.strip()):
        raise ValueError(f"unparseable definition: {definition!r}")
    pos = 0

    def peek():
        return tokens[pos] if pos < len(tokens) else None

    def take(expected=None):
        nonlocal pos
        tok = peek()
        if tok is None or (expected and tok != expected):
            raise ValueError(f"expected {expected!r} at token {pos} of {definition!r}, got {tok!r}")
        pos += 1
        return tok

    def seq():
        items = [alt()]
        while peek() == " ":
            take()
            items.append(alt())
        return ("seq", items)

    def alt():
        items = [comp()]
        while peek() == ",":
            take()
            items.append(comp())
        return ("alt", items)

    def comp():
        parts = [(take("-") == "-", unit())] if peek() == "-" else [(False, unit())]
        while peek() in ("+", "-"):
            parts.append((take() == "-", unit()))
        return ("comp", parts)

    def unit():
        if peek() == "(":
            take("(")
            node = seq()
            take(")")
            return node
        return ("id", take())

    tree = seq()
    if pos != len(tokens):
        raise ValueError(f"trailing tokens in {definition!r}")
    return tree


def satisfied(node, kos: set, modules: dict) -> bool:
    kind, body = node
    if kind == "id":
        return modules.get(body, False) if body.startswith("M") else body in kos
    if kind == "seq":
        return all(satisfied(n, kos, modules) for n in body)
    if kind == "alt":
        return any(satisfied(n, kos, modules) for n in body)
    return all(satisfied(n, kos, modules) for optional, n in body if not optional)


def evaluate(tree, kos: set, modules: dict | None = None) -> dict:
    """Block-level status of one module for a KO set; `modules` holds completeness of referenced modules."""
    blocks = tree[1]
    present = sum(satisfied(b, kos, modules or {}) for b in blocks)
    missing = len(blocks) - present
    status = {0: "complete", 1: "1 block missing", 2: "2 blocks missing"}.get(missing, "incomplete")
    return {"steps_total": len(blocks), "steps_present": present,
            "blocks_missing": missing, "complete": missing == 0, "status": status}


def evaluate_panel(trees: dict, kos: set) -> dict:
    """Evaluate every module; modules referencing other modules are evaluated after them."""
    results, pending = {}, dict(trees)
    while pending:
        ready = [m for m, t in pending.items() if not (refs(t) - results.keys())]
        if not ready:
            raise ValueError(f"unresolvable module references: {sorted(pending)}")
        for m in ready:
            results[m] = evaluate(pending.pop(m), kos, {k: v["complete"] for k, v in results.items()})
    return results


def refs(node) -> set:
    kind, body = node
    if kind == "id":
        return {body} if body.startswith("M") else set()
    return set().union(*(refs(n if kind != "comp" else n[1]) for n in body))


def kos_in(node) -> set:
    kind, body = node
    if kind == "id":
        return {body} if body.startswith("K") else set()
    return set().union(*(kos_in(n if kind != "comp" else n[1]) for n in body))
