"""Residue-overlap evaluation and paired target-level uncertainty."""
from math import comb, sqrt
import numpy as np

MODELS = ("p2rank", "learned", "learned_plm")
THRESHOLDS = tuple(round(i / 100, 2) for i in range(10, 55, 5))

def iou(predicted, reference):
    predicted, reference = set(predicted), set(reference)
    union = predicted | reference
    return len(predicted & reference) / len(union) if union else 0.0

def top_hit(pockets, reference, k=5, threshold=0.25):
    return any(p["rank"] <= k and iou(p["residues"], reference) >= threshold for p in pockets)

def wilson(hits, n):
    if not n:
        raise ValueError("The denominator must be positive")
    z = 1.959963984540054
    p = hits / n
    denominator = 1 + z*z/n
    center = (p + z*z/(2*n)) / denominator
    half = z * sqrt(p*(1-p)/n + z*z/(4*n*n)) / denominator
    return [center-half, center+half]

def paired_comparison(a, b, seed=20261023, resamples=10000):
    a, b = np.asarray(a, dtype=int), np.asarray(b, dtype=int)
    if len(a) != len(b) or not len(a):
        raise ValueError("Paired, nonempty target arrays are required")
    delta = a-b
    boot = delta[np.random.default_rng(seed).integers(0, len(a), size=(resamples, len(a)))].mean(axis=1)
    gained = int(((a == 1) & (b == 0)).sum())
    lost = int(((a == 0) & (b == 1)).sum())
    discordant = gained + lost
    p = min(1.0, 2 * sum(comb(discordant, i) for i in range(min(gained, lost)+1)) / 2**discordant) if discordant else 1.0
    return {"n":len(a), "a_hits":int(a.sum()), "b_hits":int(b.sum()), "a_only":gained,
            "b_only":lost, "both":int(((a == 1)&(b == 1)).sum()), "neither":int(((a == 0)&(b == 0)).sum()),
            "difference":float(delta.mean()), "difference_ci":np.quantile(boot,[.025,.975]).tolist(), "mcnemar_p":p}

def summarize(rows,models=MODELS):
    n = len(rows)
    if not n:
        raise ValueError("No targets to evaluate")
    totals = {}
    for model in models:
        totals[model] = {}
        for k in range(1,11):
            hits = sum(r["metrics"][model]["top_hits"][str(k)] for r in rows)
            totals[model][str(k)] = {"hits":hits,"n":n,"ci":wilson(hits,n)}
    comparisons = []
    for a,b in [("learned_plm","p2rank"),("learned_plm","learned"),("learned","p2rank")]:
        if a not in models or b not in models:continue
        for k in [1,3,5]:
            comparison = paired_comparison([r["metrics"][a]["top_hits"][str(k)] for r in rows],
                                           [r["metrics"][b]["top_hits"][str(k)] for r in rows])
            comparisons.append({"a":a,"b":b,"k":k,**comparison})
    primary = sorted([c for c in comparisons if c["a"]=="learned_plm" and c["k"]==5],key=lambda c:c["mcnemar_p"])
    previous=0
    for i,c in enumerate(primary):
        c["holm_adjusted_p"]=min(1.0,max(previous,(len(primary)-i)*c["mcnemar_p"]))
        previous=c["holm_adjusted_p"]
    return {"n":n,"totals":totals,"comparisons":comparisons,
            "thresholds":{f"{threshold:.2f}":{m:sum(r["metrics"][m]["threshold_hits"][f"{threshold:.2f}"] for r in rows) for m in models} for threshold in THRESHOLDS}}

def score_target(predictions, reference, present):
    reference, present = set(reference), set(present)
    if not reference or not reference <= present:
        raise ValueError("The nonempty reference must map to observed input residues")
    ranks=[p["rank"] for p in predictions]
    if any(not isinstance(rank,int) or rank < 1 for rank in ranks) or len(set(ranks)) != len(ranks):
        raise ValueError("Pocket ranks must be unique positive integers")
    for p in predictions:
        if not set(p["residues"]) <= present:
            raise ValueError("A predicted residue is absent from the input")
    scores=[(p["rank"],iou(p["residues"],reference)) for p in predictions]
    return {"pockets":len(predictions),"first_correct_rank":min((rank for rank,s in scores if s >= .25),default=None),
            "oracle_any_hit":any(s >= .25 for _,s in scores),
            "top_hits":{str(k):any(rank<=k and s>=.25 for rank,s in scores) for k in range(1,11)},
            "top5_best_iou":max((s for rank,s in scores if rank<=5),default=0),
            "threshold_hits":{f"{t:.2f}":any(rank<=5 and s>=t for rank,s in scores) for t in THRESHOLDS}}
