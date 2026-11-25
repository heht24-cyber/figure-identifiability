"""Decode an Escallonia boxplot and find exact feasible sample-mean bounds."""
import argparse
import csv
import json
import math
from collections import defaultdict
from itertools import combinations, combinations_with_replacement
from pathlib import Path

import pymupdf


def extract(pdf):
    """Decode Figure3.pdf with tick and drawing conventions from its R export."""
    with pymupdf.open(pdf) as document:
        page = document[0]
        drawings, words = page.get_drawings(), page.get_text("words")
    white = (1.0, 1.0, 1.0)
    lines = [d["rect"] for d in drawings if d["color"] == white and d["type"] == "s"]
    horizontal = [r.y0 for r in lines if r.width > 500 and r.height < 0.02]
    centers = sorted(r.x0 for r in lines if r.height > 250 and r.width < 0.02)
    ticks = sorted((int(w[4]), min(horizontal, key=lambda y: abs(y - (w[1] + w[3]) / 2)))
                   for w in words if w[2] < 34 and w[4].isdigit())
    if len(ticks) < 3 or not centers:
        raise ValueError("expected vector grid lines and numeric y tick labels")
    pitch = (ticks[-1][1] - ticks[0][1]) / (ticks[-1][0] - ticks[0][0])
    intercept = ticks[0][1] - pitch * ticks[0][0]
    if pitch >= 0 or any(abs(y - intercept - pitch * v) > 0.02 for v, y in ticks):
        raise ValueError("y ticks must agree on a linear calibration")

    def decode(y, step):
        numeric = (y - intercept) / pitch
        snapped = round(numeric / step) * step
        if abs(numeric - snapped) > 0.06:
            raise ValueError(f"coordinate {numeric:.5f} is off the expected lattice")
        return snapped

    by_species = {}
    for x in centers:
        labels = [w[4] for w in words if w[1] > 267 and abs((w[0] + w[2]) / 2 - x) < 1]
        name = "harrisii" if set(labels) == {"New", "species"} else " ".join(labels)
        if not name or name in by_species:
            raise ValueError("missing or duplicate species label")
        near = [d for d in drawings if abs((d["rect"].x0 + d["rect"].x1) / 2 - x) < 0.025]
        boxes = [d["rect"] for d in near if d["fill"] == white and 9 < d["rect"].width < 11]
        medians = [d["rect"] for d in near if (d["width"] or 0) > 2 and 9 < d["rect"].width < 11]
        whiskers = [d["rect"] for d in near if d["type"] == "s" and d["color"] != white
                    and d["rect"].width < 0.02 and d["rect"].y1 < 266]
        outliers = [d["rect"] for d in near if d["type"] == "fs" and d["fill"] != white
                    and [p[0] for p in d["items"]] == ["c"] * 4]
        if len(medians) != 1 or len(boxes) > 1 or len(whiskers) != 2:
            raise ValueError(f"ambiguous box geometry for {name}")
        box = boxes[0] if boxes else medians[0]
        q = (decode(box.y1, 0.25), decode(medians[0].y0, 0.25), decode(box.y0, 0.25))
        w = (decode(max(r.y1 for r in whiskers), 0.25), decode(min(r.y0 for r in whiskers), 0.25))
        o = tuple(sorted(decode((r.y0 + r.y1) / 2, 1) for r in outliers))
        by_species[name] = q, w, o
    return by_species


def summary(values):
    values = sorted(values)
    q = []
    for quartile in (1, 2, 3):
        left, rem = divmod((len(values) - 1) * quartile, 4)
        q.append(values[left] if not rem else
                 ((4 - rem) * values[left] + rem * values[left + 1]) / 4)
    fence = 1.5 * (q[2] - q[0])
    inside, outside = [], []
    for value in values:
        (inside if q[0] - fence <= value <= q[2] + fence else outside).append(value)
    return tuple(q), (min(q[0], *inside), max(q[2], *inside)), tuple(outside)


def extremum(n, box, minimize):
    """Exact optimum sum for ordered integer observations with fixed statistics."""
    q, (low, high), outliers = box
    below = [x for x in outliers if x < low]
    above = [x for x in outliers if x > high]
    if len(below) + len(above) != len(outliers) or len(outliers) >= n:
        raise ValueError("outlier positions inconsistent with whiskers")
    domains = ([(x,) for x in below] + [range(math.ceil(low), math.floor(high) + 1)] * (n - len(outliers))
               + [(x,) for x in above])
    if low < q[0]:
        domains[len(below)] = (int(low),)
    if high > q[2]:
        domains[n - len(above) - 1] = (int(high),)
    fixed, adjacent = {}, defaultdict(list)
    for k, hinge in enumerate(q, 1):
        i, rem = divmod((n - 1) * k, 4)
        if rem:
            adjacent[i + 1].append((rem, round(hinge * 4)))
        else:
            fixed[i] = hinge

    scores, parents = {0: 0}, []
    direction = 1 if minimize else -1
    for i, domain in enumerate(domains):
        new_scores, links = {}, {}
        prior = iter(scores.items())
        entry = next(prior, None)
        best = None
        for x in domain:
            if i in fixed and x != fixed[i]:
                continue
            constraints = adjacent[i]
            if constraints:
                rem, target = constraints[0]
                before, remainder = divmod(target - rem * x, 4 - rem)
                if (remainder or before > x or before not in scores
                        or any((4 - r) * before + r * x != t for r, t in constraints)):
                    continue
            else:
                while entry is not None and entry[0] <= x:
                    if best is None or direction * entry[1] < direction * best[1]:
                        best = entry
                    entry = next(prior, None)
                if best is None:
                    continue
                before = best[0]
            new_scores[x] = scores[before] + x
            links[x] = before
        if not new_scores:
            raise ValueError(f"no feasible boxplot completion at rank {i + 1}")
        scores = new_scores
        parents.append(links)

    last = min(scores, key=lambda x: direction * scores[x])
    sample = []
    for links in reversed(parents):
        sample.append(last)
        last = links[last]
    sample.reverse()
    if len(sample) != n or summary(sample) != box:
        raise ValueError("extreme witness failed its boxplot reconstruction check")
    return sample


def run(pdf, csv_path):
    source = defaultdict(list)
    with csv_path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            if row["NFLOWERS"] not in {"", "NA"}:
                source[row["species"]].append(int(row["NFLOWERS"]))
    boxes = extract(pdf)
    if boxes.keys() != source.keys():
        raise ValueError("PDF and CSV species differ")

    results = {}
    for name, box in boxes.items():
        observed = source[name]
        if summary(observed) != box:
            raise ValueError(f"{name}: PDF={box}, CSV={summary(observed)}")
        n = len(observed)
        smallest, largest = extremum(n, box, True), extremum(n, box, False)
        sums = (sum(smallest), sum(observed), sum(largest))
        results[name] = {"n": n, "quartiles": box[0], "whiskers": box[1], "outliers": box[2],
                         "sums": sums, "means": [x / n for x in sums],
                         "min_sample": smallest, "max_sample": largest}
    reversals = []
    observed_ties = 0
    for (a, x), (b, y) in combinations(results.items(), 2):
        nx, ny = x["n"], y["n"]
        observed = x["sums"][1] * ny - y["sums"][1] * nx
        if not observed:
            observed_ties += 1
            continue
        if observed > 0:
            a, b, x, y, nx, ny = b, a, y, x, ny, nx
        if x["sums"][2] * ny > y["sums"][0] * nx:
            reversals.append([a, b])
    return {"species": results, "observed_ties": observed_ties, "strict_reversals": reversals}


def self_test():
    """Compare the exact optimizer to enumeration on small integer samples."""
    checked = 0
    for n in range(2, 9):
        expect = {}
        for sample in combinations_with_replacement(range(8), n):
            box, total = summary(sample), sum(sample)
            lo, hi = expect.get(box, (total, total))
            expect[box] = min(lo, total), max(hi, total)
        for box, wanted in expect.items():
            found = tuple(sum(extremum(n, box, low)) for low in (True, False))
            if found != wanted:
                raise RuntimeError(f"brute-force mismatch for n={n}, box={box}")
            checked += 1
    return checked


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    result = run(args.pdf, args.csv)
    if args.self_test:
        result["brute_force_cases"] = self_test()
    print(json.dumps(result, indent=2))
