#!/usr/bin/env python3
"""Demo 0c, test 6 (d), optional (brief 01_DEMO0_BRIEF.md sec. 6): the lumen of a state file as a shapely union.

Each capsule is the convex hull of its two disks, both drawn as polygons inscribed in the circles with quad_segs 16, 64 and
256. Inscribed polygons lie inside the circles, so the union's area must rise with quad_segs and stay below the exact C
lumen_area, approaching it. Exits 0 with "SKIP" when shapely is not installed: it is not part of the build.

usage: python3 tests/hull_reference.py STATE_FILE
"""
import sys

try:
    from shapely.geometry import Point
    from shapely.ops import unary_union
except ImportError:
    print("SKIP 6(d) hull_reference.py: shapely is not installed")
    sys.exit(0)


def read_state(path):
    disks, area, sec = {}, None, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("["):
                sec = line
                continue
            if sec == "[disks]":
                i, x, y, r, parent, kind, born, alive = line.split()
                disks[int(i)] = (float(x), float(y), float(r), int(parent), int(alive))
            elif sec == "[areas]" and line.startswith("lumen_area"):
                area = float(line.split()[1])
    return disks, area


def main():
    disks, c_area = read_state(sys.argv[1])
    if c_area is None:
        print("FAIL 6(d): the state file has no [areas] (run with --hull=final or all)")
        return 1
    areas = []
    for q in (16, 64, 256):
        caps = []
        for i, (x, y, r, p, alive) in disks.items():
            if p < 0 or not alive or not disks[p][4]:
                continue
            px, py, pr = disks[p][0], disks[p][1], disks[p][2]
            caps.append(unary_union([Point(x, y).buffer(r, quad_segs=q), Point(px, py).buffer(pr, quad_segs=q)]).convex_hull)
        areas.append(unary_union(caps).area)
    rising = areas[0] < areas[1] < areas[2] < c_area
    gaps = [(c_area - a) / c_area for a in areas]
    print("%s 6(d) shapely union, quad_segs 16/64/256: %.6f %.6f %.6f below C %.6f by %.1e %.1e %.1e"
          % ("PASS" if rising else "FAIL", areas[0], areas[1], areas[2], c_area, gaps[0], gaps[1], gaps[2]))
    return 0 if rising else 1


if __name__ == "__main__":
    sys.exit(main())
