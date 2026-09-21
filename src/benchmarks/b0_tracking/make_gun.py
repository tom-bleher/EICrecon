#!/usr/bin/env python3
"""Deterministic synthetic B0 proton guns, in GeV/mm (not a beam-background model)."""
import argparse
import math
import random
from pathlib import Path


def generate(path, events, multiplicity, momentum, tx, ty, seed):
    rng = random.Random(seed)
    angle = -0.025
    mass = 0.93827208816
    with Path(path).open('w') as out:
        out.write('HepMC::Version 3.02.05\nHepMC::Asciiv3-START_EVENT_LISTING\n')
        for event in range(events):
            out.write(f'E {event} {multiplicity} {2*multiplicity}\nU GEV MM\n')
            for i in range(multiplicity):
                p, x, y = momentum, tx, ty
                if multiplicity > 1:
                    p = rng.uniform(8, 41)
                    theta = math.acos(rng.uniform(math.cos(.022), math.cos(.004)))
                    phi = rng.uniform(-math.pi, math.pi)
                    x, y = math.tan(theta)*math.cos(phi), math.tan(theta)*math.sin(phi)
                z = p / math.sqrt(1+x*x+y*y)
                px, py, pz = z*(x*math.cos(angle)+math.sin(angle)), z*y, z*(-x*math.sin(angle)+math.cos(angle))
                energy = math.hypot(p, mass)
                # Each independent gun proton has its own incoming parent and
                # origin vertex. Multiple guns are superposed in one event.
                parent = 2*i+1
                out.write(f'P {parent} 0 2212 {px:.17e} {py:.17e} {pz:.17e} {energy:.17e} {mass:.17e} 4\n')
                out.write(f'P {parent+1} {parent} 2212 {px:.17e} {py:.17e} {pz:.17e} {energy:.17e} {mass:.17e} 1\n')
        out.write('HepMC::Asciiv3-END_EVENT_LISTING\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--events', type=int, default=100)
    parser.add_argument('--multiplicity', type=int, default=1)
    parser.add_argument('--momentum', type=float, default=8)
    parser.add_argument('--tx', type=float, default=-.002)
    parser.add_argument('--ty', type=float, default=.012)
    parser.add_argument('--seed', type=int, default=170008)
    args = parser.parse_args()
    if args.events < 1 or args.multiplicity < 1 or args.momentum <= 0:
        parser.error('events, multiplicity and momentum must be positive')
    generate(args.output, args.events, args.multiplicity, args.momentum, args.tx, args.ty, args.seed)
