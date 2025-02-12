#!/usr/bin/env python3

import csv, operator
from itertools import chain
from functools import partial

with open('./output/1j5nenergy.csv') as original, open('./output/1j5nenergy.csv') as new:
    readers = map(lambda x: csv.reader(x, quoting=csv.QUOTE_NONNUMERIC), [original, new])
    diffs = map(partial(map, operator.sub), *readers)
    max_diff = max(map(operator.abs, chain.from_iterable(diffs)))

print('linear',max_diff)

if max_diff > 1:
    exit(1)

with open('./output/mitoenergy.csv') as original, open('./output/mitoenergy.csv') as new:
    readers = map(lambda x: csv.reader(x, quoting=csv.QUOTE_NONNUMERIC), [original, new])
    diffs = map(partial(map, operator.sub), *readers)
    max_diff = max(map(operator.abs, chain.from_iterable(diffs)))

print('mito',max_diff)

if max_diff > 1:
    exit(1)
