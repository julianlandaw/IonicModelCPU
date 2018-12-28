#!/usr/bin/env python

from Cells_C import LR1CellIto
x = LR1CellIto()
print(x.getv(0))
x.stepdt(0.1,0,0)
print(x.getv(0))
