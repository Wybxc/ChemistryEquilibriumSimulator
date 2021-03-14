# -*- coding: utf-8 -*-
from enum import Enum


class ESCState(Enum):
    '''物态，包括 SOLID（液体）, LIQUID（液体）, GAS（气体）, AQUEOUS（溶质）。
    '''
    SOLID = 's'
    LIQUID = 'l'
    GAS = 'g'
    AQUEOUS = 'aq'

    @property
    def index(self):
        return _STATE_INDEX[self]


_STATE_INDEX = {state: i for i, state in enumerate(ESCState)}
