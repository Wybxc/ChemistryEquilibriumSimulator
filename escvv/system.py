from typing import Dict, List, Tuple
from .core import ESCState, ESCEquation, ESCSystem

import numpy as np


class ESCSystemGPCV(ESCSystem):
    '''气相定容反应体系（Gas Phase Constant Volume）。
    支持气体、固体、纯液体（不考虑多相混合）。
    '''

    def __init__(self,
                 volume: float = 1.,
                 substances: Dict[str, float] = {},
                 equations: List[ESCEquation] = []):
        '''气相定容反应体系（Gas Phase Constant Volume）。
        volumn: float 体系的体积，单位为 L。
        substances: Dict[str, float] 各物质的初始量，未在此处指定者视为0。
        equations: List[ESCEquation] 体系中存在的化学方程式。        
        '''
        super().__init__(substances=substances, equations=equations)
        self.volume = volume

    def _set_state(self, name, state):
        if state == ESCState.AQUEOUS:
            raise ValueError("GPCV System doesn't support aq!")
        super()._set_state(name, state)

    def __str__(self):
        return super().__str__() + f'''
Volume = {self.volume} L'''

    def _concentration(self, amount, not_have_volume):
        amount = amount.copy()
        np.place(amount, not_have_volume, self.volume)
        return amount / self.volume

    # TODO: 拆分 ESCSystem.run


class ESCSystemGPCP(ESCSystem):
    '''气相定压反应体系（Gas Phase Constant Pressure）。
    '''

    def __init__(self,
                 volume: float = 1.,
                 substances: Dict[str, float] = {},
                 equations: List[ESCEquation] = []):
        '''气相定压反应体系（Gas Phase Constant Pressure）。
        volumn: float 体系的初始体积，单位为 L。
        substances: Dict[str, float] 各物质的初始量，未在此处指定者视为0。
        equations: List[ESCEquation] 体系中存在的化学方程式。        
        '''
        super().__init__(substances=substances, equations=equations)
        self._clapeyron_contant = 0
        self.volume = volume

    def _set_state(self, name, state):
        if state == ESCState.AQUEOUS:
            raise ValueError("GPCP System doesn't support aq!")
        super()._set_state(name, state)

    def __str__(self):
        return super().__str__() + f'''
Volume = {self.volume} L'''

    def _gas_amount(self, amount, states):
        return amount[states[ESCState.GAS.index]].sum()

    def _concentration(self, amount, states):
        # 每次取浓度时，重新计算体系的容积
        self._volume = self._gas_amount(amount, states) * self._clapeyron_contant
        amount = amount.copy()
        not_have_volume = states[ESCState.SOLID.index] | states[ESCState.AQUEOUS.index]
        np.place(amount, not_have_volume, self.volume)
        return amount / self._volume

    @property
    def volume(self):
        return self._volume

    @volume.setter
    def volume(self, volume):
        self._volume = volume
        amount = np.array(self._substances_amount)
        is_gas = self._states()[ESCState.GAS.index]
        self._clapeyron_contant = volume / self._gas_amount(amount, is_gas)
