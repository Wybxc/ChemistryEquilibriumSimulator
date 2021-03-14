# -*- coding: utf-8 -*-
from enum import Enum
from typing import Dict, List, Tuple

import numpy as np


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


class ESCEquation():
    substances: Dict[str, Tuple[float, ESCState]]
    equilibrium_K: float

    def __init__(self, substances, equilibrium_K):
        '''热化学方程式。
        substances: Dict[str, Tuple[float, ESCState]] 方程式中各物质的系数与状态
        equilibrium_K: float 标准平衡常数
        '''
        self.substances = substances
        self.equilibrium_K = equilibrium_K

    def __str__(self):
        forward, backward = [], []
        for name, (coefficient, state) in self.substances.items():
            if coefficient == int(coefficient):
                coefficient = int(coefficient)
            if coefficient == 1:
                forward.append(f'{name}({state.value})')
            elif coefficient == -1:
                backward.append(f'{name}({state.value})')
            elif coefficient > 0:
                forward.append(f'{coefficient}{name}({state.value})')
            elif coefficient < 0:
                backward.append(f'{-coefficient}{name}({state.value})')
        return f'{" + ".join(forward)} === {" + ".join(backward)}  K={self.equilibrium_K}'

    def forward(self, indexes: Dict[str, int], count: int):
        '''正反应的系数，以 ndarray 的形式返回。
        indexes: Dict[str, int] 由物质名称得到ID的表
        count: int 物质种类数
        '''
        result = np.zeros(count, dtype=np.float)
        for name, (coefficient, _) in self.substances.items():
            if coefficient > 0:
                result[indexes[name]] = coefficient
        return result

    def backward(self, indexes: Dict[str, int], count: int):
        '''逆反应的系数，以 ndarray 的形式返回。
        indexes: Dict[str, int] 由物质名称得到ID的表
        count: int 物质种类数
        '''
        result = np.zeros(count, dtype=np.float)
        for name, (coefficient, _) in self.substances.items():
            if coefficient < 0:
                result[indexes[name]] = -coefficient
        return result

    def coefficient(self, indexes: Dict[str, int], count: int):
        '''反应的系数，包括正反应和逆反应，以 ndarray 的形式返回。
        indexes: Dict[str, int] 由物质名称得到ID的表
        count: int 物质种类数
        '''
        return self.forward(indexes, count) - self.backward(indexes, count)


def _rate(substances, coefficient, rate_k=1.):
    '''计算反应速率。
    substances: 各物质的浓度。
    coefficient: 速率方程中反应物的系数。
    rate_k: 速率常数。
    '''
    return np.prod(np.float_power(substances, coefficient)) * rate_k


class ESCSystem():
    '''化学反应体系。
    '''
    equations: List[ESCEquation]

    PRECISION = 14

    def __init__(self,
                 substances: Dict[str, float] = {},
                 equations: List[ESCEquation] = []):
        '''化学反应体系。
        substances: Dict[str, float] 各物质的初始量，未在此处指定者视为0。
        equations: List[ESCEquation] 体系中存在的化学方程式。        
        '''
        self._substances_index: Dict[str, int] = {}  # 由物质名称得到ID的表
        self._substances_name: List[str] = []  # 各物质的名称
        self._substances_amount: List[float] = []  # 各物质的量
        self._substances_state: List[ESCState] = []  # 各物质的状态

        self.equations = []

        for name, concentration in substances.items():
            self[name] = concentration

        for equation in equations:
            self.add_equation(equation)

    def _set_state(self, name, state):
        self._substances_state[self._substances_index[name]] = state

    def add_equation(self, equation: ESCEquation):
        '''添加方程式。
        方程式中包含的新物质会自动加入体系中，初始量为0。
        equation: ESCEquation 方程式
        '''
        self.equations.append(equation)
        for name in equation.substances.keys():
            if name not in self._substances_name:
                self[name] = 0.
            _, state = equation.substances[name]
            self._set_state(name, state)

    def __getitem__(self, name):
        index = self._substances_index[name]
        return self._substances_amount[index]

    def __setitem__(self, name, concentration):
        if name in self._substances_name:
            index = self._substances_index[name]
            self._substances_amount[index] = concentration
        else:
            self._substances_name.append(name)
            self._substances_amount.append(concentration)
            self._substances_state.append(ESCState.GAS)
            self._substances_index[name] = len(self._substances_name) - 1

    def __len__(self):
        return len(self._substances_name)

    def __str__(self):
        equations = '  ' + '\n  '.join(map(str, self.equations))
        substances = '  ' + \
            '\n  '.join([f'{name}: {self._substances_amount[index]:.{self.PRECISION}f} mol' for name,
                         index in self._substances_index.items()])
        return f'''Equilibrium Simulator of Chemistry -- Reaction System
Equations:
{equations}
Substances: 
{substances}'''

    def _states(self):
        '''各物质的状态，以二维 ndarray 的形式返回。
        '''
        states = np.array([state.index for state in self._substances_state])
        return np.stack([states == state_index for state_index in range(len(ESCState))])

    def _concentration(self, amount, states):
        '''计算物质浓度，实现取决于具体的体系。
        '''
        amount = amount.copy()
        not_have_volume = states[ESCState.SOLID.index] | states[ESCState.AQUEOUS.index]
        np.place(amount, not_have_volume, 1.)
        return amount

    def run(self, delta=0.5, espilon=0.001, annealing=0.9):
        '''进行一次反应模拟，直到到达平衡态。
        delta: 单次反应进度的比率。
        espilon: 模拟精度，当 1-ε < Q/K < 1+ε 时认为反应平衡。
        annealing: 退火比例（每次模拟时 delta 递缩的系数）。
        '''
        coefficients = [(equation.forward(self._substances_index, len(self)),
                         equation.backward(self._substances_index, len(self)),
                         equation.equilibrium_K)
                        for equation in self.equations]

        equilibrium = False  # 是否达到平衡
        substances = np.array(self._substances_amount)
        states = self._states()
        not_have_volume = states[ESCState.SOLID.index] | states[ESCState.AQUEOUS.index]
        while not equilibrium:
            equilibrium = True
            for forward, backward, equilibrium_K in coefficients:
                concentration = self._concentration(substances, states)
                rate_forward = _rate(concentration, forward,
                                     equilibrium_K)  # 正反应速率
                rate_backward = _rate(concentration, backward, 1.)  # 逆反应速率
                rate = rate_forward - rate_backward
                if abs(rate) > espilon:  # 若未达到平衡 not(1-ε < Q/K < 1+ε)
                    step = (forward - backward) * rate
                    consume = step > 0
                    ratio = (substances[consume] / step[consume]).min()
                    if not ((substances - step * espilon) < 0).any():  # 在反应物即将耗尽时，尝试全部投入反应
                        ratio *= delta  # 否则按倍率进行
                    substances -= step * ratio
                    # 对于反应物或生成物无气体的反应
                    if not_have_volume[forward > 0].all() or not_have_volume[backward > 0].all():
                        if((substances - step * espilon / 10) < 0).any():  # 若反应物或生成物耗尽，不再计算平衡
                            continue
                    equilibrium = False
            if delta > espilon / 10:
                delta = delta * annealing
        substances = np.around(substances, self.PRECISION)
        substances[substances <= 0] = 0
        self._substances_amount = substances.tolist()
