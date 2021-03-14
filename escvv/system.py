from typing import Dict, List

import numpy as np

from .equation import ESCEquation
from .state import ESCState


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
