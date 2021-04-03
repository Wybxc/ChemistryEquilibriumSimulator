from typing import Dict, List, Tuple

import numpy as np

from ..state import ESCState


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
