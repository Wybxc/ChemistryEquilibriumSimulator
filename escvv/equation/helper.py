
import math
import re
from enum import Enum
from typing import List, Tuple

from .. import utils
from ..state import ESCState
from . import ESCEquation


class _CharType(Enum):
    '''字符类型。使用 `_CharType(c)` 获取。
    '''
    WORD = 1
    DECIMAL = 2
    DOT = 3
    PLUS = 4
    EQUAL = 7
    SPACE = 8
    OTHER = 9

    @classmethod
    def _missing_(cls, char):
        # 从具体的字符获取字符类型
        if isinstance(char, str):
            if char.isalpha() or char in '-()':
                return _CharType.WORD
            elif char.isdecimal():
                return _CharType.DECIMAL
            elif char.isspace():
                return _CharType.SPACE
            else:
                return _CHAR_TYPES.get(char, _CharType.OTHER)
        return super()._missing_(char)


_CHAR_TYPES = {
    '.': _CharType.DOT,
    '+': _CharType.PLUS,
    '=': _CharType.EQUAL,
}


class _TokenizerState(Enum):
    '''词法分析器状态。
    '''
    NORMAL = 0
    WORD = 1
    DECIMAL = 3
    DECIMAL_ = 4
    DOT = 5
    PLUS = 6
    EQUAL = 7

    def next_state(self, char: str):
        '''获取当前状态的下一个状态。
        '''
        return _TOKENIZER_STATE_TABLE.get(self, {}) \
            .get(_CharType(char), _TokenizerState.NORMAL)

    def get(self):
        '''将辅助状态合并到正常词法分析状态中。
        '''
        if self == _TokenizerState.DECIMAL_:
            return _TokenizerState.DECIMAL
        return self


_TOKENIZER_STATE_TABLE = {
    _TokenizerState.NORMAL: {
        char_type: _TokenizerState[char_type.name]
        for char_type in _CharType
        if char_type not in (_CharType.OTHER, _CharType.SPACE)
    },
    _TokenizerState.WORD: {
        _CharType.WORD: _TokenizerState.WORD,
        _CharType.DECIMAL: _TokenizerState.WORD,
        _CharType.DOT: _TokenizerState.WORD,
    },
    _TokenizerState.DECIMAL: {
        _CharType.DECIMAL: _TokenizerState.DECIMAL,
        _CharType.DOT: _TokenizerState.DECIMAL_,
    },
    _TokenizerState.DECIMAL_: {
        _CharType.DECIMAL: _TokenizerState.DECIMAL_,
    },
    _TokenizerState.EQUAL: {
        _CharType.EQUAL: _TokenizerState.EQUAL,
    },
}


def _tokenize_origin(s: str):
    '''词法分析（不处理加号）
    '''
    state = _TokenizerState.NORMAL
    token = ''
    for c in s + '\n':  # 添加换行符以终止词法分析，保证最后一个 token 可以被读取
        next_state = state.next_state(c)
        if next_state != _TokenizerState.NORMAL:  # 状态跳转
            token += c
        elif state != _TokenizerState.NORMAL:  # 规约
            yield token.strip(), state.get()
            token = c
            next_state = _TokenizerState.NORMAL.next_state(c)
        state = next_state
    # 终止密码子（大嘘）
    yield '', _TokenizerState.NORMAL
    yield '', _TokenizerState.NORMAL
    yield '', _TokenizerState.NORMAL


@utils.endless(default=('', _TokenizerState.NORMAL))
def _tokenize(s: str):
    '''词法分析（处理加号）

    WORD 后的加号只允许出现在 WORD 和 DECIMAL 之前。
    '''
    tokenize_origin = _tokenize_origin(s)
    t0 = next(tokenize_origin)
    t1, t2 = next(tokenize_origin), next(tokenize_origin)
    while t0[1] != _TokenizerState.NORMAL:
        if t0[1] == _TokenizerState.WORD and t1[1] == _TokenizerState.PLUS \
                and t2[1] not in (_TokenizerState.WORD, _TokenizerState.DECIMAL):
            t0 = (t0[0] + '+', _TokenizerState.WORD)
            t1, t2 = t2, next(tokenize_origin)
        else:
            yield t0
            t0, t1, t2 = t1, t2, next(tokenize_origin)


class _Tokenizer:
    '''词法分析器。
    '''

    def __init__(self, s: str):
        self._tokenize = _tokenize(s)
        self._tokens = []
        self._offset = 0

    def next(self):
        '''下一个 token。
        '''
        if self._offset == 0:
            token = next(self._tokenize)
            self._tokens.append(token)
            return token
        else:
            offset = self._offset
            self._offset += 1
            return self._tokens[offset]

    def previous(self):
        '''回退一个 token。
        '''
        self._offset -= 1
        return self._tokens[self._offset]


# 化学方程式文法：
# Equation   ::= Substances EQUAL Substances 'K' '=' DECIMAL
#              | Substances EQUAL Substances
# Substances ::= Material '+' Substances
#              | Material
# Material   ::= DECIMAL Name
#              | Name
# Names      ::= WORD Name
#              | WORD


def _parse(s: str, default_sufficient=1., default_state=ESCState.GAS):
    '''语法分析
    '''
    tokens = _Tokenizer(s)

    # 以下用 * 表示当前解析位置
    def equation(tokens: _Tokenizer) -> ESCEquation:
        sub1 = substances(tokens)
        # Equation   ::= *Substances EQUAL Substances 'K' '=' DECIMAL
        #              | *Substances EQUAL Substances

        token, state = tokens.next()
        # Equation   ::= Substances *EQUAL Substances 'K' '=' DECIMAL
        #              | Substances *EQUAL Substances
        if state != _TokenizerState.EQUAL:
            raise SyntaxError(
                f'Invalid syntax, expect a equal, but got {token}')

        sub2 = substances(tokens)
        # Equation   ::= Substances EQUAL *Substances 'K' '=' DECIMAL
        #              | Substances EQUAL *Substances
        token, state = tokens.next()
        # Equation   ::= Substances EQUAL Substances *'K' '=' DECIMAL
        #              | Substances EQUAL Substances *
        if state == _TokenizerState.NORMAL:
            equilibrium_K = 0.
        else:
            if not(state == _TokenizerState.WORD and token == 'K'):
                raise SyntaxError(
                    f'Invalid syntax, expect `K`, but got {token}')

            token, state = tokens.next()
            # Equation   ::= Substances EQUAL Substances 'K' *'=' DECIMAL
            if not(state == _TokenizerState.EQUAL and token == '='):
                raise SyntaxError(
                    f'Invalid syntax, expect `=`, but got {token}')

            token, state = tokens.next()
            # Equation   ::= Substances EQUAL Substances 'K' '=' *DECIMAL
            if state != _TokenizerState.DECIMAL:
                raise SyntaxError(
                    f'Invalid syntax, expect a number, but got {token}')
            equilibrium_K = float(token)

        return ESCEquation({
            **{name: (sufficient, state) for name, sufficient, state in sub1},
            **{name: (-sufficient, state) for name, sufficient, state in sub2},
        }, equilibrium_K)

    def substances(tokens: _Tokenizer) -> List[Tuple[str, float, ESCState]]:
        materials = []
        while True:  # 尾递归展开
            materials.append(material(tokens))
            # Substances ::= *Material '+' Substances
            #              | *Material
            _, state = tokens.next()
            if state == _TokenizerState.PLUS:
                # Substances ::= Material *'+' Substances
                continue
                # Substances ::= Material '+' *Substances
            else:
                # Substances ::= Material *
                tokens.previous()
                break
        return materials

    def material(tokens: _Tokenizer) -> Tuple[str, float, ESCState]:
        token, state = tokens.next()
        # Material   ::= *DECIMAL Name
        #              | *Name
        if state == _TokenizerState.DECIMAL:
            sufficient = float(token)
            # Material   ::= DECIMAL *Name
        else:
            tokens.previous()
            sufficient = default_sufficient
        # Material   ::= DECIMAL *Name
        #              | *Name
        substance_name, substance_state = name(tokens)
        return substance_name, sufficient, substance_state

    def name(tokens: _Tokenizer) -> Tuple[str, ESCState]:
        substance_name_slice = []
        while True:  # 尾递归展开
            token, state = tokens.next()
            # Names      ::= *WORD Name
            #              | *WORD
            if state == _TokenizerState.WORD:
                substance_name_slice.append(token)
                continue
                # Names      ::= WORD *Name
            else:
                # Names      ::= WORD *
                tokens.previous()
                break

        if not substance_name_slice:
            raise SyntaxError(
                f'Invalid syntax, expect a name, but got {token}')

        # 解析物质状态（从名称的最后一项中）
        substance_state_match = re.match(
            r'^(.*)\((s|l|g|aq)\)$', substance_name_slice[-1])
        if substance_state_match:
            substance_name_slice[-1], substance_state_str = substance_state_match.groups()
            substance_state = ESCState(substance_state_str)

            # 处理最后一项只有物态表示的情况
            if not substance_name_slice[-1]:
                substance_name_slice.pop()
            if not substance_name_slice:
                raise SyntaxError(
                    f'Invalid syntax, expect a name, but got ({substance_state.value})')
        else:
            substance_state = default_state

        substance_name = '_'.join(substance_name_slice)
        return substance_name, substance_state

    return equation(tokens)


def parse_equation(text: str,
                   default_sufficient: float = 1.,
                   default_state: ESCState = ESCState.GAS,
                   ) -> ESCEquation:
    '''解析文字方程式。

    按照一般的化学方程式的书写语法，例如：

    `AgNO3(s) === Ag(s) + NO2(g) + 0.5O2(g)  K=0.3924`

    需要注意的地方：
    1. 物态用 `(s)` `(l)` `(g)` `(aq)` 表示，未指定的按照默认值处理。
    2. 平衡常数可以省略。若省略，默认值为0，可以在之后重新指定（比如通过`from_gibbs`计算）。
    3. 等号可以使用一个或多个，效果相同。
    4. 物质名称必须以字母开头。名称中可以带空格，解析时空格将被替换为下划线，多个空格视为一个。
    5. 物质名称仅用于区别不同物质，不必与化学式有关。

    参数：
    text: str 需解析的文字方程式
    default_sufficient: float 未指定物质系数时的默认系数，默认为1.0
    default_state: ESCState 未指定物质状态时的默认状态，默认为`ESCState.GAS`
    '''
    return _parse(text, default_sufficient, default_state)


def from_gibbs(enthalpy: float, entropy: float, temperature: float):
    '''通过熵变和焓变计算标准平衡常数。
    enthalpy: float 焓变，单位 kJ*mol^-1
    entropy: float 熵变，单位 kJ*(mol*K)^-1
    temperature: float 温度，单位 K
    '''
    gibbs = enthalpy - entropy * temperature
    return math.exp(-gibbs / temperature / 8.314472e-3)
