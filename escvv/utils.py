from typing import Generator


def endless(default=None):
    '''无穷生成器。迭代完成后返回无限个 default。
    '''
    def wrapper(g: Generator):
        def wrapped(*args, **kwargs):
            for i in g(*args, **kwargs):
                yield i
            while True:
                yield default
        return wrapped
    return wrapper


def double(default=None):
    '''前瞻生成器，每次返回当前值和下一个值。
    '''
    def wrapper(g: Generator):
        g_endless = endless(default)(g)

        def wrapped(*args, **kwargs):
            iterator = g_endless(*args, **kwargs)
            t0 = next(iterator)
            while True:
                t1 = next(iterator)
                yield t0, t1
                t0 = t1
        return wrapped
    return wrapper
