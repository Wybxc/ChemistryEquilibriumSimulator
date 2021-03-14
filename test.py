from escvv import ESCState, ESCEquation
from escvv import ESCSystemGPCP

if __name__ == "__main__":
    system = ESCSystemGPCP(substances={
        'SO2': 1.,
        'O2': 1.,
    }, equations=[
        ESCEquation({
            'SO2': (2., ESCState.GAS),
            'O2': (1., ESCState.GAS),
            'SO3': (-2., ESCState.GAS),
        }, 0.9)
    ], volume = 0.5)
    system.run(delta=0.5, espilon=0.001, annealing=0.9)
    print(str(system))
