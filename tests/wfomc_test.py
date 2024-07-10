import pytest
import glob
from pathlib import Path

from sampling_fo2.context.wfomc_context import WFOMCContext
from sampling_fo2.parser import parse_input
from sampling_fo2.wfomc import Algo, wfomc

current_path = Path(__file__).parent.absolute()
model_files = (current_path.parent / 'models').glob('*')
algos = [Algo.FASTERv2]


def wfomc_proxy(model_file, algo):
    problem = parse_input(model_file)
    return wfomc(WFOMCContext(problem), algo)


@pytest.mark.parametrize(
    'model_file, algo',
    [(str(model_file), algo) for model_file in model_files for algo in algos]
)
def test_model(model_file, algo):
    wfomc_proxy(model_file, algo)
