from abc import ABC, abstractmethod

class Binner(ABC):
    name: str
    execs: str

    @abstractmethod
    def run_prep(self, asm_dir, bam_dir, samples, task_out_dir, **kwargs):
        ...

    @abstractmethod
    def run_binning(self, task_out_dir, options, **kwargs):
        ...

    @abstractmethod
    def construct_binning_table(self, task_name, sample, task_out_dir, ground_truth):
        ...

    @abstractmethod
    def prep_done(self, task_out_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def bin_done(self, task_out_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def post_done(self, task_out_dir, **kwargs) -> bool:
        ...

# Need to updatw refiner class to new api
class Refiner(ABC):
    name: str
    execs: str

    @abstractmethod
    def run_binning(self, task_out_dir, options, **kwargs):
        ...

    @abstractmethod
    def construct_binning_table(self, task_name, sample, task_out_dir, ground_truth):
        ...

    @abstractmethod
    def bin_done(self, task_out_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def post_done(self, task_out_dir, **kwargs) -> bool:
        ...

