from abc import ABC, abstractmethod
import gzip
from os.path import exists, join, dirname, normpath
from os import rename
import re
from ..utils import parse_fasta, rm_dir, run


class Assembler(ABC):
    name: str
    execs: str

    @abstractmethod
    def run_main(self, r1, r2, asm_dir, threads, options, **kwargs):
        """
        Runs the assembly algorithm according to its command-line interface.
        """
        ...

    @abstractmethod
    def clean_up(self, asm_dir, **kwargs):
        """
        Removes unneeded files, leaving only the assembly fasta and renames the contigs.
        """
        ...

    @abstractmethod
    def reduce_header(self, header):
        """
        Take the header string from an assembly object and removes any metadata, leaving
        only the contig name with no whitespace."""
        ...

    def main_done(self, asm_dir, **kwargs):
        return exists(join(asm_dir, f'contigs.fasta'))

    def prep_done(self, **kwargs):
        return True
    
    def reduce_fasta_header_to_contig_name(self, file, gz=False):
        fs = open(file) if not gz else gzip.open(file)
        data = list(parse_fasta(fs, gz))
        fs.close()
        with open(file,'w') as f:
            for e in data:
                e.hdr = self.reduce_header(e.hdr)
                e.write(f)


class MEGAHITAssembler(Assembler):
    name = 'MEGAHIT'
    exec = 'megahit'

    def run_main(self, r1, r2, asm_dir, threads, options, **kwargs):
        rm_dir(asm_dir, remake=True)
        cmd = f'{self.exec} -t {threads} {options} -1 {r1} -2 {r2} -f -o {asm_dir}/asm'
        run(cmd)
        self.clean_up(asm_dir)
    
    def reduce_header(self, header):
        return header.split()[0]

    def clean_up(self, asm_dir, **kwargs):
        old_fa_name = join(f'{asm_dir}/asm', 'final.contigs.fa')
        new_fa_name= join(asm_dir, f'contigs.fasta')
        if exists(old_fa_name):
            rename(old_fa_name, new_fa_name)
            rm_dir(f'{asm_dir}/asm')
            self.reduce_fasta_header_to_contig_name(new_fa_name)

        
class metaSPAdesAssembler(Assembler):
    name = 'metaSPAdes'
    exec = 'metaspades.py'

    def run_main(self, r1, r2, asm_dir, threads, options, **kwargs):
        rm_dir(asm_dir, remake=True)
        cmd = f'{self.exec} -t {threads} {options} -1 {r1} -2 {r2} -o {asm_dir}/asm'
        run(cmd)
        self.clean_up(asm_dir)

    def reduce_header(self, header):
        return re.findall('(NODE_[0-9]+)_', header)[0]

    def clean_up(self, asm_dir, **kwargs):
        old_fa_name = join(f'{asm_dir}/asm', 'contigs.fasta')
        new_fa_name= join(asm_dir, f'contigs.fasta')
        if exists(old_fa_name):
            rename(old_fa_name, new_fa_name)
            rm_dir(f'{asm_dir}/asm')
            self.reduce_fasta_header_to_contig_name(new_fa_name)