# 
# This file is part of the yeast_mito_segregation (https://github.com/statgenlmu/yeast_mito_segregation).
# Copyright (c) 2024 Dirk Metzler https://evol.bio.lmu.de/_statgen/
# 
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
## matplotlib.use('gtk3cairo')
import matplotlib.pyplot as plt
import copy

rng = np.random.default_rng(seed=42)

class cell:
    """yeast cell with mitochondrium (or mitochondrial network) 

    Attributes:
    -----------
    mito:     A list of 0s and 1s representing the alleles of the nucleoids
    mother:   The id of the mother of the cell
    startbud:     Number of nucleoids until the cell produces a daughter; default: 20; but see maxnaddexp
    nspl:     Number of fragments into which the mitochondrium is split before reshuffling reproduction; default: 6
    ndau:     Number of nucleoids that are passed on to the daughter; default: 10
    simulrep: Number of nucleoids that are reproduced before the next generation can reproduce, default 30
    maxnaddexp: geometric random variable on 0,1,2,.. with expectation value maxnaddexp is added to startbud
    id:       Cell id by which this cell is identified as a mother
    
    Methods:
    --------
    nuclrepro: choose a random nucleoid that is reproduced and puts its offspring/copy next to it
    grow:      let nucleoids multiply until there are maxn of them
    splitmito: split the mitochondrium into fragments 
    have_daughter: the mito first grows and then splits into nspl pieces half of which are passed on to daughter
    """

    def __init__(self, nuc, mother, startbud=20, nspl=6, ndau=10, simulrep=30, maxnaddexp=4):
        """Constructor to initialize the mito with a list of nucleoid alleles.

        Parameters
        ----------
        nuc :       tuple or list (perhaps of 0s and 1s) defining the nucleoid alleles
        mother:     The id of the mother of the cell
        startbud:   Number of nucleotides when the cell can produce a daughter; default:
        nspl:       Number of fragments into which the mitochondrium is split before reshuffling reproduction; default: 6
        ndau:       Number of nucloids that are passed on to the daughter; default: 10
        simulrep:   Number of nucleoids that are reproduced before the next generation can reproduce, default 30
        maxnaddexp: geometric random variable on 0,1,2,.. with expectation value maxnaddexp is added to startbud
        id:         Cell id by which this cell is identified as a mother
        """
        self.mito = list(nuc)
        self.mother = mother
        self.startbud=startbud
        self.nspl = nspl
        self.ndau = ndau
        self.simulrep = simulrep
        self.maxnaddexp = maxnaddexp
        self.maxn = self.startbud + rng.geometric(1/(self.maxnaddexp+1))-1
        self.id = 0

    def set_id(self, id):
        self.id = id
        
    def nuclrepro(self):
        """choose up to self.simulrep  random nucleoid that is reproduced an puts its offspring next to it."""
        n = min(self.simulrep, len(self.mito), self.maxn-len(self.mito))
        k = rng.choice(len(self.mito), size=n, replace=False)
        for i in sorted(k, reverse=True) :
            self.mito.insert(i + rng.integers(0,2), self.mito[i])

    def grow(self):
        """ let nucleoids multiply until there are maxn of them. """
        while len(self.mito) < self.maxn:
            self.nuclrepro()

    def splitmito(self):
        """ split the mitochondrium into fragments 

        mito is split into nspl parts that are at least one nucleoid long.
        
        Returns
        -------
        list of the fragments.
        """
        if self.nspl==1 :
            return [self.mito]
        
        min_nspl = min(self.nspl, len(self.mito)) 
        split_right = len(self.mito) - min_nspl  
        to_permute =  np.repeat((0,1), (len(self.mito) - min_nspl, min_nspl))
        split = rng.permutation(to_permute)
        a = 0
        b = 0
        i = 0
        mitolist = []
        while i < len(split) :
            b += 1
            if split[i] or i+1 == len(split) :
                mitolist.append(self.mito[a:b])
                a = b
            i += 1
        
        return mitolist

    def have_daughter(self):
        """ the cell first grows and then splits into nspl pieces half of which are passed on to daughter

        half of them (or half of nspl-1 of them) are randomly assembled and passed on to the daughter.
        The others are randomly assembled and kept.

        Returns
        -------
        the daughter cell
        """
        self.grow()
        spm = self.splitmito()
        perm = rng.permutation(range(len(spm)))
        ml = [spm[i] for i in perm]
        self.mito = list(np.concatenate(ml, axis=0))
        mito_before_division = self.mito.copy()
        # # Here it should probably be random if the daughter gets the first 
        # # half or the second half of the mito
        # inherit_first_half = rng.choice((0, 1))
        # if inherit_first_half:
        #     daumito = self.mito[0:self.ndau]
        #     self.mito = self.mito[self.ndau:self.maxn]
        # else:
        #     daumito = self.mito[-self.ndau:]
        #     self.mito = self.mito[0:self.maxn-self.ndau]
        daumito = self.mito[0:self.ndau]
        self.mito = self.mito[self.ndau:self.maxn]
        self.maxn = self.startbud + rng.geometric(1/(self.maxnaddexp+1))-1
        new_cell = cell(
            daumito, 
            self.id, 
            startbud=self.startbud, 
            nspl=self.nspl, 
            ndau=self.ndau, 
            simulrep=self.simulrep, 
            maxnaddexp=self.maxnaddexp
        )
        return new_cell

def determine_if_should_have_daughter(current_cell, p_have_daughter=1):
    """_summary_

    Parameters
    ----------
    current_cell : cell
        The current cell being iterated
    p_have_daughter : int, optional
        Probability to have a daughter. For WT is 1, meaning that WT always 
        divides. If less than 1, cell is growing slower than WT. 
        Default is 1   

    Returns
    -------
    int
        0 for no daughter, 1 for yes daughter
    """    
    if p_have_daughter >= 1:
        return 1
    
    p_choice = (1-p_have_daughter, p_have_daughter)
    does_have_daughter = rng.choice((0, 1), p=p_choice)
    return does_have_daughter

def calculate_p_have_daughter(current_cell, other_strain_growth_rate):
    """Calculate probability of current cell to have a daughter based on 
    other_strain_growth_rate

    Parameters
    ----------
    current_cell : cell
        The current cell being iterated
    other_strain_growth_rate : float
        Growth rate relative to WT

    Returns
    -------
    float
        Probability to have a daughter
    
    Notes
    -----
    The probability to have a daughter is calculated with linear interpolation 
    between (0, 1) for the WT and (1, other_strain_growth_rate) for the 
    other strain. The probability to have a daughter is then calculated as the 
    fraction of dividing cells compared to WT. 
    The status of the current cell is determined as the mean 
    of the current mito distribution of 0s and 1s
    
    """    
    mutant_ratio = np.mean(current_cell.mito)
    gr_mut = other_strain_growth_rate
    x = mutant_ratio
    q = 1.0
    m = gr_mut-1
    
    # Linear interpolation between (0, 1) for WT and (1, gr_mut) for mutant
    gr_current = m*x + q
    fract_div_cells_mut = np.power(2, gr_current) - 1
    
    p_have_daughter = fract_div_cells_mut
    
    return p_have_daughter

def generation_simulator(
        cell_family, other_strain_growth_rate=1.0
    ):
    """simulate reproduction of one generation of yeast

    Parameters
    ----------
    cell_family : Tuple
        the tuple of cells of the parental generation
    other_strain_growth_rate : float, optional
        Growth rate ratio of the other strain compared to WT, by default 1.0
    
    Returns
    -------
    the tuple of cells of the next generation
    """
    
    cfl = list(cell_family)
    n = len(cfl)
    for i in range(n):
        current_cell = cfl[i]
        p_have_daughter = calculate_p_have_daughter(
            current_cell, other_strain_growth_rate
        )
        does_have_daughter = determine_if_should_have_daughter(
            current_cell, p_have_daughter=p_have_daughter
        )
        if not does_have_daughter:
            continue
        
        k = len(cfl)
        d = current_cell.have_daughter()
        current_cell.grow()
        cfl.append(d)
        cfl[k].set_id(k)
        cfl[k].grow()
        
    return tuple(cfl)

def family_simulator(cell, ngen, start_cell_id=0, other_strain_growth_rate=1.0) :
    """ simulate reproduction of one generation of yeast

    Paramter
    --------
    cell: the ancestral cell of the family
    ngen: number of generations to be simulated

    Returns
    -------
    tuple of the tuple of cells for all simulated generations and the tuple of generation numbers
    """
    cell.grow()
    p = (cell,)
    g = (start_cell_id, )
    fam = copy.deepcopy(p)
    for i in range(ngen):
        f = generation_simulator(
            p, other_strain_growth_rate=other_strain_growth_rate
        )
        g = g + tuple(np.repeat(i+1, len(f)))
        fam = fam + copy.deepcopy(f)
        p = f
    return((fam, g))

def plot_family(fam, g, title) :
    """ Print tree of frequencies of allele 1 for a simulated family of yeast cells.
    Blue lines refer to allele frequency changes in the mother cell.

    Parameter
    ---------
    fam:     the simulated family
    g:       tuple of generation numbers for each entry of fam
    title:   the title to be written over the plot
    """
    plt.title(title)
    plt.xlabel('Generation')
    plt.ylabel('Frequency of allele 1')
    d = {}  ### dictionary assigning to (cell_id, generation) the fraction of 1 alleles
    gb = {}  ### dictionary assigning to cell_id the generation in which the cell was born
    for i in range(len(fam)) :
        d[(fam[i].id, g[i])] = sum(fam[i].mito)/len(fam[i].mito)
        if not fam[i].id in gb :
            gb[fam[i].id] = g[i]
        if g[i] > 0 :
            if gb[fam[i].id]==g[i] :
                plt.plot([g[i]-1, g[i]], [d[(fam[i].mother, g[i]-1)], d[(fam[i].id, g[i])]], 'k-')
            else :
                plt.plot([g[i]-1, g[i]], [d[(fam[i].id, g[i]-1)], d[(fam[i].id, g[i])]], 'b-')
    return(plt.show())

def freq_spect(fam, g) :
    """ Calculate the number of cells of each allele frequency bin 0-0.05, 0.05-0.15,..., 0.95-1 for each time point.

    Parameter
    ---------
    fam:     the simulated family
    g:       tuple of generation numbers for each entry of fam

    Returns
    -------
    an array in which a[i,j] is the number of cells with j 1-alleles at time i
    """
    n = g[-1]
    ## m = len(fam[0].mito)
    res = np.zeros( (n, 11) )
    for i in range(len(fam)) :
        res[g[i]-1, int(round(np.mean(fam[i].mito)*10))] += 1
    return(res)

def print_family_table(fam, g, filename, prefix="", add = True, doubling_time=1.5) :
    """ Print comma-separated essentials of simulation results to a file

    Parameter
    ---------
    fam:      the simulated family
    g:        tuple of generation numbers for each entry of fam
    filename: name of the file to print to
    prefix:   prefix to be used to make cell ids unique when
              results of several simulations are combined in a file
    add:      If false the table header is output to the file and an existing file
              would be overwritten.
    """
    birthgen = {}
    if add:
        outf = open(filename, 'a')
    else:
        outf = open(filename, 'w')
        outf.write("time,Cell_ID,Mother,h,time_start,time_alive\n")
        
    for i in range(len(fam)) :
        if fam[i].id not in birthgen :
            birthgen[fam[i].id] = g[i]
        outf.write(
            str(g[i] * doubling_time) # time
            + ',' 
            + str(prefix) + str(fam[i].id) # Cell_ID
            + ',' 
            +  str(prefix) + str(fam[i].mother) # Mother
            + ',' 
            + str(np.mean(fam[i].mito)) # h
            + ',' 
            + str(birthgen[fam[i].id]*doubling_time) # time_start
            + ',' 
            + str((g[i] - birthgen[fam[i].id])*doubling_time) + '\n' # time_alive
        )
    outf.close()

def get_family_table(fam, g, prefix="", doubling_time=1.5, index_col='') :
    """ Print comma-separated essentials of simulation results to a file

    Parameter
    ---------
    fam:      the simulated family
    g:        tuple of generation numbers for each entry of fam
    prefix:   prefix to be used to make cell ids unique when
              results of several simulations are combined in a file
    """
    birthgen = {}
    df_data = defaultdict(list)
    
    for i in range(len(fam)) :
        if fam[i].id not in birthgen :
            birthgen[fam[i].id] = g[i]
        
        df_data['time'].append(g[i] * doubling_time)
        df_data['Cell_ID'].append(f'{prefix}{fam[i].id}')
        df_data['Mother'].append(f'{prefix}{fam[i].mother}')
        df_data['h'].append(np.mean(fam[i].mito))
        df_data['time_start'].append(birthgen[fam[i].id]*doubling_time)
        df_data['time_alive'].append((g[i] - birthgen[fam[i].id])*doubling_time)
        
    df = pd.DataFrame(df_data)
    if not index_col:
        return df
    
    df = df.set_index(index_col)
    return df

def plot_freq_spect(frsp) :
    """Plots generation-wise frequency spectra

    Parameter
    ---------
    frsp:     an np.array in which frsp[i,j] is the number of cells with j 1-alleles at time i
    """
    n, m = frsp.shape
    plt.title("Blue: initial, red: generation "+str(n))
    plt.xlabel('Freqeuncy of allele 1')
    plt.ylabel('Relative frequency among cells of allele 1')
    for i in range(n) : 
        plt.plot(np.repeat(range(0,m+1),2)[1:-1],
                 np.repeat(frsp[i,]/sum(frsp[i,]), 2),
                 color = (i/(n-1), 0 ,((n-i-1)/(n-1))) )
    return(plt.show())
