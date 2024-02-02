#! /usr/bin/env python3.6
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

from yeast_mito_sim import *


##################################################################################################
## simulations 32 fixed:

for nspl in range(1, 26) :
    for ndau in range(1, 17) :
        f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
        outfilename = "init_as_data_s32p0_nspl" + str(nspl) + "_ndau" + str(ndau) + ".csv"
        print_family_table(f, g, outfilename, prefix="0_0_", add=False)
        
        for repe in range(9) :
            f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)
        
        for repe in range(10) :
            for i in range(11) :
                f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
            
            f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
            
            f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
            
            f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
            
            f, g = family_simulator( cell((1,1,1,1), -1,  startbud=32, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
            
            f, g = family_simulator( cell((0,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)

## move results to appropriate folder
        
##################################################################################################
## simulations 32 mixed, initial state more mixed:

for nspl in range(1, 26) :
    for ndau in range(1, 17) :
        f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
        outfilename = "init_as_data_mixed_s32p0_nspl" + str(nspl) + "_ndau" + str(ndau) + ".csv"
        print_family_table(f, g, outfilename, prefix="0_0_", add=False)
        
        for repe in range(9) :
            f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)
        
        for repe in range(10) :
            for i in range(11) :
                f, g = family_simulator( cell((1,0,1,0,1,0,1,0,1,0,1,0), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
            
            f, g = family_simulator( cell((0,1,0,1,0,1,0,1,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
            f, g = family_simulator( cell((0,1,0,1,0,1,0,1,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
            f, g = family_simulator( cell((0,1,0,1,0,1,0,1,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
            
            f, g = family_simulator( cell((1,1,1,1), -1,  startbud=32, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
            f, g = family_simulator( cell((0,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


###############################################################################
### Simulation with on 26 nucleoids when splitting; 10-fold number of cells

for nspl in range(1,26) :
    for ndau in range(1, 16) :
            f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            outfilename = "init_as_data_s26_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=26, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=26, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## move results to appropriate folder

###############################################################################
### Simulation with on 26 nucleoids when splitting; 10-fold number of cells, more mixed initial conditions

for nspl in range(1,26) :
    for ndau in range(1, 16) :
            f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            outfilename = "init_as_data_s26_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold_m.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((1,0,1,0,1,0,1,0,1,0,1,0), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,0,0,1,0,0,1,0), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,1,0,1,0,0,1,0), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,1,0,1,0,0,1,0), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=26, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=26, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## move results to appropriate folder

###############################################################################
### Simulation with on 38 nucleoids when splitting; 10-fold number of cells

for nspl in range(1,26) :
    for ndau in range(1, 16) :
            f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            outfilename = "init_as_data_s38_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=38, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=38, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## move results to appropriate folder

###############################################################################
### Simulation with on 38 nucleoids when splitting; 10-fold number of cells, more mixed initial conditions

for nspl in range(1,26) :
    for ndau in range(1, 16) :
            f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            outfilename = "init_as_data_s38_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold_m.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((1,0,1,0,1,0,1,0,1,0,1,0), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,0,0,1,0,0,1,0), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,1,0,1,0,0,1,0), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,1,0,1,0,0,1,0), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=38, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=38, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## move results to appropriate folder
