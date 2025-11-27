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

from yeast_mito_sim import (
    cell, 
    family_simulator, 
    print_family_table
)

f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1), 6)

print_family_table(f, g, "init_as_data.csv", prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1), 6)
    print_family_table(f, g, "init_as_data.csv", prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1), 6)
print_family_table(f, g, "init_as_data.csv", prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1), 6)
print_family_table(f, g, "init_as_data.csv", prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1), 6)
print_family_table(f, g, "init_as_data.csv", prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1), 6)
print_family_table(f, g, "init_as_data.csv", prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1), 6)
print_family_table(f, g, "init_as_data.csv", prefix="16_", add=True)

## init_as_data.csv moved to ../simulation_results

f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, nspl=2), 6)

print_family_table(f, g, "init_as_data_2.csv", prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, nspl=2), 6)
    print_family_table(f, g, "init_as_data_2.csv", prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, nspl=2), 6)
print_family_table(f, g, "init_as_data_2.csv", prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=2), 6)
print_family_table(f, g, "init_as_data_2.csv", prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=2), 6)
print_family_table(f, g, "init_as_data_2.csv", prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1, nspl=2), 6)
print_family_table(f, g, "init_as_data_2.csv", prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1, nspl=2), 6)
print_family_table(f, g, "init_as_data_2.csv", prefix="16_", add=True)

## init_as_data_2.csv moved to ../simulation_results

f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, nspl=4), 6)

print_family_table(f, g, "init_as_data_4.csv", prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, nspl=4), 6)
    print_family_table(f, g, "init_as_data_4.csv", prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, nspl=4), 6)
print_family_table(f, g, "init_as_data_4.csv", prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=4), 6)
print_family_table(f, g, "init_as_data_4.csv", prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=4), 6)
print_family_table(f, g, "init_as_data_4.csv", prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1, nspl=4), 6)
print_family_table(f, g, "init_as_data_4.csv", prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1, nspl=4), 6)
print_family_table(f, g, "init_as_data_4.csv", prefix="16_", add=True)

## init_as_data_4.csv moved to ../simulation_results

f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, nspl=8), 6)

print_family_table(f, g, "init_as_data_8.csv", prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, nspl=8), 6)
    print_family_table(f, g, "init_as_data_8.csv", prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, nspl=8), 6)
print_family_table(f, g, "init_as_data_8.csv", prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=8), 6)
print_family_table(f, g, "init_as_data_8.csv", prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=8), 6)
print_family_table(f, g, "init_as_data_8.csv", prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1, nspl=8), 6)
print_family_table(f, g, "init_as_data_8.csv", prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1, nspl=8), 6)
print_family_table(f, g, "init_as_data_8.csv", prefix="16_", add=True)

## init_as_data_8.csv moved to ../simulation_results

f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, nspl=10), 6)

print_family_table(f, g, "init_as_data_10.csv", prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, nspl=10), 6)
    print_family_table(f, g, "init_as_data_10.csv", prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, nspl=10), 6)
print_family_table(f, g, "init_as_data_10.csv", prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=10), 6)
print_family_table(f, g, "init_as_data_10.csv", prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=10), 6)
print_family_table(f, g, "init_as_data_10.csv", prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1, nspl=10), 6)
print_family_table(f, g, "init_as_data_10.csv", prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1, nspl=10), 6)
print_family_table(f, g, "init_as_data_10.csv", prefix="16_", add=True)

## init_as_data_10.csv moved to ../simulation_results


f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, ndau=5), 6)

print_family_table(f, g, "init_as_data_ndau5.csv", prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, ndau=5), 6)
    print_family_table(f, g, "init_as_data_ndau5.csv", prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, ndau=5), 6)
print_family_table(f, g, "init_as_data_ndau5.csv", prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, ndau=5), 6)
print_family_table(f, g, "init_as_data_ndau5.csv", prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, ndau=5), 6)
print_family_table(f, g, "init_as_data_ndau5.csv", prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1, ndau=5), 6)
print_family_table(f, g, "init_as_data_ndau5.csv", prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1, ndau=5), 6)
print_family_table(f, g, "init_as_data_ndau5.csv", prefix="16_", add=True)

## init_as_data_ndau5.csv moved to ../simulation_results

f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, nspl=10, ndau=5), 6)

print_family_table(f, g, "init_as_data_10_ndau5.csv", prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, nspl=10, ndau=5), 6)
    print_family_table(f, g, "init_as_data_10_ndau5.csv", prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, nspl=10, ndau=5), 6)
print_family_table(f, g, "init_as_data_10_ndau5.csv", prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=10, ndau=5), 6)
print_family_table(f, g, "init_as_data_10_ndau5.csv", prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=10, ndau=5), 6)
print_family_table(f, g, "init_as_data_10_ndau5.csv", prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1, nspl=10, ndau=5), 6)
print_family_table(f, g, "init_as_data_10_ndau5.csv", prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1, nspl=10, ndau=5), 6)
print_family_table(f, g, "init_as_data_10_ndau5.csv", prefix="16_", add=True)

## init_as_data_10_ndau5.csv moved to ../simulation_results

f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, ndau=8, simulrep=8, maxnaddexp=2), 6)

print_family_table(f, g, "init_as_data_882.csv", prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, ndau=8, simulrep=8, maxnaddexp=2), 6)
    print_family_table(f, g, "init_as_data_882.csv", prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, ndau=8, simulrep=8, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_882.csv", prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, ndau=8, simulrep=8, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_882.csv", prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, ndau=8, simulrep=8, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_882.csv", prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1, ndau=8, simulrep=8, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_882.csv", prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1, ndau=8, simulrep=8, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_882.csv", prefix="16_", add=True)

## init_as_data_10_ndau_882.csv moved to ../simulation_results

f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, ndau=6, simulrep=8, maxnaddexp=2), 6)

print_family_table(f, g, "init_as_data_682.csv", prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, ndau=6, simulrep=8, maxnaddexp=2), 6)
    print_family_table(f, g, "init_as_data_682.csv", prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, ndau=6, simulrep=8, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_682.csv", prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, ndau=6, simulrep=8, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_682.csv", prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, ndau=6, simulrep=8, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_682.csv", prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1, ndau=8, simulrep=6, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_682.csv", prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1, ndau=8, simulrep=6, maxnaddexp=2), 6)
print_family_table(f, g, "init_as_data_682.csv", prefix="16_", add=True)

## init_as_data_10_ndau_682.csv moved to ../simulation_results

for nspl in (2,4,6,8,10) :
    for ndau in (6, 8, 10) :
        f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        outfilename = "init_as_data_maxna2_nspl" + str(nspl) + "_ndau" + str(ndau) + ".csv"
        print_family_table(f, g, outfilename, prefix="0_", add=False)
        
        for i in range(11) :
            f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
            print_family_table(f, g, outfilename, prefix=str(i+1)+"_", add=True)
            
        f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="12_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="13_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="14_", add=True)
        
        f, g = family_simulator( cell((1,1,1,1), -1, nspl=nspl, ndau=ndau,  maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="15_", add=True)
        
        f, g = family_simulator( cell((0,0,0), -1, nspl=nspl, ndau=ndau,  maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="16_", add=True)

## init_as_data_....csv moved to ../simulation_results


for nspl in (2,4,6,8,10) :
    for ndau in (6, 8, 10) :
        f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        outfilename = "init_as_data_maxna2_nspl" + str(nspl) + "_ndau" + str(ndau) + "_S.csv"
        print_family_table(f, g, outfilename, prefix="0_", add=False)
        
        for i in range(11) :
            f, g = family_simulator( cell((1,0,1,0,1,0,1,0,1,0,1,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
            print_family_table(f, g, outfilename, prefix=str(i+1)+"_", add=True)
            
        f, g = family_simulator( cell((0,1,0,1,0,1,0,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="12_", add=True)
        
        f, g = family_simulator( cell((0,1,0,1,0,1,0,1,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="13_", add=True)
        
        f, g = family_simulator( cell((0,1,0,1,0,1,0,1,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="14_", add=True)
        
        f, g = family_simulator( cell((1,1,1,1), -1, nspl=nspl, ndau=ndau,  maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="15_", add=True)
        
        f, g = family_simulator( cell((0,0,0), -1, nspl=nspl, ndau=ndau,  maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="16_", add=True)

## init_as_data_....csv moved to ../simulation_results


for nspl in (2,4,6,8,10) :
    for ndau in (6, 8, 10) :
        f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=30, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        outfilename = "init_as_data_maxna2__sb_nspl" + str(nspl) + "_ndau" + str(ndau) + "_S.csv"
        print_family_table(f, g, outfilename, prefix="0_", add=False)
     
        for i in range(11) :
            f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
            print_family_table(f, g, outfilename, prefix=str(i+1)+"_", add=True)
            
        f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="12_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="13_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="14_", add=True)
        
        f, g = family_simulator( cell((1,1,1,1), -1, nspl=nspl, ndau=ndau,  maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="15_", add=True)
        
        f, g = family_simulator( cell((0,0,0), -1, nspl=nspl, ndau=ndau,  maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="16_", add=True)


## init_as_data_....csv moved to ../simulation_results



for nspl in (2,4,6,8,10) :
    for ndau in (6, 8, 10) :
        f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=25, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        outfilename = "init_as_data_maxna2_s25_nspl" + str(nspl) + "_ndau" + str(ndau) + "_S.csv"
        print_family_table(f, g, outfilename, prefix="0_", add=False)
     
        for i in range(11) :
            f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
            print_family_table(f, g, outfilename, prefix=str(i+1)+"_", add=True)
            
        f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="12_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="13_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1, nspl=nspl, ndau=ndau, maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="14_", add=True)
        
        f, g = family_simulator( cell((1,1,1,1), -1, nspl=nspl, ndau=ndau,  maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="15_", add=True)
        
        f, g = family_simulator( cell((0,0,0), -1, nspl=nspl, ndau=ndau,  maxnaddexp=2), 6)
        print_family_table(f, g, outfilename, prefix="16_", add=True)


## init_as_data_....csv moved to ../simulation_results


## simulation with best-fitting model:

f, g = family_simulator( cell((1,0,1,0,1,0,1,0,1,0,1,0), -1, nspl=4, ndau=6, maxnaddexp=2), 6)
plot_family(f,g,"")

## Simulation with best-fitting model for 16 Generations (24 h)

f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, nspl=4, ndau=10, maxnaddexp=2), 16)
outfilename = "init_as_data_maxna2_nspl" + str(4) + "_ndau" + str(10) + "_S_24.csv"
print_family_table(f, g, outfilename, prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,0,1,0,1,0,1,0,1,0,1,0), -1, nspl=4, ndau=10, maxnaddexp=2), 16)
    print_family_table(f, g, outfilename, prefix=str(i+1)+"_", add=True)
    
f, g = family_simulator( cell((0,1,0,1,0,1,0,0,0,0), -1, nspl=4, ndau=10, maxnaddexp=2), 16)
print_family_table(f, g, outfilename, prefix="12_", add=True)

f, g = family_simulator( cell((0,1,0,1,0,1,0,1,0,0), -1, nspl=4, ndau=10, maxnaddexp=2), 16)
print_family_table(f, g, outfilename, prefix="13_", add=True)

f, g = family_simulator( cell((0,1,0,1,0,1,0,1,0,0), -1, nspl=4, ndau=10, maxnaddexp=2), 16)
print_family_table(f, g, outfilename, prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1, nspl=4, ndau=10,  maxnaddexp=2), 16)
print_family_table(f, g, outfilename, prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1, nspl=4, ndau=10,  maxnaddexp=2), 16)
print_family_table(f, g, outfilename, prefix="16_", add=True)

## init_as_data_maxna2_nspl4_ndau10_S_24.csv  moved to ../simulation_results

###############################################################################
### Simulation with on average 33 nucleoids when splitting

for nspl in (1, 2,4,6,8,10,12,14,15, 50) :
    for ndau in (2, 6, 8, 9, 10, 11, 12, 14, 20) :
        f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
        outfilename = "init_as_data_maxna2_s27p6_nspl" + str(nspl) + "_ndau" + str(ndau) + "_S.csv"
        print_family_table(f, g, outfilename, prefix="0_", add=False)
     
        for i in range(11) :
            f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
            print_family_table(f, g, outfilename, prefix=str(i+1)+"_", add=True)
            
        f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
        print_family_table(f, g, outfilename, prefix="12_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
        print_family_table(f, g, outfilename, prefix="13_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
        print_family_table(f, g, outfilename, prefix="14_", add=True)
        
        f, g = family_simulator( cell((1,1,1,1), -1,  startbud=27, nspl=nspl, ndau=ndau,  maxnaddexp=6), 6)
        print_family_table(f, g, outfilename, prefix="15_", add=True)
        
        f, g = family_simulator( cell((0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau,  maxnaddexp=6), 6)
        print_family_table(f, g, outfilename, prefix="16_", add=True)


## init_as_data_....csv moved to ../simulation_results

f, g = family_simulator( cell((1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0), -1, startbud=27, nspl=8, ndau=11,  maxnaddexp=6), 6)
plot_family(f,g,"")

f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0), -1, startbud=27, nspl=8, ndau=11,  maxnaddexp=6), 6)
plot_family(f,g,"")

## simulation with best-fitting model for 24 h:

nspl=8
ndau=11
f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 16)
outfilename = "init_as_data_maxna2_s27p6_nspl" + str(nspl) + "_ndau" + str(ndau) + "_S_24.csv"
print_family_table(f, g, outfilename, prefix="0_", add=False)

for i in range(11) :
    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 16)
    print_family_table(f, g, outfilename, prefix=str(i+1)+"_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 16)
print_family_table(f, g, outfilename, prefix="12_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 16)
print_family_table(f, g, outfilename, prefix="13_", add=True)

f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 16)
print_family_table(f, g, outfilename, prefix="14_", add=True)

f, g = family_simulator( cell((1,1,1,1), -1,  startbud=27, nspl=nspl, ndau=ndau,  maxnaddexp=6), 16)
print_family_table(f, g, outfilename, prefix="15_", add=True)

f, g = family_simulator( cell((0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau,  maxnaddexp=6), 16)
print_family_table(f, g, outfilename, prefix="16_", add=True)



##################################################################################################
## simulations 33 fixed with nspl=0 or 33:

for nspl in (1, 33) :
    for ndau in (8, 10, 12, 14, 16) :
        f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=33, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
        outfilename = "init_as_data_maxna2_s33p0_nspl" + str(nspl) + "_ndau" + str(ndau) + ".csv"
        print_family_table(f, g, outfilename, prefix="0_", add=False)
        
        for i in range(11) :
            f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1,  startbud=33, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            print_family_table(f, g, outfilename, prefix=str(i+1)+"_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1,  startbud=33, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
        print_family_table(f, g, outfilename, prefix="12_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=33, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
        print_family_table(f, g, outfilename, prefix="13_", add=True)
        
        f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=33, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
        print_family_table(f, g, outfilename, prefix="14_", add=True)
        
        f, g = family_simulator( cell((1,1,1,1), -1,  startbud=33, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
        print_family_table(f, g, outfilename, prefix="15_", add=True)
        
        f, g = family_simulator( cell((0,0,0), -1,  startbud=33, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
        print_family_table(f, g, outfilename, prefix="16_", add=True)
        
## init_as_data_....csv moved to ../simulation_results



###############################################################################
### Simulation with on average 33 nucleoids when splitting; 10-fold number of cells

for nspl in range(1,26) :
    for ndau in range(1, 16) :
            f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
            outfilename = "init_as_data_maxna2_s27p6_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold_S.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((1,1,1,0,0,0,1,1,1,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,1,1,1,0,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,1,1,1,1,0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau, maxnaddexp=6), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=27, nspl=nspl, ndau=ndau,  maxnaddexp=6), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=27, nspl=nspl, ndau=ndau,  maxnaddexp=6), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## init_as_data_....csv moved to ../simulation_results/10fold/


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

## moved to ../simulation_results/S32p0/
        
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


## init_as_data_....csv moved to ../simulation_results/10fold_26/

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


## init_as_data_....csv moved to ../simulation_results/10fold_26_mixed/

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


## init_as_data_....csv moved to ../simulation_results/10fold_38/

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


## init_as_data_....csv moved to ../simulation_results/10fold_38_mixed/
###############################################################################
### Simulation with on 90 nucleoids when splitting; 10-fold number of cells, more mixed initial conditions

for nspl in range(1,26) :
    ## for ndau in range(1, 16) :
    for ndau in range(16, 46) :
            f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=90, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            outfilename = "init_as_data_s90_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold_m.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=90, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((1,0,1,0,1,0,1,0,1,0,1,0), -1,  startbud=90, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,0,0,1,0,0,1,0), -1,  startbud=90, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,1,0,1,0,0,1,0), -1,  startbud=90, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,1,0,1,0,0,1,0), -1,  startbud=90, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=90, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=90, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## init_as_data_....csv moved to ../simulation_results/10fold_90_mixed/ on kauai

###############################################################################
### Simulation with on 56 nucleoids when splitting; 10-fold number of cells, more mixed initial conditions

for nspl in range(1,26) :
    for ndau in range(1, 29) :
            f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=56, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            outfilename = "init_as_data_s56_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold_m.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1), -1, startbud=56, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((1,0,1,0,1,0,1,0,1,0,1,0), -1,  startbud=56, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,0,0,1,0,0,1,0), -1,  startbud=56, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,1,0,1,0,0,1,0), -1,  startbud=56, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,1,0,1,0,1,0,0,1,0), -1,  startbud=56, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=56, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=56, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## init_as_data_....csv moved to ../simulation_results/10fold_56_mixed/ on kauai


### INITIAL SORTED <BEGIN>

###############################################################################
### Simulation with on 26 nucleoids when splitting; 10-fold number of cells, more mixed initial conditions

for nspl in range(1,26) :
    for ndau in range(1, 16) :
            f, g = family_simulator( cell((0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), -1, startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            outfilename = "init_as_data_s26_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold_m.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), -1, startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1), -1,  startbud=26, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=26, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=26, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## init_as_data_....csv moved to ../simulation_results/10fold_26_initsort/

###############################################################################
### Simulation with on 38 nucleoids when splitting; 10-fold number of cells

for nspl in range(1,26) :
    for ndau in range(1, 16) :
            f, g = family_simulator( cell((0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), -1, startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            outfilename = "init_as_data_s38_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), -1, startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1), -1,  startbud=38, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=38, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=38, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## init_as_data_....csv moved to ../simulation_results/10fold_38_initsort/

###############################################################################
### Simulation with on 32 nucleoids when splitting; 10-fold number of cells

for nspl in range(1,26) :
    for ndau in range(1, 16) :
            f, g = family_simulator( cell((0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), -1, startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
            outfilename = "init_as_data_s32_nspl" + str(nspl) + "_ndau" + str(ndau) + "_10fold.csv"
            print_family_table(f, g, outfilename, prefix="0_0_", add=False)
            for repe in range(9) :
                f, g = family_simulator( cell((0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), -1, startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="0_"+str(repe+1)+"_", add=True)

            for repe in range(10) :
                for i in range(11) :
                    f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                    print_family_table(f, g, outfilename, prefix=str(i+1)+"_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="12_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="13_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1), -1,  startbud=32, nspl=nspl, ndau=ndau, maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="14_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((1,1,1,1), -1,  startbud=32, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="15_"+str(repe+1)+"_", add=True)
                
                f, g = family_simulator( cell((0,0,0), -1,  startbud=32, nspl=nspl, ndau=ndau,  maxnaddexp=0), 6)
                print_family_table(f, g, outfilename, prefix="16_"+str(repe+1)+"_", add=True)


## init_as_data_....csv moved to ../simulation_results/10fold_32_initsort/


### INITIAL SORTED <END>
