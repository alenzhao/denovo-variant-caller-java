#! /usr/bin/env python
#
# Copyright 2014 Google Inc. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
#
# in compliance with the License. You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
# or implied. See the License for the specific language governing permissions and limitations under
# the License.
#
import os
import sys
import numpy as np
import glob
import math
import matplotlib.pyplot as plt
from operator import itemgetter

from argparse import ArgumentParser 

def parse_file(f) : 
    d = {}
    with open(f) as fin :
        for line in fin :
            key, val = line.strip().split(":")
            key, val = key.strip(), val.strip()
            d.update([(key, val)])
    return d

def extract_arrays_threaded(listf) :
    x,y,sdev = [],[],[]
    for f in listf :
        parsedDict = parse_file(f)
        x.append(int(parsedDict.get('numThreads',1)))
        y.append(float(parsedDict.get('Mean/thread',parsedDict['Mean'])[:-2]))
        sdev.append(float(parsedDict['Sdev'][:-2])/ \
                int(parsedDict.get('numThreads',1)))
    
    x,y,sdev = zip(*[(x,y,z) for (x,y,z) in zip(x,y,sdev) if not math.isnan(y)])    
    x,y,sdev = zip(*sorted(zip(x,y,sdev), key = itemgetter(0)))
    return (x, y, sdev) 

def plot_threaded(x, y, sdev, title, figfile, offset = 10, xlab = "NumThreads",\
        ylab = "Time(ms)" ) : 
    fig, ax = plt.subplots()
    plt.errorbar(x, y, yerr = sdev)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.grid(True);
    ax.set_xlim(min(x)-offset, max(x)+offset)
    plt.title(title)
    plt.savefig(os.path.join(args.dir, figfile))
    plt.close()


if __name__ == '__main__'  :
    parser = ArgumentParser(description = "plot benchmarks")
    parser.add_argument("--dir", type = str, help = "Benchmark directory",\
            default = "~/.benchmarks")
    args = parser.parse_args()
    args.dir = os.path.expanduser(args.dir)


    # Work with varstore maxResults
    print("Files found in " +  args.dir)
    
    x,y,sdev = [],[],[]
    for f in glob.glob(args.dir+"/varstore*K") :
        parsedDict = parse_file(f)
        x.append(int(parsedDict['maxVariantResults']))
        y.append(float(parsedDict['Mean'][:-2]))
        sdev.append(float(parsedDict['Sdev'][:-2]))
        
    x,y,sdev = zip(*sorted(zip(x,y,sdev), key = itemgetter(0)))

    fig, ax = plt.subplots()
    plt.errorbar(x, y, yerr = sdev)
    ax.set_ylabel("Time(ms)")
    ax.set_xlabel("MaxResults")
    ax.grid(True);
    ax.set_xlim(min(x)-1000, max(x)+1000)
    plt.title("Varstore MaxResults")
    plt.savefig(os.path.join(args.dir, "varstoreMaxResults.png"))
    plt.close()


    # Work with varstore Threaded
    x,y,sdev = extract_arrays_threaded(glob.glob(args.dir+"/varstore10K*"))
    plot_threaded(x, y, sdev, "Varstore threaded","varstoreThreaded.png")

    # Work with ReadStore Threaded and Contig Also
    readstore_files = set(glob.glob(args.dir+"/readstore*"))
    readstore_files = set([f for f in readstore_files if 'png' not in f])
    readstore_contig_files = [f for f in readstore_files if 'Contig' in f]
    readstore_files = readstore_files.difference(readstore_contig_files)

    x,y,sdev = extract_arrays_threaded(readstore_files)
    x2,y2,sdev2 = extract_arrays_threaded(readstore_contig_files)
    
    fig, ax = plt.subplots()
    plt.errorbar(x, y, yerr = sdev, label = "random")
    plt.errorbar(x2, y2, yerr = sdev2, label = "contiguous")
    ax.set_ylabel("Time(ms)")
    ax.set_xlabel("NumThreads")
    ax.grid(True);
    ax.set_xlim(min(x)-5, max(x)+5)
    ax.legend()
    plt.title("Readstore Threaded")
    plt.savefig(os.path.join(args.dir, "readstoreThreaded.png"))
    plt.close()

