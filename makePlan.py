#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ingress Maxfield - makePlan.py

GNU Public License
http://www.gnu.org/licenses/
Copyright(C) 2016 by
Jonathan Baker; babamots@gmail.com
Trey Wenger; tvwenger@gmail.com
Christopher Goyet; ingress-maxfield@matheuxfous.xyz

Maxfield is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Maxfield is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Maxfield.  If not, see <http://www.gnu.org/licenses/>.

Original version by jpeterbaker
22 July 2014 - tvw updates csv file format
15 August 2014 - tvw updates with google API, adds -s,
                 switchted to ; delimited file
29 Sept 2014 - tvw V2.0 major update to new version
21 April 2015 - tvw V2.1 force data read to be string
26 Feb 2016 - tvw v3.0
              merged some new stuff from jpeterbaker's new version
02 Mar 2016 - tvw v3.1
              added option to skip link-by-link plots
              added timeout
11 Mai 2018 - v3.1.2 Goyet Christopher 
              rewrite source code with a lot of functions
              removed useless requests to google
              changed pickle output file in save_portals.pkle
              added get_output_file pickle file
              added TOPCOM option
              added noreduce path option
              changed optimal option, which can be really slow
              added optimize options : MU, AP, Fan, Length (of the path) or Keys (TODO: keys not really working)
            
"""

import sys
import subprocess
import os
import argparse
import networkx as nx
import numpy as np
import pandas as pd
from lib import PlanPrinterMap,geometry,agentOrder
from lib import Maxfield
from lib import Triangle as tr
import pickle
import copy
import time
from pebble import ProcessPool # to handle timeout
from concurrent.futures import TimeoutError # to handle timeout
import json

import matplotlib.pyplot as plt

# CONSTANTS
# version number
_V_ = '3.1.0'
# max portals allowed
_MAX_PORTALS_ = 1000
# number of attempts to try to get best plan
#_NUM_ATTEMPTS = 100    
_NUM_ATTEMPTS = 3
GREEN = '#3BF256' # Actual faction text colors in the app
BLUE  = '#2ABBFF'

# initiation functions 
def init_log_file(args):
    if args.log is not None:
        sys.stdout = open(args.log,'w',0)

def get_faction_color(args):
    if args.res:
        return BLUE
    else:
        return GREEN

def get_output_dir(args):
    output_directory = args.output_dir
    if not output_directory:
        output_directory = os.getcwd()
    # add ending separator
    if output_directory[-1] != os.sep:
        output_directory += os.sep
    # create directory if doesn't exist
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)
    return output_directory

def get_output_file(args):
    output_file = args.output_file
    if output_file[-4:] != '.pkl':
        output_file += ".pkl"
    return  output_file
    
def get_save_portals_file(args):
    save_portals_file = args.save_portals
    return save_portals_file
    
def get_input_file(args):
    input_file = args.input_file
    return input_file
    
def get_nagents(args):
    nagents = args.num_agents
    if nagents <= 0:
        print "Number of agents should be greater than zero"
        raise ValueError("Number of agents should be greater than zero")
    return nagents

def init_files(args):
    init_log_file(args)
    output_directory = get_output_dir(args)
    output_file = get_output_file(args)
    save_portals = get_save_portals_file(args)
    input_file = get_input_file(args)
    return output_directory, output_file, save_portals, input_file

def init_vars(args):
    color = get_faction_color(args)        
    nagents = get_nagents(args)
    optimize = args.optimize.upper()
    suboptimal = (not args.optimal)
    return color, args.google, args.api_key, nagents, args.attempts, suboptimal,\
    args.topcom , args.quiet, args.skipplot, args.noreduce, args.input_plan, optimize
    
def check_portal_formatting_input(portal):
    if len(portal) < 4:
        print "Error! Portal ",portal[0]," has a formatting problem."
        raise ValueError("Error! Portal ",portal[0]," has a formatting problem.")
    
    
def check_portals_array_init_conditions(portals):
    if len(portals) < 3:
        print "Error: Must have more than 2 portals!"
        raise ValueError("Error: Must have more than 2 portals!")
    if len(portals) > _MAX_PORTALS_:
        print "Error: Portal limit is {0}".format(_MAX_PORTALS_)
        raise ValueError("Error: Portal limit is {0}".format(_MAX_PORTALS_))
    seen = set()
    for portal in portals:
        if portal[1] in seen:
            raise ValueError("Error: multiple entries of {0}".format(portal))
        else:
            seen.add(portal[1])
    
def create_portals_array_from_inputfile(input_file):
    # You will need three pieces of information about the portals you want to include: 
    # 1) the portal name, 
    # 2) the Intel map link for this portal (see below for instructions on getting this), 
    # 3) (optional) the number of keys your teams has for this portal, and 
    # 4) (optional) SBLA if the portal is equiped with a SoftBank Link Amp.
    # Thus, each line should be name;intel_link;keys;SBLA
    portals = pd.read_table(input_file,sep=';',
                            comment='#',index_col=False,
                            names=['name','link','keys','sbla'],
                            dtype=str)
    portals = np.array(portals)
    # keep only lines with at least name and url (basestring)
    portals = np.array([portal for portal in portals if (isinstance(portal[0], basestring) and isinstance(portal[1], basestring))])
    print "Found {0} portals in portal list.".format(len(portals))
    #print portals
    check_portals_array_init_conditions(portals)
    return portals
    
def error_portal_formatting(portal):
    print "Error! Portal ",portal[0]," has a formatting problem."
    raise ValueError("Error! Portal ",portal[0]," has a formatting problem.")
    
def extract_loc_from_url(portal):
    url = portal[1]
    if not 'pll=' in url: # this is the URL
        error_portal_formatting(portal)
    coords = (url.split('pll='))
    if len(coords) != 2:
        error_portal_formatting(portal)
    coord_parts = coords[1].split(',')
    lat = int(float(coord_parts[0]) * 1.e6)
    lon = int(float(coord_parts[1]) * 1.e6)
    loc = np.array([lat,lon],dtype=float)
    return loc
    
def extract_key_number(portal):
    keys = 0
    if str(portal[2]) != 'nan':
        if str(portal[2]).lower() == 'a':
            return 1e9
        try: # this is the number of keys
            keys = int(portal[2])
        except ValueError:
            pass
    return keys
    
def extract_sbla(portal):
    sbla = False
    if str(portal[3]) != 'nan':
        try:
            sbla = (portal[3].lower() == 'sbla')
        except ValueError:
            pass
    if sbla:
        print "{0} has SBLA".format(portal[0])
    return sbla

def extract_portal_infos(a, locs, num, portal):
    # loop over columns. Four possibilities:
    # 0. First entry is always portal name
    # 1. contains "pll=" it is the Intel URL
    # 2. contains an intenger, it is the number of keys
    # 3. contains "sbla", it is an SBLA portal
    a.add_node(num)
    a.node[num]['name'] = portal[0]
    loc = extract_loc_from_url(portal)
    locs.append(loc)
    a.node[num]['keys'] = extract_key_number(portal)
    #a.node[num]['plannedlinks'] = 0
    #print "num, key: "+str(num)+", "+str(a.node[num]['keys'])
    a.node[num]['sbla'] = extract_sbla(portal)
    
def init_portals_graph(portals):
    a = nx.DiGraph() # network tool
    a.perim = []
    locs = [] # portal coordinates
    for num,portal in enumerate(portals):
        #print num, portal
        check_portal_formatting_input(portal)
        extract_portal_infos(a, locs, num, portal)
    return a, locs
    
def add_perim_graph(a, coords_xy):
    if not a.quiet:
        print("Computing convex hull ...")
    a.perim = tr.Polygon(geometry.getPerim(coords_xy), a, convex=True)
    
    
def add_coords_portals_graph(a, locs):
    locs = np.array(locs,dtype=float)
    # Convert longitude, latitude coords to  
    locs = geometry.LLtoRads(locs) #radians,
    xyz  = geometry.radstoxyz(locs)  # then to cartesian, 
    xy   = geometry.gnomonicProj(locs,xyz)  # then to gnomonic projection
    #print xy
    n = a.order() # number of nodes
    for i in xrange(n):
        #print("{0} has coords in rads {1} and in xyz {2} and gnomonic {3}".format(a.node[i]['name'], locs[i], xyz[i],xy[i]))
        a.node[i]['geo'] = locs[i]
        a.node[i]['xyz'] = xyz[i]
        a.node[i]['xy' ] = xy[i]
    add_perim_graph(a, xy)
    
def load_portals_graph(input_file, quiet):
    portals = create_portals_array_from_inputfile(input_file)
    a, locs = init_portals_graph(portals)
    a.quiet = quiet
    add_coords_portals_graph(a, locs)
    return a

# return a graph object
def load_portals_list(input_file, quiet):
    if not quiet:
        print "Loading portals list ..."
    # If the input file is a portal list, let's set things up
    if input_file[-3:] != 'pkl':
        a = load_portals_graph(input_file, quiet)
    else:
        with open(input_file,'r') as fin:
            a = pickle.load(fin)
    return a
    
def save_TOPCOM_inputfile(a, output_directory):
    n = a.order() # number of nodes
    xy = [a.node[i]['xy'] for i in xrange(n)]
    puiss = np.absolute(np.floor(np.log10(np.absolute(xy))))
    fracxy = "["
    for i in range(len(xy)):
        if i != 0:
            fracxy = fracxy + ","
        l = "["
        for j in range(len(xy[i])):
            if j != 0:
                l = l + ","
            tenpower = int(10**(puiss[i][0]+9))
            z = int(xy[i][j]*tenpower)
            #l.append(z/tenpower)
            l = l + str(z)+ "/" +str(tenpower)
        l = l + ",1]"
        fracxy = fracxy + l
    fracxy = fracxy + "]"    
    
    #print np.absolute(np.floor(np.log10(np.absolute(xy))))
    #print fracxy
    with open(output_directory+'topcom_input.txt','w') as fout:
        fout.write(fracxy)
    if not a.quiet:
        print "TOPCOM chiro computing ... "
    topcom_inputfile = open(output_directory+'topcom_input.txt')
    topcom_output = subprocess.check_output(['./topcom-0.17.8/src/points2chiro', 'v'], stdin=topcom_inputfile)
    with open(output_directory+'topcom.chiro','w') as topcom_outputfile:
        topcom_outputfile.write(topcom_output)
    
    
def pickle_object_file(a, output_directory, output_file):
    with open(output_directory+output_file,'w') as fout:
        pickle.dump(a,fout)  
    
def save_portals_list(a, output_directory, output_file):
    if output_file is not None:
        pickle_object_file(a, output_directory, output_file)
    else:
        pickle_object_file(a, output_directory, "input.pkl")

def save_best_plan(best_plan, output_directory, output_file):
    pickle_object_file(best_plan, output_directory, output_file)
        
    

def generate_plan_details_and_map(best_plan,output_directory,nagents,useGoogle,api_key,color, skipplot, quiet):
    if not quiet:
        print "Generating plan and map ..."
    best_PP = PlanPrinterMap.PlanPrinter(best_plan,output_directory,nagents,color=color)
    # generate plan details and map
    best_PP.keyPrep()
    best_PP.agentKeys()
    best_PP.planMap(useGoogle=useGoogle,api_key=api_key)
    best_PP.agentLinks()
    # These make step-by-step instructional images
    best_PP.split3instruct(useGoogle=useGoogle)
    if not skipplot:
        best_PP.animate(color, useGoogle=useGoogle)
    plt.close('all')
    return best_PP
    
def plan_optimize(a, optimize, attempts, nagents, suboptimal, start_time):
    mf = Maxfield.Maxfield(a, optimize, attempts, nagents, suboptimal, start_time)
    mf.compute()
    return mf.a, mf.best_time

def display_conclusion(best_PP, best_plan, attempts, plan_input, start_time):
    print ""
    print ""
    print ""
    if not plan_input:
        print "Found best plan after {0} iterations.".format(attempts)
    totalTime = best_plan.walktime+best_plan.linktime+best_plan.commtime
    print "Total time: {0} minutes".format(int(totalTime/60. + 0.5))
    print "Number of portals: {0}".format(best_PP.num_portals)
    print "Number of links: {0}".format(best_PP.num_links)
    print "Number of fields: {0}".format(best_PP.num_fields)
    portal_ap = (125*8 + 500 + 250)*best_PP.num_portals
    link_ap = 313 * best_PP.num_links
    field_ap = 1250 * best_PP.num_fields
    print "AP from portals capture: {0}".format(portal_ap)
    print "AP from link creation: {0}".format(link_ap)
    print "AP from field creation: {0}".format(field_ap)
    print "Total AP: {0}".format(portal_ap+link_ap+field_ap)
    tdiff = time.time() - start_time
    hrs = int(tdiff/3600.)
    mins = int((tdiff-3600.*hrs)/60.)
    secs = tdiff-3600.*hrs-60.*mins
    print "Runtime: {0:02}h {1:02}m {2:05.2f}s".format(hrs,mins,secs)

def main(args):
    start_time = time.time()
    output_directory, output_file, save_portals, input_file  = init_files(args)
    color, useGoogle, api_key, nagents, attempts, suboptimal, topcom, quiet, \
    skipplot, reduce, plan_input, optimize  = init_vars(args)
    a=load_portals_list(input_file, quiet) # graph object
    if plan_input:
        # dont compute best plan 
        best_plan = a
    else:
        save_portals_list(a, output_directory, save_portals)
        if topcom:
            save_TOPCOM_inputfile(a, output_directory)
        best_plan, best_time = plan_optimize(a, optimize, attempts, nagents, suboptimal, start_time)
        if best_plan == None:
            return 1
        if reduce:
            best_time = Maxfield.optimize_best_plan(best_plan, nagents)
        save_best_plan(best_plan, output_directory, output_file)
    
    best_PP = generate_plan_details_and_map(best_plan, output_directory, nagents, useGoogle, api_key, color, skipplot, quiet)
    display_conclusion(best_PP, best_plan, attempts, plan_input, start_time)
    return 0

def get_parsed_args():
    description=("Ingress Maxfield - Maximize the number of links "
                 "and fields, and thus AP, for a collection of "
                 "portals in the game Ingress.")
    parser = argparse.ArgumentParser(description=description,
                                     prog="makePlan.py")
    parser.add_argument('-v','--version',action='version',
                        version="Ingress Maxfield v{0}".format(_V_))
    parser.add_argument('-g','--google',action='store_true',
                        help='Make maps with google maps API. Default: False')
    parser.add_argument('-a','--api_key',default=None,type=str,
                        help='Google API key for Google maps. Default: None')
    parser.add_argument('-n','--num_agents',type=int,default='1',
                        help='Number of agents. Default: 1')
    parser.add_argument('--optimize',type=str,default="MU",
                        help='Want to optimize: MU, AP, Fan, Length (of the path) or Keys. Default: MU')
    parser.add_argument('input_file',type=str,
                        help="Input semi-colon delimited portal file")
    parser.add_argument('-p','--input_plan', action='store_true',
                        help="Input is a pickle object of a plan to analyse.")
    parser.add_argument('-d','--output_dir',default='',type=str,
                        help="Directory for results. Default: "
                        "this directory")
    parser.add_argument('-f','--output_file',default='output.pkl',
                        type=str,
                        help="Filename for pickle object of best plan founded. Default: "
                        "output.pkl")
    parser.add_argument('--save_portals', type=str,default=None,
                        help="Filename for pickle object of portals list.")
    parser.add_argument('-o','--optimal',action='store_true',
                        help='Force optimal solution. Default: False')
    parser.add_argument('--topcom',action='store_true',
                        help='Using TOPCOM external program for optimal solution. Default: False')
    parser.add_argument('-r','--res',action='store_true',
                        help='Use resistance colors. Default: False')
    parser.add_argument('--attempts',type=int,default=_NUM_ATTEMPTS,
                        help='Number of iterations to try new plans. Default: 100')
    parser.add_argument('-q','--quiet',action='store_true',
                        help='Do not display status bar. Default: False')
    parser.add_argument('-s','--skipplot',action='store_true',
                        help='Skip the step-by-step plots. Default: False')    
    parser.add_argument('--noreduce',action='store_false',
                        help='Do not reduce the length of path founded. Default: True')
    parser.add_argument('--timeout',type=float,default=None,help='Timeout in seconds. Default: None')
    parser.add_argument('--log',type=str,default=None,help='Log file. Default: print to screen.')
    return parser.parse_args()

def set_up_job_using_pebble_to_handle_timeout(args):
    with process.Pool(1) as p:
        job = p.schedule(main,args=(args,),timeout=args.timeout)
    try:
        job.get()
    except TimeoutError:
        if args.log is not None:
            sys.stdout = open(args.log,'a',0)
        print "Timeout error: your plan took longer than {0} seconds to complete. Please try submitting again and/or removing some portals.".format(args.timeout)

#Begining of the script !!!!
if __name__ == "__main__":
    args = get_parsed_args()
    # time limitation
    if args.timeout is not None:
        set_up_job_using_pebble_to_handle_timeout(args)
    else:
        #without timeout
        main(args) 
