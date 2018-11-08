#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ingress Maxfield - maxfield.py

GNU Public License
http://www.gnu.org/licenses/
Copyright(C) 2016 by
Jonathan Baker; babamots@gmail.com
Trey Wenger; tvwenger@gmail.com

General code to optimize fielding strategy.

Original version by jpeterbaker
29 Sept 2014 - tvw V2.0 major updates
26 Feb 2016 - tvw v3.0
              merged some new stuff from jpeterbaker's new version
              Added SBLA support
30 Mai 2018 - Christopher Goyet :
            adding max area optimization
            adding delaunay triangulation
            adding Fan triangulation
            adding exhaustive triangulation searching
            adding homogeneous Triangulation searching
"""
import copy
import sys
import time
import bisect
from Triangle import *
import geometry
import agentOrder
np = geometry.np
from scipy.spatial import Delaunay

'''
Some things are chosen randomly:
    Each triangle's splitting portal
    Each triangle's "final node" (unless determined by parent)
This is the number of times to randomly rebuild each first generation
triangle while attempting to get it right
'''
TRIES_PER_TRI = 1

def removeSince(a, m, t):
    # Remove all but the first m edges from a (and .edge_stck)
    # Remove all but the first t Triangules from a.triangulation
    for i in xrange(len(a.edgeStack) - m):
        p,q = a.edgeStack.pop()
        try:
            a.remove_edge(p,q)
        except Exception:
            # The only exception should be from the edge having been
            # reversed
            a.remove_edge(q,p)
            # print 'removing',p,q
            # print a.edgeStack
    while len(a.triangulation) > t:
        a.triangulation.pop()

def compute_time_attributes(a, nagents):
    # By agentOrder.getAgentOrder function ...
    m = a.size()  # number of links
    # if the ith link to be made is (p,q) then orderedEdges[i] = (p,q)
    orderedEdges = [None]*m
    for e in a.edges():
        orderedEdges[a.edges[e[0], e[1]]['order']] = e
    # movements[i][j] is the index (in orderedEdges) of agent i's jth link
    movements = agentOrder.getAgentOrder(a, nagents, orderedEdges)
    
def compute_plan_time(a, nagents):
    compute_time_attributes(a, nagents)
    return a.walktime+a.linktime+a.commtime
    
def optimize_best_plan(a, nagents):
    if not a.quiet:
        print "\nReducing the path length of founded plan ..."
    agentOrder.improveEdgeOrderMore(a) #greedy algorithm to reduce the path length.
    if not a.quiet:
        print "\033[FReducing the path length of founded plan : OK"
    # Re-run to fix the animations and stars of edges that can be done early
    # (improveEdgeOrderMore may have modified the completion order)
    try:
        #first = True
        cleanMarkedEdgesWithFields(a)
        # if not b.quiet:
        #     print "Triangulation="
        #     print [b.triangulation[i].verts for i in xrange(len(b.triangulation))]
        for t in a.triangulation:
            t.markEdgesWithFields() #clean = first)
            #first = False
    except AttributeError:
        print "Error: problem with bestgraph... no triangulation...?"
    best_time = compute_plan_time(a, nagents)
    return best_time    
    
def random_triangulate(b, perim):
    '''
    Recursively tries every triangulation in search a feasible one
        Each layer
            makes a Triangle out of three perimeter portals
            for every feasible way of max-fielding that Triangle
                try triangulating the two perimeter-polygons to the
                sides of the Triangle

    Returns True if a feasible triangulation has been made in graph a
    '''
    #print "Triangulation INPUT : "
    #print perim
    pn = len(perim)
    if pn < 3:
        # Base of recursion
        return True
    try:
        startStackLen = len(b.edgeStack)
    except AttributeError:
        startStackLen = 0
        b.edgeStack = []
    try:
        startTriLen = len(b.triangulation)
    except AttributeError:
        startTriLen = 0
        b.triangulation = []

    # odegrees = [a.out_degree(p) for p in perim
    # order = np.argsort(odegrees)

    # Try all possible first generation triangles with two edges on
    # boundary that both use node i (using i as final vertex will
    # cause no 2 first generation triangles to have same final vertex)
    rand_perim = np.random.permutation(range(0,pn))
    #print "rand perim : ",[perim[i] for i in rand_perim]
    for i in rand_perim:
        #print "i="+str(i)
        # print perim
        # print 'using %s as final'%perim[i]
        for j in xrange(TRIES_PER_TRI):
            #print "j="+str(j)
            t0 = Triangle(perim[[i,i-1,(i+1)%pn]], b, True,
                        suboptimal=False)
            #print t0.verts
            t0.findContents()
            t0.findArea()
            # t0.randSplit() # Split triangle on a random portal
            t0.nearSplit() # Split triangle on the nearest portal
            #print 'trying to build'
            if not t0.buildGraph():
                removeSince(b, startStackLen, startTriLen)
            else:
                break
        else:
            #print "else continue"
            # The loop ended "normally" so this triangle failed
            #print 'big fail'
            continue
        #print 'continuing with',perim[range(i+1-pn,i)]
        if not random_triangulate(b, perim[range(i+1-pn,i)]): # i+1 through i-1
            # remove the links formed since beginning of loop
            #print "not triangulate "+str(perim[range(i+1-pn,i)])
            removeSince(b, startStackLen, startTriLen)
            continue

        # This will be a list of the first generation triangles
        #print "triangulation append : "+str(t0.verts)
        b.triangulation.append(t0)

        # This triangle and the ones to its sides succeeded
        return True

    # Could not find a solution
    return False   


class Maxfield:
    def __init__(self, a, purpose, attempts, nagents, suboptimal, start_time, key_optimize=False):
        self.a = a
        self.purpose = purpose
        self.attempts = attempts
        self.nagents = nagents
        self.suboptimal = suboptimal
        self.start_time = start_time
        self.key_optimize = key_optimize
        
        self.best_area = None
        self.best_time = None
        self.best_triangulation = None
        self.optimized_polygon_done = []
            
    def search_in_optimize_done(self, poly):
        i = bisect.bisect_left(self.optimized_polygon_done, poly)
        if i != len(self.optimized_polygon_done) and poly == self.optimized_polygon_done[i]:
            return i
        return None
           
    def canFlip(self, degrees, keylacks, hasSBLA, p, q):
        '''
        True if reversing edge p,q is a paraeto improvement
            out-degree of q must be <8
            p must have a key surplus
            or if portal has SBLA, < 40 links, and key surplus
        '''
        case1 = (degrees[q,1] < 8) & (keylacks[p]<0)
        case2 = (hasSBLA[q]) & (degrees[q,1] < 40) & (keylacks[p]<0)
        return case1 | case2
    
    def flip(self, p, q, degrees=None, keylacks=None):
        if not self.a.edges[p,q]['reversible']:
            print '!!!! Trying to reverse a non-reversible edge !!!!'
            print p,q
        # Give the reversed edge the same properties
        #a.add_edge(q,p,a.edges[p,q])
        self.a.add_edge(q,p)
        self.a.edges[q,p].update(self.a.edges[p,q])
        self.a.remove_edge(p,q)
        if degrees is not None:
            degrees[p,0] += 1
            degrees[p,1] -= 1
            degrees[q,0] -= 1
            degrees[q,1] += 1
    
        if keylacks is not None:
            keylacks[p] += 1
            keylacks[q] -= 1
    
    def flipSome(self):
        '''
        Tries to reduce the number of keys that need to be farmed by reversing edges
        Only edges with the property reversible=True will be flipped 
        '''
        n = self.a.order()
        degrees  = np.empty([n,2],dtype=int)
        keylacks = np.empty(n,dtype=int) # negative if there's a surplus
        hasSBLA = np.empty(n,dtype=bool)
    
        # column 0 is in-degree, col 1 is out-degree
        for i in xrange(n):
            degrees[i,0] = self.a.in_degree(i)
            degrees[i,1] = self.a.out_degree(i)
            keylacks[i] = degrees[i,0]-self.a.node[i]['keys']
            hasSBLA[i] = self.a.node[i]['sbla']
    
        # This is now commented out because plans are not submitted to
        # this function without obeying the 8 outgoing links limit.
        """
        # We can never make more than 8 outogoing links. Reducing these is
        # first priority
        manyout = (degrees[:,1]>8).nonzero()[0]
        for p in manyout:
        qs = list(a.edge[p].iterkeys())
        for q in qs:
            if a.edge[p][q]['reversible'] and canFlip(degrees,keylacks,p,q):
                flip(a,p,q,degrees,keylacks)
            if degrees[p,1] <= 8:
                break
        else:
            # This runs if the for loop exits without the break
            print 'Could not reduce OUT-degree sufficiently for %s'%p
        """
    
        # It is difficult to gather more keys. Reducing key-gathering is next
        # priority
        # We'll process the ones with the greatest need first
        needkeys = (keylacks>0).nonzero()[0]
        needkeys = needkeys[np.argsort(keylacks[needkeys])][::-1]
        for q in needkeys:
            for p,q2 in self.a.in_edges(q):
                if self.a.edges[p,q]['reversible'] and self.canFlip(degrees,keylacks,hasSBLA,p,q):
                    self.flip(p,q,degrees,keylacks)
                if keylacks[q] <= 0:
                    break
            #else:
                # This runs if the for loop exits without the break
                # TODO : May be suggest a SBLA ?
                # print 'Could not reduce IN-degree sufficiently for %s'%q
            
    def display_progression_infos(self, num, attempts, best_time, newline=False):
        if not self.a.quiet:
            tdiff = time.time() - self.start_time
            hrs = int(tdiff/3600.)
            mins = int((tdiff-3600.*hrs)/60.)
            secs = tdiff-3600.*hrs-60.*mins
            sys.stdout.write("\033[K[{0:20s}] {1:2}% ({2}/{3} iterations) : {4:02}h {5:02}m {6:05.2f}s (shortest path : {7})".\
                            format('='*(20*num/attempts),
                                100*num/attempts,
                                num,attempts,
                                hrs,mins,secs, 
                                str(int(best_time/60. + 0.5))+" min" if best_time<1.e9 else "infinity"
                                ))
        if newline:
            print("")
        
    '''
    ********************************************************************************
    Default triangulation : 
        it seems they wanted to reduce the length of the path
    '''
        
    def try_another_path(self):
            b = copy.deepcopy(self.a)
            #cleanMarkedEdgesWithFields(b)
            random_triangulate(b, self.a.perim)
            agentOrder.improveEdgeOrder(b)
            totalTime = compute_plan_time(b, self.nagents)
            return b, totalTime
        
    def shorter_path_attemp(self, best_plan, best_time):
            b, totalTime = self.try_another_path()
            if not self.suboptimal:
                totalTime = optimize_best_plan(b, self.nagents)
                if not self.a.quiet:
                    sys.stdout.write("\033[K\033[F\033[K\033[F\033[K")
            if totalTime < best_time:
                #print "change"
                return b, totalTime
                #best_plan = b
                #best_PP = copy.deepcopy(PP)
                #best_time = totalTime
            return best_plan, best_time
    
    def plan_shortest_walking_distance(self):
        # try serveral plans to get shortest walking distance
        best_plan = self.a
        best_time = 1.e9
        for foobar in xrange(self.attempts):
            self.display_progression_infos(foobar, self.attempts, best_time, newline=False)
            best_plan, best_time = self.shorter_path_attemp(best_plan, best_time)
        #display_progression_infos(attempts, attempts, best_time, a.quiet, start_time, newline=True)
        final_tri = best_plan.triangulation[:]
        removeSince(best_plan, 0, 0)
        self.a = best_plan
        return final_tri, best_time   
        
    '''
    ********************************************************************************
        Triangulation where we try to maximize area and then MU
    
    '''
        
    def best_area_under(self, triangle, quiet):
        isBestSplitDone = self.search_in_optimize_done(triangle)
        if isBestSplitDone != None:
            #print "triangle already optimize :",triangle.verts, triangle
            triangle = self.optimized_polygon_done[isBestSplitDone]
            return triangle
        if triangle.area == None:
            triangle.findArea() #Already done?
        triangle.findContents()
        # if self.key_optimize:
        #     print "\nSplit to optimize keys."
        #     area0 = triangle.area + triangle.maxKEY_Split()
        # else:
        if self.suboptimal:
            triangle.multi_area = triangle.area + triangle.near_maxAREA_Split() # Split triangle on the nearest portal which maximize area
        else:
            if not quiet and not self.a.quiet:
                print "\nSearching optimal split of triangle {} :".format(triangle)
            triangle.multi_area = triangle.area + \
                triangle.maxAREA_Split_optimize(self.optimized_polygon_done, 
                                        quiet=(quiet or self.a.quiet) )
            if not quiet and not self.a.quiet:
                sys.stdout.write("\033[K\033[F\033[K\033[F")
        if len(triangle.contents) > 1 :
            #print "Optimized triangle added : ", triangle
            bisect.insort_left(self.optimized_polygon_done, triangle)
        return triangle
    
    def triangulation_maxArea_on_edges(self, triangle, polygon):
        new_perims = polygon.split_on_subpoly(triangle, convex=True)
        #print "triangle = ",triangle, "   poly = ", polygon
        #print "new_perims = ", new_perims
        for l in new_perims:
            area = triangle.multi_area
            current_tri = [triangle]
            for side_poly in l:
                #print "side_poly = ",side_poly
                sub_polygon = self.search_max_area_triangulation(side_poly, base_first=True, quiet=True)
                #print "sub_polygon = ",sub_polygon, " triangulation", sub_polygon.triangulation
                current_tri = current_tri + sub_polygon.triangulation
                area = area + sub_polygon.multi_area
            if area>polygon.multi_area:
                polygon.multi_area = area
                polygon.triangulation = current_tri
                
    def polygon_max_area_triangulation(self,  polygon, base_first=False, quiet=True):
        candidates = triangles_in_perim_sorting_by_area(self.a, polygon.verts, self.suboptimal , base_first=base_first, quiet=quiet)
        candidates = candidates[0:self.attempts]
        #print "candidates = ",candidates, " attempts = ", self.attempts
        attempts = len(candidates) # in case  attempts > len candidates
        tour = 0
        display_message(self.a, quiet, "Searching ...",  newline=True)
        for triangle in candidates:
            display_progression_triangle(self.a, self.a.quiet, tour, attempts, "attemps, optimized="+str(len(self.optimized_polygon_done)))
            triangle = self.best_area_under(triangle, quiet=False)
            display_message(self.a, self.a.quiet, "\n extending triangulation from edges ...",  newline=True)
            self.triangulation_maxArea_on_edges(triangle, polygon)
            if not self.a.quiet:
                sys.stdout.write("\033[K\033[F\033[K\033[F\033[K")
            tour += 1
        display_progression_triangle(self.a, quiet, tour, attempts, "attemps, optimized="+str(len(self.optimized_polygon_done)))
    
    def search_max_area_triangulation(self, polygon, base_first=False, quiet=True):
        polygon.multi_area = 0
        polygon.triangulation = []
        if len(polygon) <= 2:
            return polygon
        """ Is polygon known ? """
        isBestSplitDone = self.search_in_optimize_done(polygon)
        if isBestSplitDone != None:
            polygon = self.optimized_polygon_done[isBestSplitDone]
            #print "Polygon already optimize :",polygon.triangulation, polygon.multi_area
            return polygon
        #print "Polygon not optimized :",polygon
        if len(polygon) == 3:
            triangle = Triangle(polygon.verts, self.a, exterior=True, suboptimal=self.suboptimal)
            triangle.findArea()
            triangle = self.best_area_under(triangle, quiet=False)
            return triangle
        self.polygon_max_area_triangulation(polygon, base_first=base_first, quiet=quiet)
        #print "Optimized polygon added : ", polygon
        bisect.insort_left(self.optimized_polygon_done, polygon)
        return polygon
        
    '''
    ********************************************************************************
        Triangulation where we try to plan classical multis
    '''
    def best_multi_under(self, triangle, quiet):
        isBestSplitDone = self.search_in_optimize_done(triangle)
        if isBestSplitDone != None:
            #print "triangle already optimize :",triangle.verts, triangle
            triangle = self.optimized_polygon_done[isBestSplitDone]
            return triangle
        if triangle.area == None:
            triangle.findArea() #Already done?
        triangle.findContents()
        # if self.key_optimize:
        #     print "\nSplit to optimize keys."
        #     area0 = triangle.area + triangle.maxKEY_Split()
        # else:
        if self.suboptimal:
            triangle.multi_area = triangle.area + triangle.near_maxMULTI_Split() # Split triangle on the nearest portal which maximize area
        else:
            if not quiet and not self.a.quiet:
                print "\nSearching optimal split of triangle {} :".format(triangle)
            triangle.attempts = self.attempts
            triangle.multi_area = triangle.area + \
                triangle.maxMULTI_Split_optimize(self.optimized_polygon_done, 
                                        quiet=(quiet or self.a.quiet) )
            if not quiet and not self.a.quiet:
                sys.stdout.write("\033[K\033[F\033[K\033[F\033[K")
        if len(triangle.contents) > 1 :
            #print "Optimized triangle added : ", triangle
            bisect.insort_left(self.optimized_polygon_done, triangle)
        return triangle
    
    def select_base_first(self, polygon):
        polygon.findContents()
        candidates = polygon.verts + polygon.contents
        list_keys = np.array([self.a.nodes[i]['keys'] for i in candidates])
        maxkey_point = np.argsort(-list_keys)
        base_verts = [candidates[maxkey_point[0]], candidates[maxkey_point[1]] ]
        base = Polygon(base_verts, self.a)
        new_perims = polygon.split_on_subpoly(base, convex=True)
        return new_perims
        
    def search_max_MULTI_triangulation(self, polygon, base_first=False, quiet=True):
        polygon.multi_area = 0
        polygon.triangulation = []
        if len(polygon) <= 2:
            return polygon
        if len(polygon) == 3:
            triangle = Triangle(polygon.verts, self.a, exterior=True, suboptimal=self.suboptimal)
            triangle.findArea()
            triangle = self.best_multi_under(triangle, quiet=False)
            return triangle
        base_perim = self.select_base_first(polygon)
        candidates = []
        for l in base_perim:
            for poly in l:
                sub_candidates = triangles_in_perim_sorting_by_area(self.a, poly.verts, self.suboptimal, base_first=True, quiet=quiet)
                sub_candidates = sub_candidates[0:self.attempts]
                candidates += sub_candidates
        attempts = len(candidates) # in case  attempts > len candidates
        tour = 0
        display_message(self.a, quiet, "Searching ...",  newline=True)
        for triangle in candidates:
            triangle.verts = list(np.roll(triangle.verts, 1))
            display_progression_triangle(self.a, self.a.quiet, tour, attempts, "attemps, optimized="+str(len(self.optimized_polygon_done)))
            triangle = self.best_multi_under(triangle, quiet)
            if triangle.multi_area > polygon.multi_area:
                polygon.triangulation = [triangle]
                polygon.multi_area = triangle.multi_area
            tour += 1
        display_progression_triangle(self.a, self.a.quiet, tour, attempts, "attemps, optimized="+str(len(self.optimized_polygon_done)))
        return polygon
    
    '''
    ********************************************************************************
        
        Bubble Triangulation :
        Triangulation where we try to optimize keys
    '''
    
    
    
    def BUBBLE_centers(self, perim, contents, base_first=False):
        perim = perim.verts
        a = self.a
        perim_xyz = [a.nodes[i]['xyz'] for i in perim]
        remain = list(perim)+[i for i in contents \
                if (not i in perim) and geometry.spherePolygonContains(perim_xyz, a.nodes[i]['xyz']) ]
        #print "BUBBLE_centers = ",perim, remain
        # choose max keys portal
        list_keys = np.array([a.nodes[i]['keys'] for i in remain])
        maxkey_point = np.argmax(list_keys)
        #print "max keys portal : ",remain[maxkey_point]," numkey = ",list_keys[maxkey_point]
        # compute distances of other points with him
        args_closest = nearest_to_X_in_Y(a, remain[maxkey_point], remain)[0]
        #print "args_closest : ",args_closest
        # keep only the first we can link
        args_closest = args_closest[0:list_keys[maxkey_point]+8] # 8 out_link and in with keys
        #print "args_closest : ",args_closest
        closest_pts = [remain[i] for i in args_closest]
        #print "closest_pts = ",closest_pts
        # get perim of this area :
        pts = np.array([ a.node[i]['xy'] for i in closest_pts ])
        args_new_perim = np.array(geometry.getPerim(pts))
        #print "args_new_perim : ",args_new_perim
        new_perim = [closest_pts[i] for i in args_new_perim]
        if base_first and (perim[0] in new_perim) and (perim[1] in new_perim):
            for i in range(1,len(new_perim)):
                if (perim[0] == new_perim[0]) and (perim[1] == new_perim[1]):
                    break
                new_perim = new_perim[-1:] + new_perim[:-1] #rotate -1
        return np.setdiff1d(remain, closest_pts), new_perim, remain[maxkey_point]
        
        #split = split_perim(a, perim, new_perim)
        #print "split = "
        #print split
        #reste = np.setdiff1d(remain, closest_pts)
        #if len(reste)>0:
        #    BUBBLE_centers(a, np.append(perim, reste ))
    
    def fan_triangulation(self, perim, center, base_first=False):
        print "perim :", perim
        print "center : ", center
        n = len(perim)
        tri = []
        if center in perim:
            arg_center = perim.index(center)
                    #fan_perim = np.roll(perim,-arg_center)
            if arg_center>1:
                for i in xrange(n):
                    if (i != arg_center) and ((i+1)%n != arg_center ):
                        tri.append(Triangle([perim[i], perim[(i+1)%n], center],
                            self.a,True,suboptimal=self.suboptimal))
            else:
                if arg_center == 0:
                    for i in xrange(1,n-1):
                        tri.append(Triangle([center, perim[i%n], perim[(i+1)%n]],
                            self.a,True,suboptimal=self.suboptimal))
                else: #arg_center == 1:
                    for i in xrange(n-2):
                        tri.append(Triangle([perim[-i], center, perim[-i-1]],
                            self.a,True,suboptimal=self.suboptimal))
        else:
            for i in xrange(n):
                tri.append(Triangle([perim[i], perim[(i+1)%n], center],
                        self.a,True,suboptimal=self.suboptimal))
        # Test because buildgraph return false
        # if base_first and not (tri[0][0] in perim[0:2] and tri[0][1] in perim[0:2]):
        #     print "reverse"
        #     tri.reverse()
        print "tri :", [t.verts for t in tri]
        return tri
    
    def search_FAN_triangulation(self, perim, remains_pts, base_first=False):
        print "\nsearch_FAN_triangulation : ", perim
        print "\n"
        area = 0
        final_tri = []
        if len(perim)<=2:
            return final_tri, area
        new_remain, new_perim, center = self.BUBBLE_centers(perim, remains_pts, base_first=base_first)
        #print "new_perim = ", new_perim
        #print "center = ", center
        fan_tri = self.fan_triangulation(new_perim, center, base_first=base_first)
        #print "fan_triangulation =", fan_tri
        for triangle in fan_tri:
            #t0 = Triangle([triangle[0], triangle[1], triangle[2]],a,True,
            #           suboptimal=suboptimal)
            triangle.findContents()
            triangle.findArea() 
            final_tri += [triangle]
            #sub_area = triangle.near_maxAREA_Split()
            sub_area = triangle.fan_Split()
            area += triangle.area + sub_area
        # Next bubbles :
        if len(perim) == 3:
            return final_tri, area
        #print "perim : ",perim
        #print "polygon : ",new_perim
        # bubbles = split_perim(a, perim, new_perim)[0] #others splits not good
        # print "bubbles : "
        # print bubbles
        # for b in bubbles:
        #     sub_tri, sub_area = search_BUBBLE_triangulation(a, b, new_remain, attempts, suboptimal, base_first=True)
        #     #if len(sub_tri)>0:
        #     #   sub_tri.reverse()
        #     final_tri += sub_tri
        #     area += sub_area
        return final_tri, area
        
    '''
    ********************************************************************************
    
        Base on Delaunay Triangulation :
        Triangulation where we try to optimize keys
        
    
    '''
    
    def sub_Delaunay_triangluation_in_perim(self, perim):
        n = self.a.order()
        perim_xyz = [self.a.nodes[i]['xyz'] for i in perim]
        # all points in but not on perim 
        list_pts = [i for i in range(n) \
                if (not i in perim) and geometry.spherePolygonContains(perim_xyz, self.a.nodes[i]['xyz']) ]
        list_pts_keys = np.array([self.a.nodes[i]['keys'] for i in list_pts])
        args_pts = np.argsort(-list_pts_keys)
        pts_sorted = [list_pts[i] for i in args_pts]
        #print [a.nodes[list_pts[i,0]]['name'] for i in points2add]
        perim_xy = [self.a.nodes[i]['xy'] for i in perim]
        n = len(list_pts)
        #tri = [Delaunay(perim_xy)]
        tri = []
        for i in xrange(1,n):
            added_points_xy = [ self.a.nodes[ pts_sorted[j] ]['xy'] for j in range(i)]
            tri.append(Delaunay(perim_xy+added_points_xy))
        print pts_sorted
        return tri, list(perim)+pts_sorted
                
    
    def search_optimize_keys_triangulation(self, perim):
        #print "search_optimize_keys_triangulation:", perim
        best_area = 0
        final_tri = []
        if len(perim)<=2:
            return final_tri, best_area
        
        delaunay_args, delaunay_points = self.sub_Delaunay_triangluation_in_perim(perim)
        #print "preferences: ",preferences #[0:attempts]
        for tri in delaunay_args[0:attempts]:
            current_tri = []
            area = 0
            linksGraph = self.a.to_undirected()
            for triangle in tri.simplices:
                t0 = Triangle([delaunay_points[triangle[0]], delaunay_points[triangle[1]], delaunay_points[triangle[2]]],
                                self.a,True,suboptimal=self.suboptimal)
                links_Triangle(t0, linksGraph)
                t0.findContents() 
                current_tri += [t0]
            for t in current_tri:
                #print "t0.area=", t0.area
                # TODO : Not always optimal triangulation
                # TODO : add optimal option to maxAREA_Split ?
                sub_area, sub_added_linked = t.maxKEY_Split(linksGraph)
                area += t.area + sub_area
            if area>best_area or True:
                #if best_area >0:
                    #print "best found !"
                    #print "Old tri :",[tri.verts for tri in final_tri]," old area =", best_area
                    #print "New tri :",[tri.verts for tri in current_tri]," new area =", area
                best_area = area
                final_tri = current_tri
        return final_tri, best_area
        
    '''
    ********************************************************************************
        
        Homogeneous Triangulation searching
        
    '''
    
    
    def search_max_homogenous_triangulation(self, perim, base_first=False):
        final_tri = None
        best_degree = 0
        if len(perim)<=2:
            return final_tri
        candidates = triangles_homo_sorting_by_contents(self.a, perim, self.attempts, self.suboptimal, base_first=base_first)
        candidates = candidates[0:self.attempts]
        attempts = len(candidates) # case attempts < len(candidates)
        tour=0
        for triangle in candidates:
            if not self.a.quiet:
                sys.stdout.write("\r[{0:20s}] {1:2}% ({2}/{3} attemps) Best degree : {4} ".\
                                    format('#'*(20*tour/attempts),100*tour/attempts,tour,attempts,best_degree))
                sys.stdout.flush()
            degree = 0
            sys.stdout.write("\n")
            degree = triangle.split_homogenous(self.suboptimal, mindeg=best_degree, quiet=False)
            if not self.a.quiet:
                sys.stdout.write("\033[K\033[F")
            if degree > best_degree:
                best_degree = degree
                triangle.findArea()
                final_tri = triangle
            elif degree == best_degree:
                triangle.findArea()
                if (final_tri == None)  or (final_tri.area < triangle.area):
                    final_tri = triangle
            tour += 1
        sys.stdout.write("\r[{0:20s}] {1:2}% ({2}/{3} attemps) Best degree : {4}\n".\
                                    format('#'*(20*tour/attempts),100*tour/attempts,tour,attempts,best_degree))
        print "\nfinal_tri=",final_tri.verts
        return [final_tri]
        
        
    '''
    ********************************************************************************
        Exhaustive search triangulation : 
    
            it can optimize all you want, if you can wait enough ...
        
    '''
    
    def perim_edges(self, G, p):
        '''
        Returns the list of edges to make the cycle p
        with edge sorted by points num ascendant
        '''
        n = len(p)
        perim = []
        for i in range(n):
            if p[i-1]<p[i]:
                perim.append((p[i-1],p[i]))
            else:
                perim.append((p[i],p[i-1]))
        return perim
            
    def set_possibles_edges(self, G, perim):
        #unused_perim_nodes = set([n for n in G.neighbors(starting_node)].append(starting_node))
        #nodes_list = set(G.nodes).difference(unused_perim_nodes)
        #edges_list = [(starting_node, n), for n in nodes_list]
        n = len(G.nodes)
        edges_list = []
        for p in xrange(n):
            for q in xrange(p+1, n):
                if q>p:
                    edges_list.append((p, q))
                    dist = geometry.norms2(G.nodes[q]['xy']-G.nodes[p]['xy'])
                    G.add_edge(p,q,norm= dist)
        edges_set = set(edges_list).difference(set(perim))
        return edges_set
        
    def intersect_previous_link(self, G, l, pile_done):
        #print "l="+str(l)
        A = G.node[l[0]]['xy']
        B = G.node[l[1]]['xy']
        for previous_link in pile_done:
            #print str(l)+" ???? "+str(previous_link)
            if (l[0] == previous_link[0]) or (l[0] == previous_link[1]) \
            or (l[1] == previous_link[0]) or (l[1] == previous_link[1]):
                #print "extremités égales - ne croisent pas : "+str(l)+str(previous_link)
                continue
            S = G.node[previous_link[0]]['xy']
            T = G.node[previous_link[1]]['xy']
            if geometry.isIntersect(A, B, S, T):
                #print True
                return True
        #print False
        return False
            
        
    def triangulations_rec(self, G, L, perim):
        h = len(perim)
        n = len(G.nodes)
        #i=0
        #j=0
        total_liens = 3*n-2*h-3 # number edges in a fine graph
        print "total links = "+str(total_liens+h)
        triangulated_graphs = []
        pile_to_try = [list(L)]
        pile_done = []
        while len(pile_to_try) > 0:
            # if len(pile_to_try) >1:
            #     i = len(pile_to_try[1])
            # else:
            #     i= len(L)-2
            # if len(pile_to_try) >2:
            #     j = len(pile_to_try[2])
            # else: 
            #     j = len(L)-3
            #if not G.quiet:
            #    sys.stdout.write( "\r {0} - {1} - {2} / {3}".format(len(pile_to_try[0]), i, j, len(L)) )
            #    sys.stdout.flush()
            list_edges = pile_to_try.pop()
            while len(list_edges) > 0 and len(pile_done)<total_liens and len(list_edges)+len(pile_done)>=total_liens:
                l = list_edges.pop()
                if not intersect_previous_link(G, l, pile_done):
                    pile_done.append(l)
                    pile_to_try.append(copy.copy(list_edges))
                #print "pile_done = "+str(pile_done)
            if len(list_edges)== 0 and len(pile_done)<total_liens and len(list_edges)+len(pile_done)>=total_liens:
                print "condition utile !!!"
            if len(pile_done)>=total_liens:
                triangulated_graphs.append(copy.deepcopy(pile_done))
                if not G.quiet:
                    sys.stdout.write( "\r {0} triangulations so far.".format(len(triangulated_graphs)))
                    # format('
                    # sys.stdout.write( "\r [{0:20s}] {1}% - [{0:20s}] {3}% - [{0:20s}] {5}% / {6} - Founded : {7}".\
                    # format('#'*(20*(len(L)-1-len(pile_to_try[0]))/(len(L)-1)),
                    #         100*(len(L)-1-len(pile_to_try[0]))/(len(L)-1),
                    #         '#'*(20*(len(L)-2-i)/(len(L)-2)),
                    #         100*(len(L)-i-2)/(len(L)-2),
                    #         '#'*(20*(len(L)-j-3)/(len(L)-3)),
                    #         100*(len(L)-j-3)/(len(L)-3),
                    #         len(L),
                    #         len(triangulated_graphs) ))
                    sys.stdout.flush()
                #print " Triangulations founded :"+str(len(triangulated_graphs))
            if len(pile_done)!=0:   
                l = pile_done.pop()  
        print " "
        return triangulated_graphs    
        
    def triangulations_list(self, a, perim_pts):
        G = a.to_undirected()
        G.quiet = a.quiet
        #nx.set_edge_attributes(G, 0, 'norm')
        #for n in G.nodes:
        #    print str(n)+" = "+str(G.nodes[n]['xy'])
        perim_links = perim_edges(G, perim_pts)
        print "perim_links="+str(perim_links)
        possibles_edges = list(set_possibles_edges(G, perim_links))
        possibles_edges.sort( key=lambda e: G.edges[e[0], e[1]]['norm'], reverse=True)
        print "len(possibles_edges) ="+str(len(possibles_edges))
        tri_list = triangulations_rec(G, possibles_edges, perim_links)
        print "len(tri_list)="+str(len(tri_list))
        print tri_list
        return tri_list
        
    '''
    ********************************************************************************
        Computing triangulation :
        
    '''
    def plan_optimize(self):
        optimize = self.purpose
        if optimize == "KEY" or optimize == "KEYS":
            final_tri, best_area = self.search_optimize_keys_triangulation(self.a.perim)
        elif optimize == "BF":
            tri_list = self.triangulations_list(self.a.perim)
            return len(tri_list) == 0
        elif optimize == "MU":
            polygon = self.search_max_area_triangulation(self.a.perim, quiet=self.a.quiet)
            final_tri = polygon.triangulation
            best_area = polygon.multi_area
        elif optimize == "FAN" or optimize == "BUBBLE" or optimize == "BUBBLES":
            final_tri, best_area = self.search_FAN_triangulation(self.a.perim, range(len(self.a.nodes)))
        elif optimize == "LENGTH":
            final_tri, best_time = self.plan_shortest_walking_distance()
            best_area = np.sum([tri.get_Total_Area() for tri in final_tri])
        elif optimize == "HOMO" or optimize == "HOMOGENEOUS":
            final_tri  = self.search_max_homogenous_triangulation(self.a.perim)
            best_area = final_tri[0].get_Total_Area()
        elif optimize == "MULTI" or optimize == "MULTIS":
            polygon = self.search_max_MULTI_triangulation(self.a.perim, quiet=self.a.quiet)
            final_tri = polygon.triangulation
            best_area = polygon.multi_area
        else:
            print "Error : unknown optimization !"
            sys.exit(0)
        self.best_area = best_area
        if not self.a.quiet:
            print("\n\033[KBest area : {:,} m^2".format(self.best_area).replace(',', ' '))
            print("Estimated MU : {:,} MU".format(int(round(self.best_area/500))).replace(',', ' '))
        self.best_triangulation = final_tri
        return True
    
    def build_triangulate(self):
        '''
        Recursively tries every triangulation in search a feasible one
        If suboptimal, then first generation triangles have all their vertex on the convex hull 
        Returns True if a feasible triangulation has been made in graph a
        '''    
        if not self.a.quiet:
            print "Start of Triangulation Computation"
        if not self.plan_optimize():
            return False
        if not self.a.quiet:
            sys.stdout.write("\rBuild Graph : ")
            sys.stdout.flush()
        build = [t.buildGraph() for t in self.best_triangulation]
        if np.all(build):
            if not self.a.quiet:
                print "OK"
                #print [t.verts for t in self.best_triangulation]
            self.a.triangulation = self.best_triangulation
            return True
        return False
    
    
    def compute(self):
        if self.a.perim == None:
            print "ERROR : no convex hull ..."
            return False 
        if not self.build_triangulate():
            print "Impossible to build_triangulate_maxMU ?"
            return False
        self.flipSome()
        # Attach to each edge a list of fields that it completes
        # catch??? no triangulation (bad portal file?)
        for t in self.a.triangulation:
            t.markEdgesWithFields()
        return True