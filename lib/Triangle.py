#!/usr/env python
# -*- coding: utf-8 -*-
"""
Ingress Maxfield - PlanPrinterMap.py

GNU Public License
http://www.gnu.org/licenses/
Copyright(C) 2016 by
Jonathan Baker; babamots@gmail.com
Trey Wenger; tvwenger@gmail.com

Builds valid fields

original version by jpeterbaker
29 Sept 2014 - tvw V2.0 major updates
26 Feb 2016 - tvw v3.0
              merged some new stuff from jpeterbaker's new version
              Added SBLA support
30 Mai 2018 - Goyet Christopher
               Added : self.area attribut and splitOn return now sum area
               Added near_maxAREA_Split
               Added maxAREA_Split_optimize
               Added splitOn_countlinks
               Added maxKEY_Split
               Added split_homogenous
"""
import geometry
import numpy as np
import sys
import bisect

# Set to False if only perfectly optimal plans should be produced
_suboptimal = True # default value
        
def cleanMarkedEdgesWithFields(a):
    for p,q in a.edges():
        a.edges[p,q]['fields'] = []
        a.edges[p,q]['depends'] = []

def can_add_more_links_from_portal(a, p):
    has_sbla = a.nodes[p]['sbla']
    extended_by_sbla = (has_sbla and a.out_degree(p) < 40)
    return (a.out_degree(p) < 8) or extended_by_sbla

def try_reduce_out_degree(a,p):
    # Reverse as many edges out-edges of p as possible
    # now with SBLA support!
    toremove = []
    for q in a.neighbors(p):
        if can_add_more_links_from_portal(a,q):
            #print dict(a.edges[p,q])
            #print p
            #print q
            #print a.edges[p,q]['reversible']
            if a.edges[p,q]['reversible']:
                a.add_edge(q,p)
                a.edges[q, p].update(a.edges[p, q])
                toremove.append(q)
    for q in toremove:
        a.remove_edge(p,q)

def try_ordered_edge(a,p,q,reversible,suboptimal):
    #print "try_ordered_edge: p=",p,", q=",q
    # now with SBLA support
    if a.has_edge(p,q) or a.has_edge(q,p):
        #print "try_ordered_edge: p=",p,", q=",q," already done !!!"
        return True
    
    # if reversible and a.out_degree(p) > a.out_degree(q):
        # p,q = q,p

    if not can_add_more_links_from_portal(a,p):
        try_reduce_out_degree(a,p)

    if not can_add_more_links_from_portal(a,p):
    # We tried but failed to reduce the out-degree of p
        if not reversible and not suboptimal:
            print '%s already has 8 outgoing'%p
            #raise(Deadend('%s already has max outgoing'%p))
            return False
        if not can_add_more_links_from_portal(a,q):
            try_reduce_out_degree(a,q)
        if (not can_add_more_links_from_portal(a,q) and not suboptimal):
            print '%s and %s already have 8 outgoing'%(p,q)
            #raise(Deadend('%s and %s already have max outgoing'%(p,q)))
            return False
        p,q = q,p
    
    m = a.size()
    #print "add edge : ",(p,q)," order=",m
    a.add_edge(p,q,order=m,reversible=reversible,fields=[],depends=[])

    try:
        a.edgeStack.append( (p,q) )
    except AttributeError:
        a.edgeStack = [ (p,q) ]
        # print 'adding',p,q
        # print a.edgeStack
    return True
    
def links_Triangle(t, b):
    added_links = []
    if not b.has_edge(t.verts[0], t.verts[1]):
        b.add_edge(t.verts[0], t.verts[1])
        added_links.append([t.verts[0], t.verts[1]])
    if not b.has_edge(t.verts[0], t.verts[2]):
        b.add_edge(t.verts[0], t.verts[2])
        added_links.append([t.verts[0], t.verts[2]])
    if not b.has_edge(t.verts[1], t.verts[2]):
        b.add_edge(t.verts[1], t.verts[2])
        added_links.append([t.verts[1], t.verts[2]])
    return added_links

triangleContentCache = {}

def display_message(a, quiet, message,  newline=False):
    if newline:
        message += "\n"
    if not quiet and not a.quiet:
        sys.stdout.write(message)
        sys.stdout.flush()

def display_progression_triangle(a, quiet, tour, total, title_search, newline=False):
    str = "\r[{0:20s}] {1:2}% ({2}/{3} "+title_search+" ) "
    if newline:
        str += "\n"
    if not quiet and not a.quiet :
        sys.stdout.write(str.format('='*(20*tour/total),100*tour/total,tour,total))
        sys.stdout.flush()

def all_triangles_in_perim_suboptimal(a, perim, base_first=False):
    listTriangles = []
    pn = len(perim)
    if base_first:
        for r in xrange(2, pn):
            #triangle oriented in same direction than perim
            t = Triangle([perim[0], perim[1], perim[r]], a, True,
                suboptimal=True)
            listTriangles.append( t )
    else:
        for p in xrange(pn):
            for q in xrange(p+1, pn):
                for r in xrange(q+1, pn):
                    #triangle oriented in same direction than perim
                    t = Triangle([perim[p], perim[q], perim[r]], a, True,
                        suboptimal=True)
                    listTriangles.append( t )
    return listTriangles
    
def all_triangles_in_perim_optimal(a, perim, quiet, base_first=False):
    listTriangles = []
    n = a.order()
    perim_xyz = [a.nodes[i]['xyz'] for i in perim]
    # choose points on perim or  in perim
    list_pts = [i for i in range(n) if (i in perim) or geometry.spherePolygonContains(perim_xyz, a.nodes[i]['xyz']) ]
    n = len(list_pts)
    # print "perim, list_pts",perim, list_pts
    if base_first:
        tour = 0
        for r in list_pts:
            display_progression_triangle(a, quiet, tour, n, "vertex")
            if (r != perim[0]) and (r != perim[1]):
                pt3 = a.nodes[r]['geo']
                t = Triangle([perim[0], perim[1], r], a, True,
                        suboptimal=False)
                listTriangles.append( t )
            tour += 1
        display_progression_triangle(a, quiet, tour, n, "vertex", newline=True)
    else:
        for p in xrange(n):
            display_progression_triangle(a, quiet, p, n, "first vertex")
            pt1_xyz = a.nodes[list_pts[p]]['xyz']
            for q in xrange(p+1, n):
                pt2_xyz = a.nodes[list_pts[q]]['xyz']
                for r in xrange(q+1, n):
                    pt3_xyz = a.nodes[list_pts[r]]['xyz']
                    #triangle oriented in anti-clockwise direction
                    if geometry.spherePolygonContains([pt1_xyz, pt2_xyz], pt3_xyz, close=False) :
                        t = Triangle([p, q, r], a, True,
                                suboptimal=False)
                    else:
                        t = Triangle([p, r, q], a, True,
                                suboptimal=False)
                    listTriangles.append( t )
        display_progression_triangle(a, quiet, n, n, "first vertex", newline=True)
    return listTriangles
    

def all_triangles_in_perim(a, perim, suboptimal, base_first=False, quiet=True):
    if len(perim) <= 2:
        return []
    if suboptimal:
        return all_triangles_in_perim_suboptimal(a, perim, base_first=base_first)
    else:
        return all_triangles_in_perim_optimal(a, perim, quiet, base_first=base_first)

def triangles_homo_sorting_by_contents(a, perim, attempts, suboptimal, base_first=False):
    display_message(a, a.quiet, "Searching all triangles ... ", newline=True)
    candidates = all_triangles_in_perim(a, perim, suboptimal, quiet=a.quiet, base_first=base_first)
    display_message(a, a.quiet, "Finding triangles contents ... ", newline=True)
    fine = []
    n = len(candidates)
    tour = 0
    for triangle in candidates:
        display_progression_triangle(a, a.quiet, tour, n, "triangles")
        triangle.findContents() 
        if not suboptimal and len(triangle.contents) in [1093, 364, 121, 40, 13, 4, 1]:
             fine.append(triangle)
        tour += 1
    if not a.quiet:
        sys.stdout.write("\rFinding triangles contents : Fine={0}, Total={1} \n".\
                                format(len(fine),len(candidates)))
    display_message(a, a.quiet, "Triangles sorting by number of contents ... ")
    if not suboptimal:
        fine.sort(key=lambda l: len(l.contents), reverse=True)
        if not a.quiet:
            sys.stdout.write("\rTriangles sorting by number of contents : max(contents)={0}\n".\
                                format(len(fine[0].contents)))
            sys.stdout.flush()
        return fine
    candidates.sort(key=lambda l: len(l.contents), reverse=True)
    if not a.quiet:
        sys.stdout.write("\rTriangles sorting by number of contents : max(contents)={0}\n".\
                                format(len(candidates[0].contents)))
        sys.stdout.flush()
    return candidates

# the first triangle can be as you want (then keep_first=False)
# but next, we need to keep the given edge to extend the triangulation
# then keep_first must be True for all (but first triangle)
def triangles_in_perim_sorting_by_area(a, perim, suboptimal, base_first=False, quiet=True):
    if len(perim) <= 2:
        return []
    display_message(a, quiet, "Searching all triangles ... ", newline = True)
    candidates = all_triangles_in_perim(a, perim, suboptimal, base_first=base_first, quiet=quiet)
    display_message(a, quiet, "Computing all triangles area ... ", newline = True)
    tour = 0
    n = len(candidates)
    for triangle in candidates:
        display_progression_triangle(a, quiet, tour, n, "triangles")
        triangle.findArea()
        tour += 1
    display_progression_triangle(a, quiet, tour, n, "triangle", newline = True)
    display_message(a, quiet, "Sorting triangles by area ... ")
    candidates.sort(key=lambda t: t.area, reverse=True)
    display_message(a, quiet, "\rSorting triangles by area OK", newline = True)
    return candidates
    
def nearest_to_X_in_Y(a, X, Y):
    '''
        Return indexes of points in Y sorted by distances with X[i]
    '''
    
    if type(X) is list:
        n = len(X)
    else:
        n = 1
        X = [X]
    coords_X = np.array([a.node[p]['geo'] for p in X])
    coords_Y = np.array([a.node[p]['geo'] for p in Y])
    dists = [geometry.greatArcAng(coords_Y, coords_X[i] ) for i in range(n)]  
    closest = [np.argsort(dists[i]) for i in range(n)]
    return closest
    
class Polygon:
    def __init__(self, verts, a, convex=False):
        self.a = a
        self.verts = list(verts)
        self.pts = np.array([a.node[p]['xyz'] for p in verts])
        self.contents = []
        self.attempts = None
        self.triangulation = None
        self.multi_area = None
        
    def __eq__(self, other):
        if isinstance(other, Polygon):
            if len(self) != len(other):
                return False
            a = np.array(self.verts)
            for i in range(len(self)):
                a = np.roll(a, 1)
                if list(a) == list(other.verts):
                    return True
            # return all(x in self.verts for x in other.verts) and \
            #     all(x in other.verts for x in self.verts)
            return False
        return NotImplemented
    
    def __ne__(self, other):
        x = self.__eq__(other)
        if x is not NotImplemented:
            return not x
        return NotImplemented
    
    def __lt__(self, other):
        a = sorted(self.verts)
        b = sorted(other.verts)
        return a < b
        
    def __len__(self):
        return len(self.verts)
        
    def __repr__(self):
        return "{}".format(self.verts)
        
    def tostr(self):
        # Just a string representation of the triangle
        return str([self.a.node[i]['name'] for i in self.verts])
        
    def findContents(self,candidates=None):
        if candidates == None:
            candidates = xrange(self.a.order())
        for p in candidates:
            if p in self.verts:
                continue
            if geometry.spherePolygonContains(self.pts,self.a.node[p]['xyz']):
                self.contents.append(p)
        
    def split_on_subpoly(self, subpoly, convex=True, max=None):
        new_perims = geometry.split_perim(self.a, self.verts, subpoly.verts, max=max)
        for l in new_perims:
            for i in range(len(l)):
                # TODO : remove when len(l[i])<=2
                l[i] = Polygon(l[i], self.a, convex=convex)
        return new_perims
        

class Triangle(Polygon):
    def __init__(self,verts,a,exterior=False,suboptimal=_suboptimal):
        '''
        verts should be a 3-list of Portals
        verts[0] should be the final one used in linking
        exterior should be set to true if this triangle has no triangle parent
            the orientation of the outer edges of exterior Triangles do not matter
        '''
        # If this portal is exterior, the final vertex doesn't matter
        Polygon.__init__(self, verts, a)
        #self.verts = list(verts)
        #self.a = a
        self.area = None
        self.init_verts = verts
        self.exterior = exterior
        self.triangulation = [self]

        # This randomizes the Portal used for the jet link. I am
        # experimenting with having maxfield.triangulate and
        # Triangle.split choose this portal carefully, so don't
        # randomize
        """
        if exterior:
            # Randomizing should help prevent perimeter nodes from getting too many links
            final = np.random.randint(3)
            tmp = self.verts[final]
            self.verts[final] = self.verts[0]
            self.verts[0] = tmp
        """
        self.children = []
        self.center = None
        self.suboptimal = suboptimal
        #self.best_spliting = None
    

    def findArea(self):
        self.area = geometry.sphereTriArea( np.array([self.a.node[p]['geo'] for p in self.verts]) )
        
    def get_Total_Area(self):
        #if self.area == None:
        #    self.findArea()
        self.multi_area = self.area
        for i in range(len(self.children)):
            self.multi_area += self.children[i].get_Total_Area()
        return self.multi_area
    
    def findContents(self,candidates=None):
        if candidates == None:
            candidates = xrange(self.a.order())
        triangleKey = sum([1 << int(p) for p in self.verts])
        if triangleKey in triangleContentCache:
            self.contents.extend(triangleContentCache[triangleKey])
        else:
            for p in candidates:
                if p in self.verts:
                    continue
                if geometry.sphereTriContains(self.pts,self.a.node[p]['xyz']):
                    self.contents.append(p)
            triangleContentCache[triangleKey] = self.contents

    def near_points_from_verts(self):
        arg_nearVerts = nearest_to_X_in_Y(self.a, self.verts, self.contents)
        n = len(arg_nearVerts) # n = 3
        # keep only firsts :
        nearVerts = [self.contents[arg_nearVerts[i][0]] for i in range(n)]
        # contentLL = np.array([self.a.node[p]['geo'] for p in self.contents])
        # dists = [geometry.greatArcAng(contentLL, self.a.node[self.verts[i]]['geo'] ) for i in range(3)]        
        # arg_nearVerts = [np.argmin(dists[i]) for i in range(3)]
        # nearVerts = [self.contents[i] for i in arg_nearVerts]
        return nearVerts
        
    def near_points_from_edges(self): 
        contentPts = np.array([self.a.node[p]['xyz'] for p in self.contents])
        verts_xyz = [self.a.node[self.verts[i]]['xyz'] for i in range(3)]
        edges_xyz = verts_xyz - np.roll(verts_xyz, -1)
        norm_edges = np.sum(edges_xyz**2,1)
        vectContent = [contentPts-self.a.node[self.verts[i]]['xyz'] for i in range(3)]
        dists = [np.divide(np.sum(np.cross(vectContent[i], edges_xyz[i])**2,1),norm_edges[i]) for i in range(3)]
        nears = [self.contents[np.argmin(dists[i])] for i in range(3)]
        # AB = self.a.node[self.verts[1]]['xyz']- self.a.node[self.verts[0]]['xyz']
        # distAB = np.sum(np.cross(contentPts-self.a.node[self.verts[0]]['xyz'], AB)**2,1)
        # nearAB = self.contents[np.argmin(distAB)]
        # BC = self.a.node[self.verts[2]]['xyz']- self.a.node[self.verts[1]]['xyz']
        # distBC = np.sum(np.cross(contentPts-self.a.node[self.verts[1]]['xyz'], BC)**2,1)
        # nearBC = self.contents[np.argmin(distBC)]
        # AC = self.a.node[self.verts[2]]['xyz']- self.a.node[self.verts[0]]['xyz']
        # distAC = np.sum(np.cross(contentPts- self.a.node[self.verts[0]]['xyz'], AC)**2,1)
        # nearAC = self.contents[np.argmin(distAC)]
        # return [nearAB, nearBC, nearAC]
        return nears

    def randSplit(self):
        if len(self.contents) == 0:
            return 0
        p = self.contents[np.random.randint(len(self.contents))]        
        area = self.splitOn(p)
        for child in self.children:
            area += child.randSplit()
        return area

    def nearSplit(self):
        # Split on the node closest to final
        # Recursive function
        if len(self.contents) == 0:
            return 0
        closest = self.near_points_from_verts()
        area = self.splitOn(closest[0])
        for child in self.children:
            area += child.nearSplit()
        return area
    
    def fan_Split(self):
        # Split on the node closest to final
        # Recursive function
        if len(self.contents) == 0:
            return 0
        closest = self.near_points_from_edges()
        area = self.splitOn(closest[2])
        for child in self.children:
            area += child.fan_Split()
        return area
            
    def near_maxAREA_Split(self):
        # Recursive function
        if len(self.contents) == 0:
            return 0
        if len(self.contents) == 1:
            return self.splitOn(self.contents[0])
        # take the point closest to a vertex for which 
        # the big triangle are max area
        closest = self.near_points_from_verts()     
        areas = []
        n = len(closest)
        for i in range(n):
            t = Triangle([closest[i], self.verts[(i+1)%n], self.verts[(i+2)%n]], 
                        self.a, True, suboptimal=self.suboptimal)
            t.findArea()
            #pt1 = closest[i]
            #pt2 = self.a.node[self.verts[(i+1)%n]]['xyz']
            #pt3 = self.a.node[self.verts[(i+2)%n]]['xyz']
            #areas.append(geometry.planeTriArea([pt1, pt2, pt3]))
            areas.append(t.area)
        pref = np.argmax(areas)
        area = self.splitOn(closest[pref])
        for child in self.children:
            area += child.near_maxAREA_Split()
        return area
        
    def maxAREA_Split_optimize_on_children(self, optimized_polygon_done, p, tour, attempts, degree, quiet):
        if not quiet and not self.a.quiet:
            sys.stdout.write("\r[{0:20s}] {1:3}% ({2}/{3} attemps on degree {4})".\
                                format('='*(20*tour/attempts),100*tour/attempts,tour,attempts,degree))
            sys.stdout.flush()
        #print "maxAREA_Split_optimize_on_children :", self.init_verts, " with p=",p, " and tour =", tour
        self.verts = list(np.roll(self.init_verts, -(tour%3)))
        p_area = self.splitOn(p)
        for i in range(len(self.children)):
            child = self.children[i]
            if i == 0 : #only for opposite children
                j = bisect.bisect_left(optimized_polygon_done, child)
                if j != len(optimized_polygon_done) and child == optimized_polygon_done[j]:
                    optimized_child = optimized_polygon_done[j]
                    # print "Polygon already optimize :",optimized_child.verts, " by ", optimized_child.triangulation
                    # shift = optimized_child.verts.index(child.verts[0])
                    # print "shift = ", shift
                    # optimized_child.verts = list(np.roll(optimized_child.verts, -shift))
                    self.children[i] = optimized_child
                    #print "replacing children by ", self.children[i].verts
                    p_area += (self.children[i].multi_area - self.children[i].area)
                    continue
                #print "children not optimized :",child
            if not quiet and not self.a.quiet:
                print("")
            child.multi_area = child.maxAREA_Split_optimize(optimized_polygon_done, degree=degree+1, quiet=quiet)
            p_area += child.multi_area
            child.multi_area += child.area
            if i==0 and len(child.contents)>1:#only for opposite children
                #print "Optimized child added : ", child
                bisect.insort_left(optimized_polygon_done, child)
            if not quiet and not self.a.quiet:
                sys.stdout.write("\033[K\033[F")
                sys.stdout.flush()
        if not quiet and not self.a.quiet:
            sys.stdout.write("\033[K")
            sys.stdout.flush()
        return p_area
        
    def maxArea_best_candidate(self, optimized_polygon_done, candidates, degree, quiet):
        best_split = -1
        best_area = 0
        ordered_verts = self.verts
        tour = 0
        attempts = len(candidates)
        seen = set()
        for p in candidates:
            if p>0 and not p in seen:
                seen.add(p)
                p_area = self.maxAREA_Split_optimize_on_children(optimized_polygon_done, p, tour, attempts, degree, quiet)
                if p_area > best_area:
                    best_area = p_area
                    best_split = p
                    best_children = self.children
                    ordered_verts = self.verts
            tour += 1
        self.verts = ordered_verts
        return best_area, best_split, best_children
        
    def maxAREA_Split_optimize(self, optimized_polygon_done, degree=1, quiet=True):
        # Recursive function
        if len(self.contents) == 0:
            #if not quiet and not self.a.quiet:
            #    sys.stdout.write("\033[F")
            return 0
        if len(self.contents) == 1:
            #if not quiet and not self.a.quiet:
            #    sys.stdout.write("\033[F")
            return self.splitOn(self.contents[0])
        #print "triangle :", self.verts
        # candidates are the points nears a vertex or an edge of the triangle
        nearVerts = self.near_points_from_verts()
        if not self.exterior:
            #interdiction d'avoir le finalSplit comme base
            nearVerts[1] = 0
            nearVerts[2] = 0
        #print "nearVerts =", nearVerts, " in ", self.verts
        #nearEdges = self.near_points_from_edges()
        #print "nearEdges =", nearEdges, " in ", self.verts
        #nearEdges = list(np.roll(nearEdges, -1))
        candidates = nearVerts #+nearEdges
        best_area, best_split, best_children = self.maxArea_best_candidate(optimized_polygon_done, candidates, degree, quiet)
        self.children = best_children
        self.center = best_split
        return best_area
        
    def near_maxMULTI_Split(self):
        # Recursive function
        if len(self.contents) == 0:
            return 0
        if len(self.contents) == 1:
            return self.splitOn(self.contents[0])
        arg_nearVerts = nearest_to_X_in_Y(self.a, self.verts[0], self.contents)
        closest = self.contents[arg_nearVerts[0][0]]   
        area = self.splitOn(closest)
        area += self.children[0].near_maxMULTI_Split()
        return area
        
    def maxMULTI_Split_optimize_on_children(self, optimized_polygon_done, p, tour, attempts, degree, quiet):
        if not quiet and not self.a.quiet:
            sys.stdout.write("\r[{0:20s}] {1:3}% ({2}/{3} attemps on degree {4})".\
                                format('='*(20*tour/attempts),100*tour/attempts,tour,attempts,degree))
            sys.stdout.flush()
        p_area = self.splitOn(p)
        for i in range(1):
            child = self.children[i]
            j = bisect.bisect_left(optimized_polygon_done, child)
            if j != len(optimized_polygon_done) and child == optimized_polygon_done[j]:
                self.children[i] = optimized_polygon_done[j]
                #print "Polygon already optimize :",self.children[i].triangulation, self.children[i].multi_area
                p_area += (self.children[i].multi_area - self.children[i].area)
                continue
            #print "children not optimized :",child
            if not quiet and not self.a.quiet:
                print("")
            child.multi_area = child.maxMULTI_Split_optimize(optimized_polygon_done, degree=degree+1, quiet=quiet)
            p_area += child.multi_area
            child.multi_area += child.area
            if len(child.contents)>1:
                #print "Optimized child added : ", child
                bisect.insort_left(optimized_polygon_done, child)
            if not quiet and not self.a.quiet:
                sys.stdout.write("\033[K\033[F")
                sys.stdout.flush()
        if not quiet and not self.a.quiet:
            sys.stdout.write("\033[K")
            sys.stdout.flush()
        return p_area
        
    def maxMULTI_best_candidate(self, optimized_polygon_done,  candidates, degree, quiet):
        best_split = -1
        best_area = 0 
        tour = 0
        attempts = len(candidates)
        for p in candidates:
            p_area = self.maxMULTI_Split_optimize_on_children(optimized_polygon_done, p, tour, attempts, degree, quiet)
            if p_area > best_area:
                best_area = p_area
                best_split = p
                best_children = self.children
            tour += 1
        return best_area, best_split, best_children
        
    def maxMULTI_Split_optimize(self, optimized_polygon_done, degree=1, quiet=True):
        # Recursive function
        if len(self.contents) == 0:
            #if not quiet and not self.a.quiet:
            #    sys.stdout.write("\033[F")
            return 0
        if len(self.contents) == 1:
            #if not quiet and not self.a.quiet:
            #    sys.stdout.write("\033[F")
            return self.splitOn(self.contents[0])
        #print "triangle :", self.verts
        # candidates are the points nears a vertex or an edge of the triangle
        arg_nearVerts = nearest_to_X_in_Y(self.a, self.verts[0], self.contents)
        candidates = [self.contents[i] for i in arg_nearVerts[0]][:self.attempts]
        best_area, best_split, best_children = self.maxMULTI_best_candidate(optimized_polygon_done, candidates, degree, quiet)
        self.children = best_children
        self.center = best_split
        return best_area
        
        
    def splitOn_countlinks(self,p, b):
        opposite  =  Triangle([p,self.verts[1],
                               self.verts[2]],self.a,True)
        # The other two children must also use my final as their final
        adjacents = [\
                     Triangle([self.verts[0],\
                               self.verts[2],p],self.a),\
                     Triangle([self.verts[0],\
                               self.verts[1],p],self.a)\
                    ]
        
        self.children = [opposite]+adjacents
        self.center = p
        
        added_links = links_Triangle(opposite, b)
        added_links += links_Triangle(adjacents[0], b)
        added_links += links_Triangle(adjacents[1], b)

        for child in self.children:
            child.findContents(self.contents)
            child.findArea()
        return np.sum([child.area for child in self.children]), added_links
    
    def maxKEY_Split(self, b):
        # Recursive function
        if len(self.contents) == 0:
            return 0, []
        if len(self.contents) == 1:
            return self.splitOn_countlinks(self.contents[0], b)
        
        # First : looking to the vertex number of links
        # if keys all used,  split near of vertex with max links
        n = len(self.verts) # = 3
        verts_keys = np.array([self.a.node[self.verts[i]]['keys'] for i in range(n)])
        verts_degree = np.array([b.degree(self.verts[i]) for i in range(n)])
        max_linked = np.argmin(verts_keys-verts_degree)
        if verts_keys[max_linked]+0 < verts_degree[max_linked] :
            contentPts = np.array([self.a.node[p]['xyz'] for p in self.contents])
            displaces = contentPts - self.a.node[self.verts[max_linked]]['xyz']
            dists = np.sum(displaces**2,1)
            closest = np.argmin(dists)
            area, added_linked = self.splitOn_countlinks(self.contents[closest], b)
        else:
            # Look only on portals with max keys           
            contentKeys = np.array([self.a.node[p]['keys'] for p in self.contents])
            maxkeynumber = np.max(contentKeys)
            #print "maxkeynumber=",maxkeynumber
            contentwithMaxKeys = np.array([p for p in self.contents if self.a.node[p]['keys']==maxkeynumber])
            if len(contentwithMaxKeys)==1:
                area, added_linked = self.splitOn_countlinks(contentwithMaxKeys[0], b)
            else:
                # several portals have a maximum of keys
                # We must choose among the points that have the same number of keys
                # -> Split away from vertex with max key
                argVertMaxkey = np.argmax(verts_keys)
                candidates_vertex = self.verts[:]
                candidates_vertex.pop(argVertMaxkey) 
                n = len(candidates_vertex) # = 2
                
                contentPts = np.array([self.a.node[p]['xyz'] for p in contentwithMaxKeys])
                displaces = [contentPts - self.a.node[candidates_vertex[i]]['xyz'] for i in range(n)]
                dists = [np.sum(displaces[i]**2,1) for i in range(n)]
                closest = [np.argmin(dists[i]) for i in range(n)]
                
                keys_candidates_vertex = [self.a.node[candidates_vertex[i]]['keys'] for i in range(n)]
                if not np.all(keys_candidates_vertex == keys_candidates_vertex[0]):
                    # split closest to the vertex with min keys
                    pref = np.argmin(keys_candidates_vertex)
                else:
                    # take the point closest to a vertex for which 
                    # the big triangle are max area
                    areas = []
                    for i in range(n):
                        pt1 = contentPts[closest[i]]
                        pt2 = self.a.node[candidates_vertex[(i+1)%n]]['xyz']
                        pt3 = self.a.node[self.verts[argVertMaxkey]]['xyz']
                        areas.append(geometry.planeTriArea([pt1, pt2, pt3]))
                    #print "areas="+str(areas)
                    pref = np.argmax(areas)
                    
                area, added_linked = self.splitOn_countlinks(contentwithMaxKeys[closest[pref]], b)

        for child in self.children:
            child_area, child_added_linked = child.maxKEY_Split(b)
            area += child_area
            added_linked += child_added_linked
        return area, added_linked
        
    def split_homogenous(self, suboptimal, mindeg=1, maxdeg=8, quiet=True):
        #print "triangle : ", self.verts,"   mindeg = ",mindeg
        # Recursive function
        n = len(self.contents)
        if (maxdeg<=1) or n == 0:
            return 1
        if n == 1:
            self.splitOn(self.contents[0])
            return 2
        best_degree = 0
        best_split = None
        num = 0
        for p in self.contents:
            # if not quiet:
            #     sys.stdout.write("\r[{0:20s}] {1:2}% ({2}/{3} contents) Best degree : {4} ".\
            #                     format('='*(20*num/n),100*num/n,num,n,best_degree))
            #     sys.stdout.flush()
            self.splitOn(p)
            len_list = [len(self.children[i].contents) for i in range(3)]
            argchild = np.argsort(len_list)
            if not suboptimal and not np.all(len_list == len_list[0]):
                #print "not all ",len_list
                continue
            #else:
                #print "FOUND :",len_list
            d0 = self.children[argchild[0]].split_homogenous(suboptimal, mindeg=mindeg-1, maxdeg=maxdeg-1, quiet=quiet)
            if d0 < mindeg or d0 < best_degree:
                continue
            d1 = self.children[argchild[1]].split_homogenous(suboptimal, mindeg=mindeg-1, maxdeg=d0, quiet=quiet) # if use all points : mindeg=d0
            if d1 < mindeg:
                continue
            d2 = self.children[argchild[2]].split_homogenous(suboptimal, mindeg=mindeg-1, maxdeg=d0, quiet=quiet)
            if d2 < mindeg:
                continue
            degree = np.min([d0, d1, d2])
            if degree > best_degree:
                #print "best degree found :",degree
                best_degree = degree
                best_split = p
                best_children = self.children
            num += 1
        if best_degree > 0:
            self.children = best_children
            self.center = best_split
            return best_degree+1
        return 1
        
    def splitOn(self, p, first_vertex=0):
        # Splits this Triangle to produce 3 children using portal p
        # p is passed as the first vertex parameter in the
        # construction of 'opposite', so it will be opposite's
        # 'final vertex' unless randomization is used
        # 'opposite' is the child that does not share the final vertex
        # Because of the build order, it's safe for this triangle to
        # believe it is exterior
        #print "splitOn ", self.verts," on p=",p
        
        opposite  =  Triangle([p,self.verts[1],
                               self.verts[2]],self.a,True)
        # The other two children must also use my final as their final
        adjacents = [\
                     Triangle([self.verts[0],\
                               p, self.verts[2]],self.a),\
                     Triangle([self.verts[0],\
                               self.verts[1],p],self.a)\
                    ]
        
        self.children = [opposite]+adjacents
        #print "splitOn make children :", self.children
        self.center = p

        for child in self.children:
            child.findContents(self.contents)
            child.findArea()
        return np.sum([child.area for child in self.children])

    def buildFinal(self):
        #print "buildFinal for tri:",self.verts
        #print 'building final',self.tostr()
        if self.exterior:
            # Avoid making the final the link origin when possible
            #print self.verts,' is exterior'
            #print self.tostr(),'is exterior'
            if not (try_ordered_edge(self.a,self.verts[1],\
                               self.verts[0],self.exterior,self.suboptimal) \
                and try_ordered_edge(self.a,self.verts[2],\
                               self.verts[0],self.exterior,self.suboptimal) \
                ):
                return False
        else:
            #print self.verts,' is NOT exterior'
            #print self.tostr(),'is NOT exterior'
            if not (try_ordered_edge(self.a,self.verts[0],\
                               self.verts[1],self.exterior,self.suboptimal) \
                and try_ordered_edge(self.a,self.verts[0],\
                               self.verts[2],self.exterior,self.suboptimal) \
                ):
                return False

        if len(self.children) > 0:
            for i in [1,2]:
                if not self.children[i].buildFinal():
                    return False
        return True

    def buildExceptFinal(self):
        #print "buildExceptFinal for tri:",self.verts
        #print 'building EXCEPT final',self.tostr()
        if len(self.children) == 0:
            #print 'no children'
            p,q = self.verts[2] , self.verts[1]
            return try_ordered_edge(self.a,p,q,True,self.suboptimal)
            

        # Child 0 is guaranteed to be the one opposite final
        #print 'Child 0 of ',self.verts
        if not self.children[0].buildGraph():
            return False

        for child in self.children[1:3]:
            if not child.buildExceptFinal():
                return False
        return True

    def buildGraph(self):
        # print 'building',self.tostr()
        '''
        TODO
        A first generation triangle could have its final vertex's
        edges already completed by neighbors.
        This will cause the first generation to be completed when
        the opposite edge is added which complicates completing inside
        descendants.
        This could be solved by choosing a new final vertex (or
        carefully choosing the order of completion of first generation
        triangles).
        '''
        NoError = True
        #print "buildGraph for tri : ",self.verts
        if (                                                \
            self.a.has_edge(self.verts[0],self.verts[1]) or \
            self.a.has_edge(self.verts[1],self.verts[0])    \
           ) and                                            \
           (                                                \
            self.a.has_edge(self.verts[0],self.verts[2]) or \
            self.a.has_edge(self.verts[2],self.verts[0])    \
           # ) and                                            \
           # (                                                \
           #  self.a.has_edge(self.verts[1],self.verts[2]) or \
           #  self.a.has_edge(self.verts[2],self.verts[1])    \
           ):
            # print 'ERROR : Final vertex completed !!! ',self.verts
            # if (self.a.has_edge(self.verts[1],self.verts[2]) or \
            #     self.a.has_edge(self.verts[2],self.verts[1]) ):
            #     print "All vertex completed :'( "
            # else:
            #     print "Verts 1-2 not completed ??? "
            NoError = False
            #return False
        if not self.buildExceptFinal():
            print "buildExceptFinal Failed"
            return False
        #print "buildExceptFinal OK for tri : ",self.verts
        if not self.buildFinal():
            print "buildFinal Failed"
            return False
        return NoError

    def contains(self,pt):
        return np.all(np.sum(self.orths*(pt-self.pts),1) < 0)
        
    def display_edges_error(self, p, q):
        print 'a does NOT have edge',p,q
        print 'there is a programming error'
        print 'You should first call Build method ! '
        print 'a only has the edges:'
        for p,q in self.a.edges:
            print p,q
        print 'a has %s 1st gen triangles:'%len(self.a.triangulation)
        for t in self.a.triangulation:
            print t.verts
            
    def get_edges(self):
        edges = [(0,0)]*3
        for i in range(3):
            p = self.verts[i-1]
            q = self.verts[i-2]
            if not self.a.has_edge(p,q):
                p,q = q,p
            # The graph should have been completed by now, so the edge p,q exists
            edges[i] = (p,q)
            if not self.a.has_edge(p,q):
                self.display_edges_error(p, q)
        return edges
        
    def get_last_edge(self):
        edges = self.get_edges()
        edgeOrders = [self.a.edges[p,q]['order'] for p,q in edges]
        lastInd = np.argmax(edgeOrders)
        # The edge that completes this triangle
        p,q = edges[lastInd]
        return p, q, lastInd

    # Attach to each edge a list of fields that it completes
    def markEdgesWithFields(self, clean=False):
        if clean:
            cleanMarkedEdgesWithFields(self.a)
        edges = self.get_edges()
        p, q, lastInd = self.get_last_edge()
        self.a.edges[p,q]['fields'].append(self.verts)
        if not self.exterior:
            # the last edge depends on the other two
            del edges[lastInd]
            self.a.edges[p, q]['depends'].extend(edges)
        else:
            # in an exterior triangle that has children, only the edge
            # on the opposite side of the "final" vertex is a dependency;
            # childless exterior triangles can be built in any order
            if len(self.children) > 0:
                self.a.edges[p,q]['depends'].append(edges[0])
        for child in self.children:
            child.markEdgesWithFields()
        # all edges starting from inside this triangle have to be completed before it
        for c in self.contents:
            self.a.edges[p,q]['depends'].append(c)
        #print("edge %d-%d depends on: %s" % (p, q, self.a.edges[p,q]['depends']))

    def edgesByDepth(self,depth):
        # Return list of edges of triangles at given depth
        # 0 means edges of this very triangle
        # 1 means edges splitting this triangle
        # 2 means edges splitting this triangle's children 
        # etc.
        if depth == 0:
            return [ (self.verts[i],self.verts[i-1]) for i in range(3) ]
        if depth == 1:
            if self.center == None:
                return []
            return [ (self.verts[i],self.center) for i in range(3) ]
        return [e for child in self.children\
                  for e in child.edgesByDepth(depth-1)]
