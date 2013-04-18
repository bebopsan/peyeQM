# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 09:30:35 2013

Library for definition of elements to use as gmsh input

@author: santiago
"""


class Point(object):
    """Represents a point in 2-D space."""
    def __init__(self, id_tag, x=0, y = 0):
        id_tag = id_tag+1
        self.x = x
        self.y = y
        self.id_tag = id_tag
    def __str__(self):
        return 'Point(%d) = {%f, %f, 0, 1.0};\n' %(self.id_tag, self.x, self.y)
    def __add__(self, point):
        id_tag = point.id_tag 
        new_point = Point(id_tag, self.x + point.x, self.y + point.y)
        return new_point
    def distance(self, other):
        """Computes the distance from this location to other."""
        import math
        dx = self.x - other.x
        dy = self.y - other.y
        return math.sqrt(dx*dx + dy*dy)
    def isLeft(self, other): return self.x < other.x
    def isRight(self, other): return self.x > other.x
    def isAbove(self, other): return self.y > other.y
    def isBelow(self, other): return self.y < other.y
class Line(object):
    """
    Represents a line in a 2D space. Takes point instances as arguments.
    Parameters:
    -----------
    id_tag (int):   Current state of the counter in lines, loops or surfaces.
    point1,point2 (Point()):  Instances of Point class, represent the initial
                             and final coordinates of a straight line
    phys_tag(int): pysical entity to which the line blongs
    transfinite:   Tells if the line is going to be transfinite 
    (int, float)   Default value is False, but the way to define 
                   a transfinite line is to give a set where the first value 
                   is the number of divisions and the second the nonlinear
                   parameter explain in gmsh documentation.
                                
    """
    def __init__(self, id_tag, point1, point2, phys_tag = None, transfinite = False):
        self.id_tag = id_tag + 1
        self.phys_tag = phys_tag
        self.points = (point1, point2)
        self.lenght = point1.distance(point2)
        self.transfinite = transfinite
    def __str__(self):
        line_text ='Line(%d) = {%d, %d};\n'\
                %(self.id_tag, self.points[0].id_tag, self.points[1].id_tag)
        if self.transfinite == False:
            return line_text
        else:
            tran_text = 'Transfinite Line {%d} = %d Using Progression %s;\n'\
                        %(self.id_tag, self.transfinite[0], self.transfinite[1])
            return line_text + tran_text
class Circle_arc():
    """
    Defines circle arcs to be used in the construction of circular shapes in gmsh.
    Parameters:
    -----------
    id_tag (int):  Current state of the counter in lines, loops or surfaces.
    points = [point1, center, point2 :  
                    list of Instances of Point class, represent the 
                    initial point, the center of the circle and
                    and final point.
    phys_tag(int): pysical entity to which the line blongs
    transfinite:   Tells if the line is going to be transfinite 
    (int, float)   Default value is False, but the way to define 
                   a transfinite line is to give a set where the first value 
                   is the number of divisions and the second the nonlinear
                   parameter explain in gmsh documentation.
    """
    def __init__(self, id_tag, points, phys_tag = None, transfinite = False):
        self.id_tag = id_tag + 1
        self.phys_tag = phys_tag
        self.points = points
        self.transfinite = transfinite
    def __str__(self):
        arc_text ='Circle(%d) = {%d, %d, %d};\n'\
                %(self.id_tag, self.points[0].id_tag,self.points[1].id_tag,\
                  self.points[2].id_tag)
        if self.transfinite == False:
            return arc_text
        else:
            tran_text = 'Transfinite Line {%d} = %d Using Progression %s;\n'\
                        %(self.id_tag, self.transfinite[0], self.transfinite[1])
            return arc_text + tran_text
class Shape(object):
    """The parent class that all shapes inherit from."""
class Circle(Shape):
    def __init__(self, id_tag, center, radius, ratio = 2.5, phys_tag = None,\
                transfinite = False, recombine = True):
        self.center = center
        self.radius = radius
        self.id_tag = id_tag
        self.transfinite = transfinite
        square_hside = radius/ratio 
        dcorner = square_hside
        square_corner = Point(center.id_tag, center.x - dcorner, center.y - dcorner)
        self.square = Rectangle(self.id_tag, square_corner, 2*square_hside,\
                                2*square_hside, phys_tag = phys_tag, \
                                transfinite = self.transfinite)
                #------------ Assemble point instances of the circle ------------
        self.points = []
        self.points.append(self.center)
        self.points.extend(self.square.points)
        dp = radius * 0.70710678118654757# distance to point component
        self.points.append(Point(self.points[-1].id_tag, center.x - dp, center.y - dp))
        self.points.append(Point(self.points[-1].id_tag, center.x + dp, center.y - dp))
        self.points.append(Point(self.points[-1].id_tag, center.x + dp, center.y + dp))
        self.points.append(Point(self.points[-1].id_tag, center.x - dp, center.y + dp))
        #------------ Assemble line instances of the rectangle ------------
        self.lines = []
        self.lines.extend(self.square.lines)
        arc_1 = Circle_arc(self.lines[-1].id_tag, [self.points[-4], self.center,\
                      self.points[-3]], transfinite = self.transfinite)
        self.lines.append(arc_1)
        arc_2 = Circle_arc(self.lines[-1].id_tag, [self.points[-3], self.center,\
                      self.points[-2]], transfinite = self.transfinite)
        self.lines.append(arc_2)
        arc_3 = Circle_arc(self.lines[-1].id_tag, [self.points[-2], self.center,\
                      self.points[-1]], transfinite = self.transfinite)
        self.lines.append(arc_3)
        arc_4 = Circle_arc(self.lines[-1].id_tag, [self.points[-1], self.center,\
                      self.points[-4]], transfinite = self.transfinite)
        self.lines.append(arc_4)
        for i in range(4):
            diag = Line(self.lines[-1].id_tag, self.square.points[i],\
                         self.points[-4+i], transfinite = self.transfinite)
            self.lines.append(diag)                     
        #---------- Define line loops and surfaces ---------------------
        lids = [] # line ids
        for line in self.lines:
            lids.append(line.id_tag)
        self.perimeters = []
        self.perimeters.append(self.square.perimeter)
        self.surfaces = []
        self.surfaces.append(self.square.surface)
        for i in range(4):
            if i == 3:
                self.perimeters.append(Line_loop(self.perimeters[-1].id_tag, \
                            [-lids[i],lids[i+ 8], lids[i+4], -lids[i+5]]))
            else:
                self.perimeters.append(Line_loop(self.perimeters[-1].id_tag, \
                               [-lids[i],lids[i+ 8], lids[i+4], -lids[i+ 9]]))
            
            self.surfaces.append( Plane_Surface(self.perimeters[-1].id_tag,\
                          self.perimeters[-1], transfinite = self.transfinite))
        self.phys_tag = phys_tag
        if self.phys_tag:
            self.physical_surface = Physical_Surface(self.phys_tag, self.surfaces)
        self.recombine = recombine          
    def __str__(self):
        circle_string = ''
        for point in self.points:
            circle_string = circle_string + point.__str__()
        for line in self.lines:
            circle_string = circle_string + line.__str__()
        for perimeter in self.perimeters:
            circle_string = circle_string + perimeter.__str__()
        for surface in self.surfaces:
            circle_string = circle_string + surface.__str__()
        if self.phys_tag:
            circle_string = circle_string + self.physical_surface.__str__()
            if self.recombine:
                reco = 'Recombine Surface {%s};\n'%(surface.id_tag)
                circle_string = circle_string + reco
        return circle_string
    def contains(self, loc):
        """Returns True if Location loc is inside this circle
        (including the boundary)."""       
        return loc.distance(self.center) <= self.radius

class Line_loop():
    """
    This is used to define a surface. A line loop contains line tags
    which are signed if the direction of the line is against the direction
    of the loop
    
    Parameters:
    ----------
    id_tag (int): Current state of the counter in lines, loops or surfaces.
    lines (list): list of the signed id_tags of each line that conforms 
                  a loop
    """
    def __init__(self, id_tag, lines):
        self.id_tag = id_tag + 1
        self.lines = lines
    def __str__(self):
        str_list = ''
        for line in self.lines:
            if line == self.lines[-1]:
                str_list = str_list + str(line)
            else:
                str_list = str_list + str(line) + ', '
        return 'Line Loop(%s) = {%s};\n'%(self.id_tag, str_list)
class Plane_Surface():
    """
    A surface is the way to tell gmsh what to mesh when doing 2D meshes.
    Parameters:
    ----------
    id_tag (int): Current state of the counter in lines, loops or surfaces.
    line_loop
    (LineLoop()): An instance of Line_loop() class. it marks the closed 
                  boundary of one surface
    """
    def __init__(self, id_tag, line_loop, transfinite = False):
        self.id_tag = id_tag + 1
        self.line_loop = line_loop
        self.transfinite = transfinite
    def __str__(self):
        if not self.transfinite:
            return 'Plane Surface(%s) = {%s};\n' %(self.id_tag, self.line_loop.id_tag)
        else:
            surf = 'Plane Surface(%s) = {%s};\n' %(self.id_tag, self.line_loop.id_tag)
            transf = 'Transfinite Surface {%s};\n'%(self.id_tag)
            return surf + transf
class Physical_Surface():
    """
    Physical print for a surface
    Parameters:
    ----------
    id_tag (int): Current state of the counter in lines, loops or surfaces.
    plane_surface: instance of Plane_Surface().
    """        
    def __init__(self, id_tag, plane_surfaces):
        self.id_tag = id_tag
        self.plane_surface = plane_surfaces
        self.tags = ''
        for ps in plane_surfaces:
            if ps == plane_surfaces[-1]:
                self.tags += str(ps.id_tag)+ ' '
            else:
                self.tags += str(ps.id_tag)+ ', '
    def __str__(self):
        return 'Physical Surface(%s) = {%s};\n' %(self.id_tag, self.tags)        
        
class Rectangle(Shape):
    """
    Represents a Rectangle shape. Does not have rotations implemented.
    Parameters:
    ----------
    id_tag (int):   Current state of the counter in lines, loops or surfaces.
    Origin:     Instance of Point class, represent the bottom left corner
    width:           
    height:   can be integer or float. 
    phys_tag(int): pysical entity to which the line blongs
    transfinite:   Tells if the lines are going to be transfinite 
    (int, float)   Default value is False. see transfinite parameter in Line 
                    class documentation.
    """
    def __init__(self, id_tag, origin, width, height, phys_tag= None,\
                transfinite = False, recombine = True):
        from copy import deepcopy
        self.id_tag = id_tag + 1
        self.origin = deepcopy(origin)
        self.width = width
        self.height = height
        self.transfinite = transfinite
        #------------ Assemble point instances of the rectangle ------------
        self.points = []
        self.points.append(origin)
        self.points.append(Point(self.points[0].id_tag, origin.x+ width,origin.y))
        self.points.append(Point(self.points[1].id_tag, origin.x+ width,\
                                origin.y + height))
        self.points.append(Point(self.points[2].id_tag, origin.x,\
                                origin.y + height))
        self.center = Point(self.points[3].id_tag, origin.x + self.width/2.0, origin.y + self.height/2.0,)                    
        #------------ Assemble line instances of the rectangle ------------
        self.lines = []
        line_1 =  Line(self.id_tag, self.points[0],self.points[1],\
                      transfinite = self.transfinite)
        self.lines.append(line_1)
        line_2 =  Line(line_1.id_tag, self.points[1],self.points[2],\
                      transfinite = self.transfinite)
        self.lines.append(line_2)
        line_3 =  Line(line_2.id_tag, self.points[2],self.points[3],\
                      transfinite = self.transfinite)
        self.lines.append(line_3)
        line_4 =  Line(line_3.id_tag, self.points[3],self.points[0],\
                      transfinite = self.transfinite)
        self.lines.append(line_4)
        #---------- Define line loops and surface ---------------------
        line_ids = []
        for line in self.lines:
            line_ids.append(line.id_tag)
        self.perimeter = Line_loop(self.lines[-1].id_tag, line_ids)
        self.surface = Plane_Surface(self.perimeter.id_tag, self.perimeter, transfinite = self.transfinite)
        self.last_tag = self.surface.id_tag
        self.phys_tag = phys_tag
        if self.phys_tag:
            self.physical_surface = Physical_Surface(self.surface.id_tag, [self.surface])
            self.last_tag = self.physical_surface.id_tag
        self.recombine = recombine      
    def __str__(self):
        rectangle_string = ''
        for point in self.points:
            rectangle_string = rectangle_string + point.__str__()
        for line in self.lines:
            rectangle_string = rectangle_string + line.__str__()
        rectangle_string = rectangle_string + self.perimeter.__str__()
        rectangle_string = rectangle_string + self.surface.__str__()
        if self.phys_tag:
            rectangle_string = rectangle_string + self.physical_surface.__str__()
        if self.recombine:
            reco = 'Recombine Surface {%s};\n'%(self.surface.id_tag)
            rectangle_string = rectangle_string + reco
        return rectangle_string
        
    def contains(self, loc):
        """return True if Location loc is inside this rectangle
        (including the boundary)"""
        
        if loc.isLeft(self.points[3]): return False
        if loc.isRight(self.points[1]): return False
        if loc.isAbove(self.points[3]): return False
        if loc.isBelow(self.points[1]): return False
        return True
    def contains_circ(self, circle):
        assert isinstance(circle,Circle), "Print this is not a circle"
        p = []
        r = circle.radius
        p.append( Point(0,circle.center.x + r))
        p.append( Point(0,circle.center.x - r))
        p.append( Point(0,circle.center.y + r))
        p.append( Point(0,circle.center.y - r))
        for point in p:
            if not self.contains(point): return False
        return True
        
class Unitary_cell(Rectangle):
    def add_circular_inclussion(self, id_tag, loc, radius, phys_tag = None):
        """
        loc:      Is a location inside the unitary cell that is given 
        (Point()) relative to the corner. This will mark the place where the
                  inclussion is going to be placed
        phys_tag: Will be the physical entity tag for the inclussion.
        """
        #assert isinstance(loc, Point)
        print len(self.points)
        loc.x += self.origin.x/2.
        loc.y += self.origin.x/2.
        loc.id_tag = self.points[-1].id_tag + 1
        print loc
        self.circle = Circle(self.last_tag, loc, radius, phys_tag = phys_tag, transfinite = self.transfinite)
        assert self.contains_circ(self.circle), "Circle is not inside rectangular"\
        " shape"
        self.points.extend(self.circle.points)
        self.lines.extend(self.circle.lines)
        for i in range(4):
            diag = Line(self.lines[-1].id_tag, self.circle.points[-4+i],\
                         self.points[i], transfinite = self.transfinite)
            self.lines.append(diag)      
        #---------- Define line loops and surfaces ---------------------
        self.perimeters = []
        self.perimeters.extend(self.circle.perimeters)
        self.surfaces = []
        self.surfaces.extend(self.circle.surfaces)
                
        lids = [] # line ids
        for line in self.lines:
            lids.append(line.id_tag)
            
        for i in range(4):
           if i == 3:
               self.perimeters.append(Line_loop(self.perimeters[-1].id_tag, \
                           [-lids[i+ 8],lids[i+ 16], lids[i], -lids[i+ 13]]))
           else:
               self.perimeters.append(Line_loop(self.perimeters[-1].id_tag, \
                              [-lids[i+ 8],lids[i+ 16], lids[i], -lids[i+ 17]]))
           self.surfaces.append( Plane_Surface(self.perimeters[-1].id_tag,\
                          self.perimeters[-1], transfinite = self.transfinite))
        if self.phys_tag:
            self.physical_surfaces = []
            if self.circle.phys_tag:
                self.physical_surfaces.extend(self.circle.physical_surfaces)
            self.physical_surfaces.append( Physical_Surface(self.phys_tag, self.surfaces[5:]))
    def __str__(self):
        unit_cell_string = ''
        for point in self.points:
            unit_cell_string = unit_cell_string + point.__str__()
        for line in self.lines:
            unit_cell_string = unit_cell_string + line.__str__()
        for perimeter in self.perimeters:
            unit_cell_string = unit_cell_string + perimeter.__str__()
        for surface in self.surfaces:
            unit_cell_string = unit_cell_string + surface.__str__()
        if self.phys_tag:
            for phys_surf in self.physical_surfaces:
                unit_cell_string = unit_cell_string + phys_surf.__str__()
        for surface in self.surfaces:
            if self.recombine:
                reco = 'Recombine Surface {%s};\n'%(surface.id_tag)
                unit_cell_string = unit_cell_string + reco
        return unit_cell_string
