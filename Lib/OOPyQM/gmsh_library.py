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
        self.points = [point1, point2]
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
    plane_surfaces:  list of instances of Plane_Surface().
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
class Physical_Line():
    """"
    Representation of a physical line in a surface.
    Parameters:
    id_tag (int)   Current state of the counter in lines, loops or surfaces.
    lines:       List of instances of Line() Class.
    """
    def __init__(self, id_tag, lines):
        self.id_tag = id_tag + 1
        self.lines = lines
        self.tags = ''
        for ln in lines:
            if ln == lines[-1]: #is the last
                self.tags += str(ln.id_tag)+ ' '
            else:
                self.tags += str(ln.id_tag)+ ', '
    def __str__(self):
        return 'Physical Line(%s) = {%s};\n' %(self.id_tag, self.tags)        
        
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
        self.kind = 'empty'
        self.id_tag = id_tag + 1
        self.origin = deepcopy(origin)
        self.width = width
        self.height = height
        self.transfinite = transfinite
        self.phys_lines = []
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
        self.perimeters = [self.perimeter]
        self.surface = Plane_Surface(self.perimeter.id_tag, self.perimeter, transfinite = self.transfinite)
        self.surfaces = [self.surface]
        self.last_tag = self.surface.id_tag
        self.phys_tag = phys_tag
        if self.phys_tag:
            self.physical_surface = Physical_Surface(self.surface.id_tag, [self.surface])
            self.physical_surfaces = [self.physical_surface]
            self.last_tag = self.physical_surface.id_tag
        self.recombine = recombine   
        from copy import deepcopy
        self.last_point = deepcopy(self.points[-1].id_tag)
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
        p.append( Point(0,circle.center.x + r, circle.center.y))
        p.append( Point(0,circle.center.x - r, circle.center.y))
        p.append( Point(0,circle.center.x, circle.center.y + r))
        p.append( Point(0,circle.center.x, circle.center.y - r))
        for point in p:
            if not self.contains(point): return False
        return True
        
class Unitary_cell(Rectangle):
    """
    Represents a Rectangle shape to be used as unit cell
    Does not have rotations implemented.

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
    def add_circular_inclussion(self, id_tag, loc, radius, phys_tag = None):
        """
        loc:      Is a location inside the unitary cell that is given 
        (Point()) relative to the corner. This will mark the place where the
                  inclussion is going to be placed
        phys_tag: Will be the physical entity tag for the inclussion.
        """
        self.kind = 'circle'
        self.id_tag = id_tag + 1
        #assert isinstance(loc, Point)
        loc.x += self.origin.x
        loc.y += self.origin.y
        loc.id_tag = self.points[-1].id_tag + 1
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
        self.last_tag = self.surfaces[-1].id_tag
        if self.phys_tag:
            self.physical_surfaces = []
            if self.circle.phys_tag:
                self.physical_surfaces.append(self.circle.physical_surface)
            self.physical_surfaces.append( Physical_Surface(self.phys_tag, self.surfaces[5:]))
            self.last_tag = self.physical_surfaces[-1].id_tag
        from copy import deepcopy
        self.last_point = deepcopy(self.points[-1].id_tag)
    
    def __str__(self):
        if len(self.points) > 5:
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
        else:
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
class Grid():
    """
    A grid is a class that is meant to represent a group of unitary cells
    in a semicrystal. It can be initiated by giving a matrix of 
    integer values each of them representing a kind of unitary cell, be it a
    homogeneus cell or a cell with inclussions.
    
    Parameters:
    ----------
    
    sketch:   2D Array of integer values. Each of which represents a kind of
              unitary cell  
    """
    
    def __init__(self, origin, a, sketch, dicts, phys_tags = [None, None], transfinite = False):
        from copy import deepcopy        
        self.sketch = sketch
        self.dicts = dicts
        self.transfinite = transfinite
        self.sk_shape = [len(sketch),len(sketch[0])]
        self.cells = deepcopy(self.sketch)
        assert self.check_rectangular(), "Grid is not rectangular"
        self.origin = origin
        self.phys_tags = phys_tags
        self.phys_lines = []
        ri = 0 # Row id counter
        for row in self.cells:
            ci = 0 # Column id counter
            for cell in row:
                p = deepcopy(dicts[cell]) # Properties
                #--------- What to do for the first unitary cell------------
                if ri == 0 and ci == 0:                    
                    cell = Unitary_cell(0, origin, a, a, transfinite = self.transfinite)
                    if p.kind == 'circle':
                        cc = deepcopy(p.origin) # circle center
                        cell.add_circular_inclussion(1, cc, p.r)
                    elif p.kind == 'square':
                        raise NotImplementedError('Pending implementation')
                    elif p.kind == 'empty':
                        pass
                    else:
                        raise NotImplementedError("No more kinds of cells yet")
                #-- What to do for the first row excluding the first cell---
                elif ri == 0 and ci != 0:
                    last_tag = row[ci-1].last_tag + 1
                    last_point = deepcopy(row[ci-1].last_point)+1
                    st_point = deepcopy(row[ci-1].points[1]) #Starting point
                    st_point.id_tag = last_point
                    cell = Unitary_cell(last_tag, st_point, a, a, transfinite = self.transfinite)
                    cell.lines[3] = row[ci-1].lines[1]
                    # Insert inclusion
                    if p.kind == 'circle':
                        cc = deepcopy(p.origin) # circle center
                        cc.tag_id = cell.last_tag
                        cell.add_circular_inclussion(cell.last_tag, cc, p.r)
                    elif p.kind == 'square':
                        raise NotImplementedError('Pending implementation')
                    elif p.kind == 'empty':
                        cell.perimeters[-1].lines[3] = cell.lines[3].id_tag
                        pass
                    else:
                        raise NotImplementedError("No more kinds of cells yet")
                    # Redefine shared points in line definitions
                    for line in cell.lines:
                        if cell.points[0] in line.points:
                            line.points[line.points.index(cell.points[0])] = row[ci-1].points[1]
                        elif cell.points[3] in line.points:
                            line.points[line.points.index(cell.points[3])] = row[ci-1].points[2]
                    # Replaces common points
                    cell.points[0] = row[ci-1].points[1]
                    cell.points[3] = row[ci-1].points[2]
                    # This fixes the direction of the shared side
                    if p.kind == 'empty':
                        cell.perimeters[-1].lines[3] = - cell.perimeters[-1].lines[3]
                    else:
                        cell.perimeters[-1].lines[2] = - cell.perimeters[-1].lines[2]
                    
                #------------ What to do for the rest of cells--------------
                else:
                    if ci == 0:
                        last_tag = self.cells[ri-1][-1].last_tag + 1
                        last_point = deepcopy(self.cells[ri-1][-1].last_point)+1
                        st_point = deepcopy(self.cells[ri-1][ci].points[3]) #Starting point
                    else: 
                        last_tag = row[ci-1].last_tag + 1
                        last_point = deepcopy(row[ci-1].last_point)+1
                        st_point = deepcopy(row[ci-1].points[1]) #Starting point
                    st_point.id_tag = last_point
                    cell = Unitary_cell(last_tag, st_point, a, a, transfinite = self.transfinite)
                    # Replace vertical shared line, if the cell is not on the 
                    # left boundary                    
                    if ci != 0:
                        cell.lines[3] = row[ci-1].lines[1]  
                    # Replace shared bottom line with the top line of the 
                    # row bellow
                    cell.lines[0] = self.cells[ri-1][ci].lines[2]
                    # Insert inclusion
                    if p.kind == 'circle':
                        cc = deepcopy(p.origin) # circle center
                        cc.tag_id = cell.last_tag
                        cell.add_circular_inclussion(cell.last_tag, cc, p.r)
                    elif p.kind == 'square':
                        raise NotImplementedError('Pending implementation')
                    elif p.kind == 'empty':
                        cell.perimeters[-1].lines[3] = cell.lines[3].id_tag
                        cell.perimeters[-1].lines[0] = cell.lines[0].id_tag
                        pass
                    else:
                        raise NotImplementedError("No more kinds of cells yet")
                    # Redefine shared points in line definitions
                    if ci != 0:
                        for line in cell.lines:
                            if cell.points[0] in line.points:
                                line.points[line.points.index(cell.points[0])] = row[ci-1].points[1]
                            elif cell.points[3] in line.points:
                                line.points[line.points.index(cell.points[3])] = row[ci-1].points[2]
                            elif cell.points[1] in line.points:
                                line.points[line.points.index(cell.points[1])] = self.cells[ri-1][ci].points[2]
                        # Replaces common points
                        cell.points[0] = row[ci-1].points[1]
                        cell.points[3] = row[ci-1].points[2]
                        cell.points[1] = self.cells[ri-1][ci].points[2]
                        # This fixes the direction of the shared sides
                        if p.kind == 'empty':
                            cell.perimeters[-1].lines[3] = - cell.perimeters[-1].lines[3]
                            cell.perimeters[-1].lines[0] = - cell.perimeters[-1].lines[0]
                        else:
                            cell.perimeters[-1].lines[2] = - cell.perimeters[-1].lines[2]
                            cell.perimeters[-4].lines[2] = - cell.perimeters[-4].lines[2]
                    else:
                        for line in cell.lines:
                             if cell.points[0] in line.points:
                                line.points[line.points.index(cell.points[0])] = self.cells[ri-1][ci].points[3]
                             elif cell.points[1] in line.points:
                                line.points[line.points.index(cell.points[1])] = self.cells[ri-1][ci].points[2]
                        # Replaces common points
                        cell.points[0] = self.cells[ri-1][ci].points[3]
                        cell.points[1] = self.cells[ri-1][ci].points[2]
                        # This fixes the direction of the shared sides
                        if p.kind == 'empty':
                            cell.perimeters[-1].lines[0] = - cell.perimeters[-1].lines[0]
                        else:
                            cell.perimeters[-4].lines[2] = - cell.perimeters[-4].lines[2]
                        
                row[ci] = cell
                
                ci += 1
            ri += 1 
    def check_rectangular(self):
        """
        Check if the grid is rectangular.
        Returns:
            True, if it is.
        """
        sizes = []
        i = 0
        for row in self.cells:
            if i > 0:
                if len(row) not in sizes:
                    return False
            sizes.append(len(row))
            i += 1
        return True
    def extract_bottom_lines(self, tag = None):
        """
        This method returns a set or list of all lines that conform the bottom.
        
        """
        lines = []
        for cell in self.cells[0]:
            lines.append(cell.lines[0])
        if tag == None:
            phys_bot_lines = Physical_Line(self.cells[-1][-1].last_tag, lines)
        else:
            phys_bot_lines = Physical_Line(tag, lines)
        self.cells[-1][-1].last_tag += 1 
        self.phys_lines.append(phys_bot_lines)
    def extract_top_lines(self, tag = None):
        """
        This method returns a set or list of all lines that conform the top.
        
        """
        lines = []
        for cell in self.cells[-1]:
            lines.append(cell.lines[2])
        if tag == None:
            phys_top_lines = Physical_Line(self.cells[-1][-1].last_tag, lines)
        else:
            phys_top_lines = Physical_Line(tag, lines)
        self.cells[-1][-1].last_tag += 1
        self.phys_lines.append(phys_top_lines)
    def extract_right_lines(self, tag = None):
        """
        This method returns a set or list of all lines that conform the
        right side.
        """
        lines = []
        for row in self.cells:
            lines.append(row[-1].lines[1])
        if tag == None:
            phys_right_lines = Physical_Line(self.cells[-1][-1].last_tag, lines)
        else:
            phys_right_lines = Physical_Line(tag, lines)
        self.cells[-1][-1].last_tag += 1
        self.phys_lines.append(phys_right_lines)
    def extract_left_lines(self, tag = None):
        """
        This method returns a set or list of all lines that conform the
        right side.
        """
        lines = []
        for row in self.cells:
            lines.append(row[0].lines[3])
        if tag == None:
            phys_left_lines = Physical_Line(self.cells[-1][-1].last_tag, lines)
        else:
            phys_left_lines = Physical_Line(tag, lines)
        self.cells[-1][-1].last_tag += 1
        self.phys_lines.append(phys_left_lines)
    def define_boundary(self, tags = [None,None,None,None]):
        """
        This method extracts lines in a clockwise manner.
        If no tags are given it will ad tags from the last tag defined.
        """
        assert len(tags) == 4, "tags must be a 4 element list"
        self.extract_bottom_lines(tag = tags[0])
        self.extract_right_lines(tag = tags[1])
        self.extract_top_lines(tag = tags[2])
        self.extract_left_lines(tag = tags[3])
    
    def __str__(self):      
        grid_str = ''
        points = ''
        lines = ''
        perimeters = ''
        surfaces = ''
        phys_tags =''
        recombine =''
        phys_lines = ''
        rectangle_surfaces = []
        inclussion_surfaces = []
        repeated = dict()      
        for row in self.cells:
            for cell in row:
                for point in cell.points:
                    if point.__str__() not in repeated:
                        repeated[point.__str__()] = 1
                        points +=  point.__str__()
                    else:
                        repeated[point.__str__()] += 1
                for line in cell.lines:
                    if line.__str__() not in repeated:
                        repeated[line.__str__()] = 1
                        lines += line.__str__()
                    else:
                        repeated[line.__str__()] += 1
                for perimeter in cell.perimeters:
                    if perimeter.__str__() not in repeated:
                        repeated[perimeter.__str__()] = 1
                        perimeters += perimeter.__str__()
                    else:
                        repeated[perimeter.__str__()] += 1
                for surface in cell.surfaces:
                    if surface.__str__() not in repeated:
                        repeated[surface.__str__()] = 1
                        surfaces +=  surface.__str__()
                    else:
                        repeated[surface.__str__()] += 1
                if cell.kind == 'circle':
                    rectangle_surfaces.extend(cell.surfaces[5:])
                    inclussion_surfaces.extend(cell.surfaces[:5])
                elif cell.kind == 'square':
                    raise NotImplementedError('No squares yet')
                elif cell.kind == 'empty':
                    rectangle_surfaces.extend(cell.surfaces)
                else:
                    raise NotImplementedError('No more inclusion types available')
            for surface in cell.surfaces:
                if cell.recombine:
                    reco = 'Recombine Surface {%s};\n'%(surface.id_tag)
                    recombine += reco
        if self.phys_lines != []:
            for pl in self.phys_lines:
                phys_lines += pl.__str__()
            
        phys_tags +=  Physical_Surface(self.phys_tags[0], rectangle_surfaces).__str__()
        phys_tags +=  Physical_Surface(self.phys_tags[1], inclussion_surfaces).__str__()
        grid_str += points + lines + perimeters + surfaces + phys_tags + recombine + phys_lines
        return grid_str
class Properties_bag():
    """
    This class is suposed to contain properties as to make the construction 
    of Grids easier

    Paramenters:
    kind
        - 'circle'
        - 'square'        
        - 'empty'
    origin  Instance of class Point
    a       Width of rectangular inclusion
    b       Height of rectangular inclusion
    r       Radius of circular inclusion
    """
    def __init__(self, kind, origin = None, a = None, b = None, r = None):
        self.kind = kind
        self.origin = origin
        self.a = a
        self.b = b
        self.r = r
    def __str__(self):
        if self.kind == 'circle':
            return " Circle with center at %s has radius %s." %(self.origin, self.r)
        elif self.kind == 'rectangle':
            return " Rectangle with corner at %s has width %s and height %s." %(self.origin, self.a,self.b)
        elif self.kind == 'empty':
            return " Empty cell"
        else:
            raise NotImplementedError("No more inclussion shapes allowed.")
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        