import numpy as np
import re
import pysvg
import pysvg.parser
from itertools import chain

MAX_LINE_WIDTH = 30.0 # Set the maximal line width for stroke 
FLUX_CUTOFF = 20.0

class Markers(object):
    
    def __init__(self):
        pass

    def getXML(self):
        return \
            """
            <marker id="end_triangle" viewBox="0 0 12 12" refX="0" refY="6" markerUnits="strokeWidth" markerWidth="3" markerHeight="3" orient="auto">
            <path d="M 0 0 L 12 6 L 0 12 z" fill="black" stroke="black"/>
            </marker>
            <marker id="start_triangle" viewBox="0 0 12 12" refX="12" refY="6" markerUnits="strokeWidth" markerWidth="3" markerHeight="3" orient="auto">
            <path d="M 12 0 L 0 6 L 12 12 z" fill="black" stroke="black"/>
            </marker>
            """

class DrawFlux(object):

    def __init__(self, svg_fname):
        # parse the input SVG
        self.svg        = pysvg.parser.parse(svg_fname)
        
        main_group = self.svg.getAllElements()[-2] # the main G structure is one before last
        defs_elem  = main_group.getElementsByType(pysvg.structure.Defs)[0]
        defs_elem.insertElementAt(Markers(), 1)
        
        layer_rxn  = main_group.getElementByID(u'Layer_rxn')[0] # get the reaction layer
        reactions  = layer_rxn.getElementsByType(pysvg.structure.G) # get all the reaction G structures
        
        # make a dictionary from reaction names to the SVG elements of their arrows
        self.rxn_svg_dict = {}
        self.rxn_marker_start_dict = {}
        self.rxn_marker_end_dict = {}

        for r in reactions:
            rxn_id = DrawFlux.normalize_rxn_name(r.get_id())
            self.rxn_svg_dict.setdefault(rxn_id, []).append(r)
            self.rxn_marker_start_dict.setdefault(rxn_id, [])
            self.rxn_marker_end_dict.setdefault(rxn_id, [])

            # gather all the arrowheads that are originally set in the SVG
            # since only they should be changed. The paths also contain other lines
            # that should not have any arrowheads.
            linking = r.getElementsByType(pysvg.linking.A)[0]
            for path in linking.getElementsByType(pysvg.shape.Path):
                    
                if path.get_marker_start() != None:
                    self.rxn_marker_start_dict[rxn_id].append(path)
                if path.get_marker_end() != None:
                    self.rxn_marker_end_dict[rxn_id].append(path)
        
    @staticmethod
    def normalize_rxn_name(s):
        # strip all non-alphanumerics characters from the gene names in the model
        s = re.sub('[^\w]', '', s.lower())
        s = s.replace('_', '')
        return s
    
    def Reset(self):
        for r in chain.from_iterable(self.rxn_svg_dict.itervalues()):
            r.set_stroke('black')
            r.set_stroke_opacity(0.2)
            r.set_stroke_width(2)
        
        for p in chain.from_iterable(self.rxn_marker_start_dict.itervalues()):
            p.set_marker_start(None)

        for p in chain.from_iterable(self.rxn_marker_end_dict.itervalues()):
            p.set_marker_end(None)

    def ToSVG(self, model, solution, html):
        self.Reset()

        model_rxn_names = [DrawFlux.normalize_rxn_name(r.name)
                           for r in model.reactions]
        
        fluxes = np.array(solution.x, dtype=float, ndmin=1)
        direction = np.sign(fluxes)
        fluxes = np.abs(fluxes)
        
        above_cutoff = np.nonzero(fluxes > FLUX_CUTOFF)[0]
        fluxes[above_cutoff] = 0
        
        widths = 1 + (fluxes / np.max(fluxes)) * MAX_LINE_WIDTH
        
        # find the model reactions in the SVG according to their name and change the
        # color and width of the arrow according to the solution flux
        for rxn_name in set(model_rxn_names):
            i = model_rxn_names.index(rxn_name) # names are not unique, always choose the first index
            if rxn_name not in self.rxn_svg_dict:
                #print 'Cannot find %s in the SVG file (%d)' % (rxn_name, width[i])
                continue

            for r in self.rxn_svg_dict[rxn_name]:
                if i in above_cutoff:
                    r.set_stroke('red')
                    r.set_stroke_opacity(0.5)
                    r.set_stroke_width(10)
                elif direction[i] != 0:
                    r.set_stroke('green')
                    r.set_stroke_opacity(1)
                    r.set_stroke_width(widths[i])
                else:
                    r.set_stroke('black')
                    r.set_stroke_opacity(0.2)
                    r.set_stroke_width(2)
                
                if direction[i] == -1:
                    for p in self.rxn_marker_start_dict[rxn_name]:
                        p.set_marker_start('url(#start_triangle)')
                elif direction[i] == 1:
                    for p in self.rxn_marker_end_dict[rxn_name]:
                        p.set_marker_end('url(#end_triangle)')
                    
        html.write(self.svg.getXML())
