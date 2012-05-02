#!/usr/bin/python

import Image
import ImageDraw
import ImageFont
import PSDraw
import numpy
import copy
import sys
import os
from math import ceil
from operator import itemgetter
from collections import defaultdict

######   CHANGE THESE    ##########
SORT_ORDER='freq'

NON_JN_COVERAGE_KEY = "Non Junction"
JUNCT_COLORS=( ( 255, 0, 0 ), ( 255, 165, 0 ), ( 255, 255, 0 ), \
                  ( 0, 128, 0 ), ( 0, 0, 255 ), ( 128, 0, 128 ) )

junction_color_map = { NON_JN_COVERAGE_KEY: (0,0,0) }
def find_junction_color( junct_type, junct_index=None ):
    if junction_color_map.has_key( junct_type ):
        return junction_color_map[ junct_type ]
    
    if junct_type == NON_JN_COVERAGE_KEY: 
        return ( 0, 0, 0 )
    
    if junct_index == None:
        return ( 0, 0, 0 )
    
    color_index = junct_index % len( JUNCT_COLORS )
    junction_color_map[ junct_type ] = JUNCT_COLORS[ color_index ]
    return JUNCT_COLORS[ color_index ]


WIDTH=1800
HEIGHT=200
HEIGHT_GROWTH_VALUE = 200

UPPER_BORDER=100
LOWER_BORDER=20
LEFT_BORDER=20
RIGHT_BORDER=20

CDS_HEIGHT = 12
EXON_HEIGHT = 12
EXON_COLOR='black'
CONNECTED_EXON_COLOR='grey'
EXON_SPACING = 2

HISTOGRAM_SIZE=80

NUM_BINS = 100

###### DONT CHANGE THESE ###########
VERTICAL_SPACE = HEIGHT - UPPER_BORDER - LOWER_BORDER
HORIZONTAL_SPACE = WIDTH - LEFT_BORDER - RIGHT_BORDER

font_base_dir = os.path.join( os.path.split( os.path.abspath( __file__ ) )[0], './fonts/' )
arial_45 = ImageFont.truetype(os.path.join( font_base_dir, "arial.ttf" ), 45)
arial_bd_10 = ImageFont.truetype( os.path.join( font_base_dir, "arialbd.ttf"), 10)
arial_bd_12 = ImageFont.truetype( os.path.join( font_base_dir, "arialbd.ttf"), 12)

###### END DONT CHANGE THESE #######

class GenePlot( object ):
    def _draw_header( self ):
        font = arial_45
        title = "Gene: " + self.gene.name + " (%s strand)" % self.gene.strand
        self.draw.text( (LEFT_BORDER,20), title, fill='black', font=font )

    resize_time = 0
    def resize( self, new_size ):
        """Resize the image ( self.im ) to have new_size.

        """
        new_im = Image.new( "RGB", new_size, (255, 255, 255, 0 ) )
        new_im.paste( self.im, (0,0) )
        self.im = new_im
        self.resize_time += 1
        self.draw = ImageDraw.Draw( self.im )
    
    def _get_curr_upper_edge( self ):
        return self._curr_upper_edge

    def _set_curr_upper_edge( self, value ):
        if value + 100 > self.im.size[1]:
            new_height = self.im.size[1] + HEIGHT_GROWTH_VALUE
            self.resize( (WIDTH, new_height) )
        self._curr_upper_edge = value
       
    curr_upper_edge = property( _get_curr_upper_edge, _set_curr_upper_edge )
    
    def __init__( self, gene, connected_exons, reads=None ):
        self.gene = gene
        self.reads = reads
        self.connected_exons = connected_exons
        
        self.left_gene_bndry = gene.boundaries.start
        self.right_gene_bndry = gene.boundaries.stop + 1
        self.gene_len = self.right_gene_bndry - self.left_gene_bndry

        self.exon_bndrys = self.get_corrected_exon_boundaries( \
            self.gene.exon_bndrys )
        self.exon_percs = self.convert_exon_pos_to_percentages()
        
        self.exon_font = arial_bd_10
        
        self.im = Image.new( "RGB", ( WIDTH, HEIGHT ), (255, 255, 255, 0 ) )
        self.draw = ImageDraw.Draw( self.im )
        
        self.curr_upper_edge = UPPER_BORDER
        
        self._draw_header()
        
        if self.reads != None:
            self.plot_read_coverage( self.reads )
        
        self.draw_exons()
                
        self.draw_connected_exons()
        
    def save( self, output_prefix="" ):
        output_fname = output_prefix + self.gene.name + ".jpg"
        self.im.save( output_fname, "JPEG")
    
    def show( self ):
        self.im.show()

    def get_corrected_exon_boundaries( self, gene_exon_bndrys ):
        """Convert exons to correct format for use in the rest of plotting

        """
        return numpy.array([ (start, stop + 1) for start, stop in gene_exon_bndrys ])
    
    def convert_exon_pos_to_percentages( self ):
        """Convert exon positions, in genomic coordinates,
        to percentages of the gene. To aid in plotting.
        
        """
        assert self.exon_bndrys.shape[1] == 2
        return (numpy.array(self.exon_bndrys, dtype=float) - self.left_gene_bndry)/self.gene_len

    @staticmethod
    def get_exon_coords( ( start, stop ), upper_offset ):
        # the upper left corner is just the upper border
        # plus the level times the exon height + 1
        upper_edge = upper_offset
        lower_edge = upper_edge + EXON_HEIGHT
        left_edge = LEFT_BORDER + int(start*HORIZONTAL_SPACE)
        right_edge = LEFT_BORDER + int(stop*HORIZONTAL_SPACE)
        return ( left_edge, upper_edge, right_edge, lower_edge )
    
    def draw_exons( self ):
        """Draw ( possibly overlapping ) exons. 
        
        Returns the lower coordinate of the space used ( all horizontal space is used )
        """
        exon_percs = self.exon_percs

        # we need to be careful because exons may overlap
        undrawn_exons = list( enumerate(exon_percs ) )
        # the following makes sure that we always put exons that overlap
        # into different levels, while not using more vertical space than 
        # we need to. The algorithm is pretty simple: we initialize the level
        # to 1, then we draw as many exons as we can at that level, keeping
        # track of the ones we can't draw. Then, we increment the level, and
        # draw as many as we can. We repeat this until there is nothing else 
        # to be drawn.
        while len( undrawn_exons ) > 0:
            prev_drawn_exon_coords = None
            exons_to_draw = copy.deepcopy( undrawn_exons )
            del undrawn_exons[:]
            # draw all of the undrawn exons that we can
            for index, coords in exons_to_draw:
                # if this doesn't overlap the previously drawn exon, draw it
                if prev_drawn_exon_coords == None \
                        or coords[0] > prev_drawn_exon_coords[1]:
                    exon_coords = self.get_exon_coords(coords, self.curr_upper_edge )
                    self.draw.rectangle( exon_coords, fill=EXON_COLOR )
                    self.draw.text( ( exon_coords[0]+2, exon_coords[1] ), \
                               str(index), fill='white', font=self.exon_font )
                    prev_drawn_exon_coords = coords
                # otherwise, save it for later
                else:
                    undrawn_exons.append( ( index, coords ) )

            # move onto the next level
            self.curr_upper_edge += ( EXON_HEIGHT + EXON_SPACING )
        
        return

    def draw_connected_exons( self ):
        """Draw ( possibly overlapping ) connected exons. 
        
        Returns the lower coordinate of the space used ( all horizontal space is used )
        """
        def draw_connected_exon( coords, color=CONNECTED_EXON_COLOR ):
            upper_edge = self.curr_upper_edge
            exon1_coords = self.get_exon_coords(coords[:2], upper_edge)
            self.draw.rectangle( exon1_coords, fill=color )
            # draw.text( ( exon1_coords[0]+2, exon2_coords[1] ), 
            #     str(index), fill='white', font=self.exon_font )
            exon2_coords = self.get_exon_coords(coords[2:], upper_edge)
            self.draw.rectangle( exon2_coords, fill=color )
            # draw.text( ( exon1_coords[0]+2, exon2_coords[1] ), 
            #     str(index), fill='white', font=self.exon_font )
            top = exon1_coords[1] + 4
            bot = exon1_coords[3] - 4
            coords = ( exon1_coords[2], top, exon2_coords[0], bot )
            self.draw.rectangle( coords, fill=color )
        
        exon_percs = self.exon_percs

        # draw header
        font = arial_bd_12
        self.draw.text( (LEFT_BORDER,self.curr_upper_edge+6), \
                        "Connected Exons:", fill='black', font=font )
        self.curr_upper_edge += 20
        
        # build connected exon coords
        connected_exon_coords = []
        for start_exon, stop_exon in self.connected_exons:
            junct_type = ( self.exon_bndrys[ start_exon ][1],
                              self.exon_bndrys[ stop_exon ][0] )
            color = find_junction_color( junct_type )
            connected_exon_coords.append( \
                ( ( exon_percs[start_exon][0], exon_percs[start_exon][1], \
                    exon_percs[stop_exon][0], exon_percs[stop_exon][1] ), \
                  color ) )

        connected_exon_coords.sort()
        connected_exon_percs = numpy.array( [ i[0] for i in connected_exon_coords ] )
        connected_exon_colors = [ i[1] for i in connected_exon_coords ]
        # we need to be careful because exons may overlap
        curr_level = 0
        undrawn_exons = list( enumerate(zip( connected_exon_percs, connected_exon_colors )) )
        # the following makes sure that we always put exons that overlap
        # into different levels, while not using more vertical space than 
        # we need to. The algorithm is pretty simple: we initialize the level
        # to 1, then we draw as many exons as we can at that level, keeping
        # track of the ones we can't draw. Then, we increment the level, and
        # draw as many as we can. We repeat this until there is nothing else 
        # to be drawn.
        while len( undrawn_exons ) > 0:
            prev_drawn_exon_coords = None
            exons_to_draw = copy.deepcopy( undrawn_exons )
            del undrawn_exons[:]
            # draw all of the undrawn exons that we can
            for index, ( coords, color ) in exons_to_draw:
                # if this doesn't overlap the previously drawn exon, draw it
                if prev_drawn_exon_coords == None \
                    or coords[0] > prev_drawn_exon_coords[-1]:
                    draw_connected_exon( coords, color )
                    prev_drawn_exon_coords = coords
                # otherwise, save it for later
                else:
                    undrawn_exons.append( ( index, ( coords, color ) ) )
            # move onto the next level
            curr_level += 1
            self.curr_upper_edge += ( EXON_HEIGHT + EXON_SPACING )
        
        return

    def draw_transcript_header( self, md ):
        font = arial_bd_12
        if SORT_ORDER == 'lambda':
            text = "Lambda: %e      Freq:   %e        LogPrb:   %e" \
                % ( md.lasso_lambda, md.freq, md.improbability )
        elif SORT_ORDER == 'freq':
            text = "Freq:   %e      Lambda: %e        LogPrb:   %e" \
                % ( md.freq, md.lasso_lambda, md.improbability )
        
        self.draw.text( (LEFT_BORDER, self.curr_upper_edge+6), \
                        text, fill='black', font=font )
        
        self.curr_upper_edge += 20
    
    def draw_transcript( self, exon_indices, color=EXON_COLOR ):
        exon_percs = self.exon_percs
        
        for exon_index in exon_indices:
            exon_coords = self.get_exon_coords( \
                exon_percs[ exon_index ], self.curr_upper_edge )
            self.draw.rectangle( exon_coords, fill=color )
        
        self.curr_upper_edge += EXON_HEIGHT + EXON_SPACING
        
        return
    
    def draw_transcripts( self, transcripts, min_lambda=0.0, min_improb=0.05 ):
        prev_value = None
        transcripts.sort(SORT_ORDER)
        for exon_indices, md in transcripts.iter_transcripts_and_metadata():
            if SORT_ORDER == 'lambda': value = md.lasso_lambda
            elif SORT_ORDER == 'freq': value = md.freq
            else: assert False

            #if value != prev_value:
            self.draw_transcript_header( md )
            
            if md.improbability < min_improb: color='purple'
            elif  md.lasso_lambda < min_lambda: color = 'red'
            elif md.freq > 0: color = 'green'
            else: color = 'black'
            
            self.draw_transcript( exon_indices, color )
            
            if SORT_ORDER == 'lambda': prev_value = md.lasso_lambda
            elif SORT_ORDER == 'freq': prev_value = md.freq
            else: assert False
        
    def plot_read_coverage( self, reads ):
        def get_read_coverage( reads ):
            # find the starts and stops for verification of junction reads
            starts = set()
            stops = set()
            for start, stop in self.exon_bndrys:
                starts.add( start )
                stops.add( stop )
            
            def gene_len_zeros_array_factory():
                return numpy.zeros( self.gene_len )

            read_coverages = defaultdict(gene_len_zeros_array_factory)
            
            def find_covered_bases_and_jn_type( read ):
                # find the bases that this read covers.
                pos = read.pos
                junct_type = []
                base_indices_to_add_coverage = []
                # make sure that cigar is a sequence of matching(exon) and 
                # skipping(intron) regions and then store junctions as tuple 
                # of exon_stop and exon_start positions also add positions 
                # within each exon at which to add coverage for this junction type
                for index, ( region_type, region_size ) in enumerate( read.cigar ):
                    # make sure that the first region_type is not an intron 
                    assert not( region_type == 3 and index == 0 )

                    # If this is a matching region
                    if region_type == 0:
                        # check that pos+matching region is an exon stop position 
                        # for internal matches. This ensures that the junction 
                        # corresponds to an exon, except at the end of 
                        # the read.
                        if index < len( read.cigar ) - 1 \
                                and pos + region_size + 1 not in stops: return ([], None)

                        # add the base positions that we want to plot to 
                        # base_indices_to_add_coverage.
                        bases = range( pos - ( self.left_gene_bndry-1 ), \
                                       pos + region_size + 1 \
                                       - (self.left_gene_bndry-1) - 1 )
                        
                        base_indices_to_add_coverage.extend( bases )
                        
                        # pos shifts to an exon stop
                        pos += region_size + 1

                    # if this is a skipped region
                    elif region_type == 3:
                        pos += region_size - 1
                        # make sure this correpsonds with an intron
                        if pos + 1 not in starts: return ([], None)
                    # this is neither a skipped region nor a match, 
                    # so we ignore this read.
                    # STUB for dealing with other region types
                    else:
                        return ([], None)

                    # add position to junction type unless it is the end of the read
                    if index < len( read.cigar ) - 1:
                        if region_type == 3:
                            junct_type.append( pos + 1 )
                        elif region_type == 0:
                            junct_type.append( pos )
                        else:
                            assert False

                return base_indices_to_add_coverage, tuple( junct_type )

            def update_read_coverage( read ):
                if len( read.cigar ) == 1:
                    start = read.pos - self.left_gene_bndry + 1
                    stop = start + read.alen
                    read_coverages[NON_JN_COVERAGE_KEY][ start:stop ] += 1

                elif len( read.cigar ) > 1:
                    # check that read contains odd number of cigar regions possibly 
                    # corresponding to alternating exon-intron... regions starting 
                    # and ending with exon regions
                    if len( read.cigar ) % 2 == 0:
                        return 
                    covered_bases, junct_type = find_covered_bases_and_jn_type( read )
                    # if we dont know the junction or this is an invalid junction read 
                    if junct_type == None: return
                    
                    # add the bases
                    for index in covered_bases:
                        read_coverages[ junct_type ][ index ] += 1
                return
            
            for read1, read2 in reads.iter_paired_reads( self.gene.boundaries ):
                update_read_coverage( read1 )
                update_read_coverage( read2 )
            
            return read_coverages

        def get_binned_reads( all_read_type_basepair_coverages ):
            def find_fractional_cnt_for_pixel( basepair_coverages, lower_bp, upper_bp ):
                frac_cnt = 0
                # if this pixel is entirely within a basepair, then add the 
                # fraction of the base that that pixel occupies. This is exactly 
                # like a standard histogram: we want the sum over the base 
                # coverage over pixels to add up to the total base coverage
                if int( upper_bp ) == int( lower_bp ):
                    frac_cnt += (upper_bp - lower_bp) \
                                * basepair_coverages[ int( lower_bp ) ]
                else:
                    # get the lower fraction
                    frac_cnt +=  ( ceil( lower_bp ) - lower_bp ) \
                                 * basepair_coverages[ int( lower_bp ) ]
                    # get the middle basepairs
                    for cnt in basepair_coverages[ ceil( lower_bp ) : int( upper_bp ) ]:
                        frac_cnt += cnt
                    # get the upper fraction
                    if upper_bp - 1e-6 > int(upper_bp):
                        frac_cnt +=  ( upper_bp - int(upper_bp) ) \
                                     * basepair_coverages[ int( upper_bp ) ]

                return frac_cnt

            def num_pixels_zeros_array_factory():
                return numpy.zeros( num_pixels )

            # create binned reads arrays for non-junct and each junct type
            pixels_coverages = defaultdict( num_pixels_zeros_array_factory )
            max_bin_height = 0
            # loop through each pixel, and find the coverage taking into account
            # fractional basepair coverage
            for pixel_index in xrange(num_pixels):
                lower_bp = pixel_index * bps_per_pixel
                upper_bp = lower_bp + bps_per_pixel
                
                total_pixel_frac_cnt = 0
                frac_cnts = {}
                # build the fractional count arrays
                for read_type, basepair_coverages in \
                        all_read_type_basepair_coverages.iteritems():
                    frac_cnt = find_fractional_cnt_for_pixel( \
                        basepair_coverages, lower_bp, upper_bp )
                    frac_cnts[ read_type ] = frac_cnt
                    total_pixel_frac_cnt += frac_cnt
                    
                # find the maximum height over all pixels
                if total_pixel_frac_cnt > max_bin_height:
                    max_bin_height = total_pixel_frac_cnt
                
                for read_type, frac_cnt in frac_cnts.iteritems():
                    pixels_coverages[ read_type ][ pixel_index ] += frac_cnt
            
            # normalize the binned reads
            for read_type, pixels_coverage in pixels_coverages.iteritems():
                pixels_coverage /= max_bin_height
            
            return pixels_coverages
        
        def plot_binned_reads( pixels_coverages ):
            boundaries = set()
            for start, stop in self.exon_bndrys:
                boundaries.add( start - self.left_gene_bndry )
                boundaries.add( stop - self.left_gene_bndry )
            boundaries = sorted( boundaries )
            boundaries = [ int(bndry / bps_per_pixel) for bndry in boundaries  ]

            # TODO color bins to match clustered average heights
            # draw average bin coverage
            for start, stop in zip( boundaries[:-1], boundaries[1:] ):
                if stop - start == 0: continue
                avg_frac = pixels_coverages[NON_JN_COVERAGE_KEY][ start : stop ].mean()
                if avg_frac == 0: continue
                left_coords = ( \
                    LEFT_BORDER + start, \
                        self.curr_upper_edge + (1 - avg_frac) * HISTOGRAM_SIZE) 
                right_coords = ( \
                    LEFT_BORDER + stop, \
                        self.curr_upper_edge + HISTOGRAM_SIZE )
                self.draw.rectangle( ( left_coords, right_coords), \
                                         fill=(128, 128, 128, 128) )
            
            # loop through pixel indices and plot each read_type at this pixel
            for pixel_index in xrange( num_pixels ):
                curr_hist_height = 0
                # draw read coverage for non-junct first and then junction reads on top
                for read_index, ( read_type, pixels_coverage ) in \
                        enumerate( sorted( pixels_coverages.iteritems() ) ):
                    junct_frac = pixels_coverage[ pixel_index ]
                    # if there is nothing to plot, skip this read type
                    if junct_frac == 0: continue
                    
                    # stack coverage on top of previously plotted coverages
                    left_coords = ( \
                        LEFT_BORDER + pixel_index, \
                            self.curr_upper_edge + ((1 - junct_frac) * \
                            HISTOGRAM_SIZE) - curr_hist_height )
                    right_coords = ( \
                        LEFT_BORDER + pixel_index, \
                            self.curr_upper_edge + HISTOGRAM_SIZE - curr_hist_height ) 
                    # choose color so that adjacent (ordered by junction coordinates)
                    # junctions have different colors
                    color = find_junction_color( read_type, read_index )
                    self.draw.rectangle( ( left_coords, right_coords), fill=color )
                    curr_hist_height += junct_frac * HISTOGRAM_SIZE
            
            # draw baoundary markers
            for bndry in boundaries:
                left_coords = ( LEFT_BORDER + bndry, self.curr_upper_edge ) 
                right_coords = ( LEFT_BORDER + bndry, self.curr_upper_edge + \
                                     HISTOGRAM_SIZE )
                self.draw.rectangle( ( left_coords, right_coords), \
                                         fill=( 128, 128, 255, 128) )
            return
        
        
        # get the read starts histogram
        read_coverages = get_read_coverage( reads )

        # convert basepairs into pixels for use in next two functions
        num_pixels  = HORIZONTAL_SPACE
        bps_per_pixel = float( self.gene_len )/num_pixels                
            
        pixels_coverages = get_binned_reads( read_coverages )
        plot_binned_reads( pixels_coverages  )
        #self.im.paste( im, (LEFT_BORDER, self.curr_upper_edge) )
        self.curr_upper_edge += HISTOGRAM_SIZE + EXON_SPACING
        
        return
    
if __name__ == '__main__':
    sys.path.insert(0, os.path.abspath(".") )
    from reads import Reads
    from gene_models import GeneModel, parse_junctions_file
    
    #from slide import Reads, GeneModel, parse_junctions_file
    from slide import build_objects
    genes, all_transcripts = build_objects( open( sys.argv[1] ) )
    for gene_name in genes:
        gene = genes[ gene_name ]
        transcripts = all_transcripts[ gene_name ]
        connected_exons = transcripts.iter_connected_exons()

        # if we included reads so that we can plot read density, parse
        # the file.
        if len( sys.argv ) == 3: reads = Reads( sys.argv[2], "rb" )
        else: reads = None
        
        plot = GenePlot( gene, connected_exons, reads )
        plot.draw_transcripts( transcripts, 0, 0 )
        plot.save()
        plot.show()
    
