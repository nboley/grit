import sys, os
from collections import defaultdict

def meta_data( string ):
    return dict( x.split("=") for x in string.split(";" ) )

def main():
    # make a dict of CDS transcripts
    with open( sys.argv[1] ) as fp:
        cds_transcripts = set()
        mRNAs = set()
        transcripts = defaultdict( list )
        for line in fp:
            if line.startswith( "#" ): continue
            data = line.split()
            try:
                if data[2] == 'mRNA':
                    mRNAs.add( meta_data( data[8] )['ID'] )
                elif data[2] == 'protein':
                    cds_transcripts.add( meta_data( data[8] )['Derives_from'] )
                elif data[2] in ('CDS',  'exon'):
                    parents = meta_data( data[8] )['Parent'].split( "," )
                    for parent in parents:
                        data[-1] = parent
                        transcripts[parent].append( "\t".join( data ) )
            except IndexError:
                continue
        
        for transcript_id, exons in transcripts.iteritems():
            if transcript_id not in mRNAs:
                continue
            strands = [ exon.split()[6] for exon in exons ]
            if len( set( strands ) ) > 1: continue
            #if transcript_id in cds_transcripts:
            #    continue
            for line in exons:
                if line.startswith( "dmel_mitochondrion_genome" ):
                    continue
                print "chr" + line
main()
