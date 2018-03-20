#!/usr/bin/env python3

import argparse
import cairo
import re

parser = argparse.ArgumentParser(description='Input a fasta file and a motifs files --> \
Output a vector image of each gene with motif locations')
parser.add_argument("-f", "--filename", required=True, help='Fasta file containing \
single exons (capital) and flanking introns (lower-case)')
parser.add_argument("-m", "--motifs", required=True, help='A file containing one \
query motif per line (accepts degenerate symbols)')
args = parser.parse_args()
        

#Dictionaries to translate between regex and sequence
regex = {'A':'[Aa]','T':'[TtUu]','G':'[Gg]','C':'[Cc]','U':'[TtUu]',
        'W':'[AaTtUu]','S':'[GgCc]','M':'[AaCc]','K':'[GgTtUu]','R':'[AaGg]','Y':'[CcTtUu]',
        'B':'[CcGgTtUu]','D':'[AaGgTtUu]','H':'[AaCcTtUu]','V':'[AaCcGg]',
        'N':'[GgCcAaTtUu]'}

# Convert input motifs to regex
with open(args.motifs, 'r') as fh:
    regex_motifs=[] # create empty list to output converted motifs
    for line in fh:
        regex_motif=''
        motif=line.strip().upper()
        for each in motif:
            regex_motif+=regex[each]
        regex_motifs.append(regex_motif)

# Dict to store gene lengths and exon position
genes_index={}

# Dict to store motif positions
motifs_index={}

def main_coords(sequence):
    """Input single fasta sequence containing capital EXON region.
    Outputs list of [gene length, exon position, exon length]."""
    gene_len = len(seq)
    exon = re.search('[A-Z]',sequence)
    exon_start = exon.start()
    # Trim off upstream intron
    sequence = sequence[exon_start:]
    exon = re.search('[a-z]',sequence)
    exon_len = exon.start()
    out = [gene_len, exon_start, exon_len]
    return out

def motif_coords(sequence):
    """Input single fasta sequence containing capital EXON region.
    Save positions of desired motifs."""
    out=[]
    for motif in regex_motifs:
        mot = re.finditer(motif,sequence)
        for i in mot:
            out.append(i.start())
    return out

# Main loop
with open(args.filename,'r') as fasta:
    NR=0 #counts lines
    for line in fasta:
        NR+=1
        # Save first header line
        if NR == 1:
            header=line.strip()
            seq=''
        # Save sequence
        elif line.startswith('>') == False:
            # Combines sequence lines if there are new line characters
            seq=seq+line.strip()
        # Save important info below
        elif line.startswith('>'):
            genes_index[header[1:]] = main_coords(seq)
            motifs_index[header[1:]] = motif_coords(seq)
            # Save next header and empty seq
            header=line.strip()
            seq=''
            
    # Final fasta record saved below
    genes_index[header[1:]] = main_coords(seq)
    motifs_index[header[1:]] = motif_coords(seq)

# Save dictionary key containing max gene length
max_width = max(genes_index, key=genes_index.get)
# Save max gene length
max_width = genes_index[max_width][0]

# Set width based on max gene length and height based on gene count
width,height = int(max_width*1.25), int(len(genes_index)*120)

surface = cairo.SVGSurface("plot.svg", width, height)
context = cairo.Context(surface)

# Useful coordinates
x_start = width*0.1  # x position to start drawing each gene
y_div = height/(len(genes_index)+1) # coefficient for y position of each gene

# First draw genes as lines with solid rectangles for exons
gene_cnt=1
for val in genes_index.values():
    gene_len=val[0] # length of gene
    ex_start=val[1] # starting pos of exon
    ex_end=ex_start+val[2] # ending pos of exon
    
    #draw gene line
    context.set_line_width(2)
    context.move_to(x_start, y_div*gene_cnt)
    context.line_to(x_start+gene_len, y_div*gene_cnt)
    context.stroke()
    
    #draw exon box
    context.set_line_width(10)
    context.move_to(x_start+ex_start, y_div*gene_cnt)
    context.line_to(x_start+ex_end, y_div*gene_cnt)
    context.stroke()
    # increment index for next gene
    gene_cnt+=1


# Add lines indicating locations of motifs along genes
gene_cnt=1
for m_list in motifs_index.values():
    for pos in m_list:
        context.set_line_width(18)
        context.set_source_rgb(250,0,0)
        context.move_to(x_start+pos, y_div*gene_cnt)
        context.line_to(x_start+pos+2, y_div*gene_cnt)
        context.stroke()
    gene_cnt+=1

gene_cnt=1
for name in genes_index:
    #text label above gene
    context.move_to(10, (y_div*gene_cnt)-20)
    context.set_source_rgb(0,0,0)
    context.select_font_face("Open Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    context.set_font_size(16)
    context.show_text(name)
    gene_cnt+=1

# Add helpful info at the bottom
context.move_to(10,height-15)
context.set_font_size(12)
context.show_text("*Black boxes represent exons. Red lines \
indicate locations of input motifs.")

surface.finish()
