
import os
import re
import sys
import pysam
import getopt

__authors__  = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Apply changes in an assembly file to a fasta'


PYTHONVERSION = sys.version_info[:2]
if PYTHONVERSION < (3,0):
    range = xrange
    import string
    NUCLEOTIDE_COMPLEMENT = \
        string.maketrans(
            'acgturyswkmbdhvnACGTURYSWKMBDHVN',
            'tgcaayrwsmkvhdbnTGCAAYRWSMKVHDBN'
        )
else:
    NUCLEOTIDE_COMPLEMENT = {
        65: 'T',  97: 't',
        66: 'V',  98: 'v',
        67: 'G',  99: 'g',
        68: 'H', 100: 'h',
        71: 'C', 103: 'c',
        72: 'D', 104: 'd',
        75: 'M', 107: 'm',
        77: 'K', 109: 'k',
        78: 'N', 110: 'n',
        82: 'Y', 114: 'y',
        83: 'W', 115: 'w',
        84: 'A', 116: 'a',
        85: 'A', 117: 'a',
        86: 'B', 118: 'b',
        87: 'S', 119: 's',
        89: 'R', 121: 'r',
    }


class MalformedAssemblyFileError(BaseException):
    pass


num = len
empty = ''

# The following search patterns should ALWAYS return a match,
# which may just be empty strings (intended):
trim_5p_Ns = re.compile("^([nN]*)") 
trim_3p_Ns = re.compile("([nN]*)$")


def _fold(linear_string, width=None):
    if width is None:
        return linear_string

    # TODO: profile whether this, as implemented, is any faster
    #   than stuffing a fixed-length list object

    folded_string = []
    linear_length = len(linear_string)
    for offset in range(0, linear_length, width):
        if offset+width < linear_length:
            folded_string.append(linear_string[offset:offset+width])
        else:
            folded_string.append(linear_string[offset:])
            
    return '\n'.join(folded_string)


def format_fasta(name, sequence, width=None):
    return ">%s\n%s\n" % (name, _fold(sequence, width=width))


def reverse_sequence(sequence):
    if sequence is None:
        return None
    return sequence[::-1].translate(NUCLEOTIDE_COMPLEMENT)


def usage(message=None, exitcode=1, stream=sys.stderr):
    message = '' if message is None else '\nERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <input.assembly> <input.fasta> <out-prefix>\n" % (
        __program__))
    stream.write("\n")
    #                0        10        20        30        40        50        60        70        80
    #                |---------|---------|---------|---------|---------|---------|---------|---------|
    stream.write("Options:\n")
    stream.write("    -c,--cprops-names           Use literal names in .assembly header, don't parse\n")
    stream.write("    -g,--gap-length <uint>      Specify the fixed-length scaffolding gap size [100]\n")
    stream.write("    -N,--trim-terminal-Ns       Trim Ns from contig ends before joining into scaffolds\n")
    stream.write("    -s,--scaffold-prefix <str>  The prefix string for output sequence names [Scaffold]\n")
    stream.write("    -z,--zero-pad-length <uint> Zero-pad scaffold names to <uint> digits [0]\n")
    stream.write("    -h                          This help message\n")
    stream.write("\n%s" % message)
    stream.write("\n")
    sys.exit(exitcode)
    

def main(argv):
    long_options = (
        'help',
        'cprops-names',
        'gap-length=',
        'scaffold-prefix=',
        'zero-pad-length=',
        'trim-terminal-Ns'
    )
    try:
        options, arguments = getopt.getopt(argv, 'hcNg:s:z:', long_options)
    except getopt.GetoptError as error:
        usage(error)

    cprops_names = False
    gap_length = 100
    zero_pad_length = 1
    scaffold_prefix = "Scaffold"
    trim_terminal_Ns = False
    try:
        for flag, value in options:
            if   flag in {'-h','--help'}: usage(exitcode=0)
            elif flag in {'-N','--trim-terminal-Ns'}: trim_terminal_Ns = True
            elif flag in {'-g','--gap-length'}: gap_length = int(value)
            elif flag in {'-s','--scaffold-prefix'}: scaffold_prefix = value
            elif flag in {'-z','--zero-pad-length'}: zero_pad_length = int(value)
            elif flag in {'-c','--cprops-names'}: cprops_names = True
    except ValueError:
        usage("Invalid argument to '%s': %s" % (flag, value))
        
    if gap_length < 0:
        usage("--gap-size argument must be an unsigned integer")
    if zero_pad_length < 0:
        usage("--zero-pad-length argument must be an unsigned integer\n")
    if num(arguments) != 3:
        usage('Unexpected number of arguments')
        
    assembly_file = open(arguments[0], 'r')
    ifasta_file = pysam.FastaFile(arguments[1])
    ofasta_file = open(arguments[2]+'.fasta', 'w')
    ochain_file = open(arguments[2]+'.chain', 'w')

    ifasta_lengths = dict(zip(ifasta_file.references, ifasta_file.lengths))

    gap_sequence = 'N' * gap_length

    prev_name = None
    prev_start = 0
    seen_head = False
    seen_body = False
    cprops_index = dict()
    chain_number = 0
    scaffold_number = 0
    format_scaffold_name = "{0:s}{1:0{2:d}d}".format
    for line in assembly_file:
        if line.isspace() or \
           line.startswith('#'):
            continue

        line = line.strip()
        
        if line.startswith('>'):
            curr_name, curr_index, curr_length = line.lstrip('>').split()
            curr_index = int(curr_index)
            curr_length = int(curr_length)

            if cprops_names:
                orig_name = (curr_name,)
            else:
                orig_name = curr_name.split(':::')

            if orig_name[0] != prev_name:
                prev_start = 0

            cprops_index[curr_index] = (orig_name[0], prev_start, prev_start+curr_length, curr_name)

            prev_name = orig_name[0]
            prev_start += curr_length

            seen_head = True
            
        elif not seen_head:
            raise MalformedAssemblyFileError("No cprops-style assembly file-header detected")

        else:
            scaffold_chain = []
            scaffold_number += 1            
            scaffold_sequence = []
            scaffold_start = 0
            scaffold_end = 0
            scaffold_name = format_scaffold_name(scaffold_prefix, scaffold_number, zero_pad_length)
            scaffold_contig_count = 0
            for vector in line.split():
                vector = int(vector)
                index  = abs(vector)
                strand = vector // index
                
                contig_name, contig_start, contig_end, cprops_name = cprops_index[index]
                
                contig_sequence = ifasta_file.fetch(contig_name, contig_start, contig_end)
                contig_length = contig_end - contig_start
                
                if trim_terminal_Ns:
                    # The following pattern searches should ALWAYS return a match
                    # if not, it's the result of a bug:
                    match_5p = trim_5p_Ns.search(contig_sequence)
                    match_3p = trim_3p_Ns.search(contig_sequence)
                    assert match_5p is not None
                    assert match_3p is not None
                    adjust_5p = len(match_5p.group())
                    adjust_3p = len(match_3p.group())
                    if ((adjust_5p + adjust_3p) >= contig_length):
                        # contig is all Ns
                        continue
                    contig_sequence = contig_sequence[adjust_5p:contig_length-adjust_3p]
                    contig_length -= adjust_5p + adjust_3p
                    contig_start += adjust_5p
                    contig_end -= adjust_3p
                
                if strand < 0:
                    contig_sequence = reverse_sequence(contig_sequence)

                scaffold_end = scaffold_start + contig_end - contig_start

                scaffold_sequence.append(contig_sequence)
                scaffold_chain.append([
                    'chain',
                    0,
                    contig_name,
                    ifasta_lengths[contig_name],
                    '+',
                    contig_start,
                    contig_end,
                    scaffold_name,
                    0,                    
                    '-' if strand < 0 else '+',
                    scaffold_start,
                    scaffold_end,
                    chain_number+1
                ])
                scaffold_contig_count += 1
                scaffold_start = scaffold_end + gap_length
                chain_number += 1
                seen_body = True
                
            scaffold_length = scaffold_end
            for chain in scaffold_chain:
                chain[8] = scaffold_length
                if chain[9] == '-':
                    chain[10], chain[11] = scaffold_length - chain[11], scaffold_length - chain[10]

                ochain_file.write(' '.join(map(str, chain)) + '\n')
                ochain_file.write('%d\n\n' % (chain[6]-chain[5]))
                
            ofasta_file.write(
                format_fasta(
                    scaffold_name,
                    gap_sequence.join(scaffold_sequence),
                    width=100
                )
            )

    if not seen_body:
        raise MalformedAssemblyFileError("No asm-style assembly file-body detected")
            
    assembly_file.close()
    ofasta_file.close()
    ochain_file.close()

main(sys.argv[1:])

