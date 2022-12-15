#!/usr/bin/python3

"""

DESCRIPTION:
    Template code for the Dynamic Programming assignment in the Algorithms in Sequence Analysis course at the VU.
    

run : python3 align.py -v -s (global,local,semiglobal) -m (pam250,blosum62,identity) -g 2 test.fasta output.txt

AUTHOR:
    <Name/Surname :THEODOROS FOSKOLOS, studentID: 2768082>
"""



import argparse
import pickle



def parse_args():
    "Parses inputs from commandline and returns them as a Namespace object."

    parser = argparse.ArgumentParser(prog = 'python3 align.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Aligns the first two sequences in a specified FASTA\n'
        '  file with a chosen strategy and parameters.\n'
        '\n'
        'defaults:\n'
        '  strategy = global\n'
        '  substitution matrix = pam250\n'
        '  gap penalty = 2')
        
    parser.add_argument('fasta', help='path to a FASTA formatted input file')
    parser.add_argument('output', nargs='*', 
        help='path to an output file where the alignment is saved\n'
             '  (if a second output file is given,\n'
             '   save the score matrix in there)')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
        help='print the score matrix and alignment on screen', default=False)
    parser.add_argument('-s', '--strategy', dest='strategy',
        choices=['global','semiglobal','local'], default="global")
    parser.add_argument('-m', '--matrix', dest='substitution_matrix',
        choices=['pam250','blosum62','identity'], default='pam250')
    parser.add_argument('-g', '--gap_penalty', dest='gap_penalty', type=int,
        help='must be a positive integer', default=2)

    args = parser.parse_args()

    args.align_out = args.output[0] if args.output else False
    args.matrix_out = args.output[1] if len(args.output) >= 2 else False
                      # Fancy inline if-else statements. Use cautiously!
                      
    if args.gap_penalty <= 0:
        parser.error('gap penalty must be a positive integer')

    return args



def load_substitution_matrix(name):
    "Loads and returns the specified substitution matrix from a pickle (.pkl) file."
    # Substitution matrices have been prepared as nested dictionaries:
    # the score of substituting A for Z can be found with subst['A']['Z']
    # NOTE: Only works if working directory contains the correct folder and file!
    
    with open('substitution_matrices/%s.pkl' % name, 'rb') as f:
        subst = pickle.load(f)
    return subst
    
    

def load_sequences(filepath):
    "Reads a FASTA file and returns the first two sequences it contains."
    
    seq1 = []
    seq2 = []
    with open(filepath,'r') as f:
        for line in f:
            if line.startswith('>'):
                if not seq1:
                    current_seq = seq1
                elif not seq2:
                    current_seq = seq2
                else:
                    break # Stop if a 3rd sequence is encountered
            else:
                current_seq.append(line.strip())
    
    if not seq2:
        raise Exception('Error: Not enough sequences in specified FASTA file.')
    
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2



def align(seq1, seq2, strategy, substitution_matrix, gap_penalty):
    "Do pairwise alignment using the specified strategy and parameters."
    # This function consists of 3 parts:
        
    #   1) Initialize a score matrix as a "list of lists" of the appropriate length.
    #      Fill in the correct values for the first row and column given the strategy.
    #        (local / semiglobal = 0  --  global = stacking gap penalties)
    #   2) Fill in the rest of the score matrix using Dynamic Programming, accounting
    #      for the selected alignment strategy, substitution matrix and gap penalty.
    #   3) Perform the correct traceback routine on your filled in score matrix.
    #
    # Both the resulting alignment (sequences with gaps and the corresponding score)
    # and the filled in score matrix are returned as outputs.
    #
    # NOTE: You are strongly encouraged to think about how you can reuse (parts of)
    #       your code between steps 2 and 3 for the different strategies!
    
    
    ### 1: Initialize
    M = len(seq1)+1
    N = len(seq2)+1
    score_matrix = []
    
    # put zeros into the matrix
    for i in range(M):
        row = []
        score_matrix.append(row)
        for j in range(N):
            row.append(0)
    
    #####################
    # START CODING HERE #
    #####################
    # initialize the global matrix 1st row/column with gap penalties 
    if strategy == 'global':
        
        for i in range(1,M):
            score_matrix[i][0] = score_matrix[i-1][0] - gap_penalty
            
        for j in range(1,N):
                        
            score_matrix[0][j] = score_matrix[0][j-1] - gap_penalty
        
        #####################
        #  END CODING HERE  #
        #####################

    
    
    ### 2: Fill in Score Matrix
 
    #####################
    # START CODING HERE #
    #####################
    def dp_function():
        
        # Calculate the different scores from different positions
        diag = (substitution_matrix[seq1[i-1]][seq2[j-1]]) + (score_matrix[i-1][j-1])
        vertical = score_matrix[i-1][j] - gap_penalty
        horizontal = score_matrix[i][j-1] - gap_penalty
        
        if strategy == 'global':
            # Get the max score for global
            score = max(diag, horizontal, vertical)
               
        elif strategy == 'local':
            
            # local strategy make negative numbers 0 and then take the max for the score
            
            if diag < 0:
           
                diag = 0
                
            if horizontal < 0:
               
                horizontal = 0
                
            if vertical <0:
                
                vertical = 0 
            
            score = max(diag, horizontal, vertical)
            
        
        elif strategy == 'semiglobal':
            score =  max(diag, horizontal, vertical)
            

        return score
    

    # Call the function and fill the matrix with the correct values
    for i in range(1,M):
        for j in range(1,N):
            score_matrix[i][j] = dp_function()
       
      
    #####################
    #  END CODING HERE  #
    #####################   
    
    
    ### 3: Traceback
    
    #####################
    # START CODING HERE #
    ##################### 
    
    aligned_seq1 = ""
    aligned_seq2 = ""
    align_score = 0 
    
    if strategy == 'global':
        # Bottom right value gives us the score for Global
        align_score = score_matrix[M-1][N-1]
        i = M-1
        j = N-1
        starter = score_matrix[M-1][N-1]
       
        # GLOBAL Traceback/ alignments
        while  j!=0 or i!=0:
            
           # Vertical traceback/ gap
            if starter == score_matrix[i-1][j] - gap_penalty and i != 0:
                aligned_seq1 =seq1[i-1] + aligned_seq1 # -1 because M and N len seq + 1
                aligned_seq2 =  "-" + aligned_seq2 
                starter = score_matrix[i-1][j]
                i = i - 1
                
                            
            # Diagonial traceback
            elif starter == substitution_matrix[seq1[i-1]][seq2[j-1]] + score_matrix[i-1][j-1] and (i and j)!= 0:
                aligned_seq1 = seq1[i-1] + aligned_seq1 
                aligned_seq2 = seq2[j-1] + aligned_seq2 
                starter = score_matrix[i-1][j-1]
                i = i - 1
                j = j - 1
          
                          
            # Horizontal traceback/ gap   
            elif starter == score_matrix[i][j-1] - gap_penalty and j!= 0 :
                aligned_seq1 =  "-" + aligned_seq1 
                aligned_seq2 = seq2[j-1] + aligned_seq2 
                starter = score_matrix[i][j-1]
                j = j-1
              
                
            
    
    if strategy == 'semiglobal':  # Find the max in the most right column or last line of the matrix/ SEMIGLOBAL
        starter = score_matrix[M-1][N-1]
        maxim= -1
        positioncol , positionrow= 0, 0
        pos = ""
        
        # Search for the max in last row                    
        for j in range (1, N):
            
            if maxim <= score_matrix[M-1][j]:
                maxim = score_matrix[M-1][j]
                positioncol = j
                positionrow = M-1
                pos = "lastrow"
                
        # Search for the max in last column
        for i in range (M-1,-1,-1):
                
            if maxim <= score_matrix[i][N-1] :
                maxim = score_matrix[i][N-1]
                positionrow = i
                positioncol = N-1
                pos = "lastcol"
                
                
        # If the max is on the last column, put gaps from bottom right till you reach it
        if pos == "lastcol":
        
            for gap in range (M-1, positionrow, -1):
                aligned_seq1 =  seq1[gap-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
     
        # If the max is on the last row, put gaps from bottom right till you reach it           
        if pos == "lastrow":
            for gap in range (N-1, positioncol, -1):
                aligned_seq1 =  "-" + aligned_seq1
                aligned_seq2 = seq2[gap-1] + aligned_seq2
                
    
        align_score = maxim # max value gives us the alignment score

        
       # Ends when you reach far left, starting from the max/ SEMIGLOBAL TRACEBACK
        while  (positioncol !=0) or (positionrow != 0) :
            
              # If both row and column are higher than 0 do regular traceback
              if  positioncol !=0 and positionrow!= 0:
                  
                  
                  # Vertical traceback / gap
                  if maxim == score_matrix [positionrow -1][positioncol] - gap_penalty and positionrow != 0:
                      aligned_seq1 = seq1 [positionrow-1] + aligned_seq1 
                      aligned_seq2 =  "-" + aligned_seq2 
                      maxim = score_matrix[positionrow-1][positioncol]
                      positionrow = positionrow - 1
    
                                   
                  # Diagonial traceback / gap
                  elif maxim == substitution_matrix[seq1[positionrow-1]][seq2[positioncol-1]] + score_matrix [positionrow -1][positioncol-1] :
                      aligned_seq1 = seq1[positionrow-1] + aligned_seq1 
                      aligned_seq2 = seq2[positioncol-1] + aligned_seq2 
                      maxim = score_matrix[positionrow-1][positioncol-1]
                      positionrow = positionrow - 1
                      positioncol = positioncol - 1
                  
                                  
                  # Horizontal traceback / gap    
                  elif maxim == score_matrix[positionrow][positioncol-1] - gap_penalty and positioncol != 0:
                      aligned_seq1 =  "-" + aligned_seq1 
                      aligned_seq2 = seq2[positioncol-1] + aligned_seq2 
                      maxim = score_matrix[positionrow][positioncol-1]
                      positioncol = positioncol-1 
                  
                    
              # When we reach  the 0th row or column just gap till you reach the upper left scorebox   and then break the loop            
              elif positionrow == 0:
                  for gap in range (positioncol, 0, -1):
                      aligned_seq1 =  "-" + aligned_seq1
                      aligned_seq2 = seq2[gap-1] + aligned_seq2
                      
                  break
                        
                       
              elif positioncol == 0:
                  for gap in range (positionrow, 0, -1):
                      aligned_seq1 = seq1[gap-1] + aligned_seq1
                      aligned_seq2 = "-" + aligned_seq2
                      
                  break
                
                      
    # Find the max anywhere in the matrix / LOCAL TRACEBACK
    if strategy == 'local':
    
        maxim= 0
        maxrow = 0
        maxcol = 0
        # Search for the max
        for row in range (1,M):
            for column in range (1,N):
                if maxim < score_matrix[row][column]:
                    maxim = score_matrix[row][column]
                    maxrow = row
                    maxcol = column
                elif maxim == score_matrix[row][column] and column > maxcol:
                    maxim = score_matrix[row][column]
                    maxrow = row
                    maxcol = column
                    
        # The max value equals  the alignment score         
        align_score = maxim
        
        # Ends if it is in the first row/column with zeros        
        while  (maxrow and maxcol) !=0:
            
            # Vertical traceback/ gap
            if maxim == score_matrix[maxrow-1][maxcol] - gap_penalty :
                aligned_seq1 =seq1[maxrow-1] + aligned_seq1 
                aligned_seq2 =  "-" + aligned_seq2 
                maxim = score_matrix[maxrow-1][maxcol]
                maxrow = maxrow - 1
                
                           
           # Diagonial traceback/ if it is 0 i need to check if it comes from a real 0 and not from a negative number otherwise i have a break on the last else
            elif maxim == substitution_matrix[seq1[maxrow-1]][seq2[maxcol-1]] + score_matrix[maxrow-1][maxcol-1] or (maxim ==0 and substitution_matrix[seq1[maxrow-1]][seq2[maxcol-1]] + score_matrix[maxrow-1][maxcol-1]==0) : #or maxim == 0
                aligned_seq1 = seq1[maxrow-1] + aligned_seq1 
                aligned_seq2 = seq2[maxcol-1] + aligned_seq2 
                maxim = score_matrix[maxrow-1][maxcol-1]
                maxrow = maxrow - 1
                maxcol = maxcol - 1
                
                                                   
            # Horizontal traceback / gap        
            elif maxim == score_matrix[maxrow][maxcol-1] - gap_penalty:
                aligned_seq1 =  "-" + aligned_seq1 
                aligned_seq2 = seq2[maxcol-1] + aligned_seq2 
                maxim = score_matrix[maxrow][maxcol-1]
                maxcol = maxcol-1
                
            # Break the loop if none of the above exist
            else: 
                break

    #####################
    #  END CODING HERE  #
    #####################   

    
    alignment = (aligned_seq1, aligned_seq2, align_score)
    
    return (alignment, score_matrix)



def print_score_matrix(s1,s2,mat):
    "Pretty print function for a score matrix."
    
    # Prepend filler characters to seq1 and seq2
    s1 = '-' + s1
    s2 = ' -' + s2
    
    # Print them around the score matrix, in columns of 5 characters
    print(''.join(['%5s' % aa for aa in s2])) # Convert s2 to a list of length 5 strings, then join it back into a string
    for i,row in enumerate(mat):               # Iterate through the rows of your score matrix (and keep count with 'i').
        vals = ['%5i' % val for val in row]    # Convert this row's scores to a list of strings.
        vals.insert(0,'%5s' % s1[i])           # Add this row's character from s2 to the front of the list
        print(''.join(vals))                   # Join the list elements into a single string, and print the line.



def print_alignment(a):
    "Pretty print function for an alignment (and alignment score)."
    
    # Unpack the alignment tuple
    seq1 = a[0]
    seq2 = a[1]
    score = a[2]
    
    # Check which positions are identical
    match = ''
    for i in range(len(seq1)): # Remember: Aligned sequences have the same length!
        match += '|' if seq1[i] == seq2[i] else ' ' # Fancy inline if-else statement. Use cautiously!
            
    # Concatenate lines into a list, and join them together with newline characters.
    print('\n'.join([seq1,match,seq2,'','Score = %i' % score]))



def save_alignment(a,f):
    "Saves two aligned sequences and their alignment score to a file."
    with open(f,'w') as out:
        out.write(a[0] + '\n') # Aligned sequence 1
        out.write(a[1] + '\n') # Aligned sequence 2
        out.write('Score: %i' % a[2]) # Alignment score


    
def save_score_matrix(m,f):
    "Saves a score matrix to a file in tab-separated format."
    with open(f,'w') as out:
        for row in m:
            vals = [str(val) for val in row]
            out.write('\t'.join(vals)+'\n')
    


def main(args = False):
    # Process arguments and load required data
    if not args: args = parse_args()
    
    sub_mat = load_substitution_matrix(args.substitution_matrix)
    seq1, seq2 = load_sequences(args.fasta)

    # Perform specified alignment
    strat = args.strategy
    gp = args.gap_penalty
    alignment, score_matrix = align(seq1, seq2, strat, sub_mat, gp)

    # If running in "verbose" mode, print additional output
    if args.verbose:
        print_score_matrix(seq1,seq2,score_matrix)
        print('') # Insert a blank line in between
        print_alignment(alignment)
    
    # Save results
    if args.align_out: save_alignment(alignment, args.align_out)
    if args.matrix_out: save_score_matrix(score_matrix, args.matrix_out)



if __name__ == '__main__':
    main()