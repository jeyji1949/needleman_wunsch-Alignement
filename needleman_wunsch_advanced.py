import numpy as np

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    
    # Create score matrix with headers
    # Create a matrix with extra row and column for gaps
    S = [[0]*(m+1) for _ in range(n+1)]
    
    # Initialize first row and column with gap penalties
    for i in range(1, n+1):
        S[i][0] = S[i-1][0] + gap
    for j in range(1, m+1):
        S[0][j] = S[0][j-1] + gap
    
    # Track backtrace directions (for visualization)
    # D: diagonal, U: up, L: left
    trace = [['' for _ in range(m+1)] for _ in range(n+1)]
    trace[0][0] = '0'
    for i in range(1, n+1):
        trace[i][0] = 'U'
    for j in range(1, m+1):
        trace[0][j] = 'L'
    
    # Fill the matrix
    print("\n" + "="*60)
    print("STEP-BY-STEP MATRIX FILLING:")
    print("="*60)
    
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = S[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up = S[i-1][j] + gap
            left = S[i][j-1] + gap
            
            # Find max value
            max_val = max(diag, up, left)
            S[i][j] = max_val
            
            # Record trace direction(s) - for simplicity, we record the first max found
            if max_val == diag:
                trace[i][j] = 'D'
            elif max_val == up:
                trace[i][j] = 'U'
            else:
                trace[i][j] = 'L'
            
            # Print current step
            print(f"\nCell [{i},{j}] - Comparing '{seq1[i-1]}' with '{seq2[j-1]}':")
            print(f"  Diagonal ({seq1[i-1]},{seq2[j-1]}): {S[i-1][j-1]} + {match if seq1[i-1]==seq2[j-1] else mismatch} = {diag}")
            print(f"  Up (gap in seq2): {S[i-1][j]} + {gap} = {up}")
            print(f"  Left (gap in seq1): {S[i][j-1]} + {gap} = {left}")
            print(f"  → Selected: {max_val} (direction: {trace[i][j]})")
    
    # Print final matrix with sequences as headers
    print("\n" + "="*60)
    print("FINAL SCORE MATRIX:")
    print("="*60)
    
    # Create headers
    headers = [' '] + list(seq2)
    
    # Print column headers
    print("\n     " + " ".join(f"{h:4}" for h in headers))
    
    # Print each row with row header
    print("     " + "-" * (5 * (m+1)))
    for i in range(n+1):
        if i == 0:
            row_header = ' '
        else:
            row_header = seq1[i-1]
        print(f"{row_header:2} | " + " ".join(f"{S[i][j]:4}" for j in range(m+1)))
    
    # Print trace matrix
    print("\nTRACE MATRIX (D: Diagonal, U: Up, L: Left, 0: Start):")
    print("     " + " ".join(f"{h:4}" for h in headers))
    print("     " + "-" * (5 * (m+1)))
    for i in range(n+1):
        if i == 0:
            row_header = ' '
        else:
            row_header = seq1[i-1]
        print(f"{row_header:2} | " + " ".join(f"{trace[i][j]:4}" for j in range(m+1)))
    
    # Backtrace to find optimal alignment
    print("\n" + "="*60)
    print("BACKTRACE PATH:")
    print("="*60)
    
    i, j = n, m
    path = []
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace[i][j] == 'D':
            path.append(f"({i},{j}) - Diagonal: {seq1[i-1]} ↔ {seq2[j-1]}")
            i -= 1
            j -= 1
        elif i > 0 and trace[i][j] == 'U':
            path.append(f"({i},{j}) - Up: {seq1[i-1]} ↔ -")
            i -= 1
        else:
            path.append(f"({i},{j}) - Left: - ↔ {seq2[j-1]}")
            j -= 1
    
    # Reverse to show from start to end
    for step in reversed(path):
        print(step)
    
    # Reconstruct alignment
    i, j = n, m
    align1, align2 = [], []
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace[i][j] == 'D':
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and trace[i][j] == 'U':
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1
    
    # Reverse alignments
    align1 = ''.join(reversed(align1))
    align2 = ''.join(reversed(align2))
    
    # Calculate score (should match S[n][m])
    alignment_score = S[n][m]
    
    return alignment_score, align1, align2, S, path

def display_alignment(seq1, seq2, score, match=1, mismatch=-1, gap=-2):
    """Display the alignment with match/mismatch indicators"""
    print("\n" + "="*60)
    print("OPTIMAL ALIGNMENT:")
    print("="*60)
    
    # Create match/mismatch line
    match_line = []
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            match_line.append('|')  # Match
        elif seq1[i] == '-' or seq2[i] == '-':
            match_line.append(' ')  # Gap
        else:
            match_line.append('*')  # Mismatch
    
    # Print alignment
    print(f"Seq1: {seq1}")
    print(f"      {''.join(match_line)}")
    print(f"Seq2: {seq2}")
    print(f"\nAlignment Score: {score}")
    
    # Calculate statistics
    matches = match_line.count('|')
    mismatches = match_line.count('*')
    gaps = seq1.count('-') + seq2.count('-')
    
    print(f"\nStatistics:")
    print(f"  Length: {len(seq1)}")
    print(f"  Matches: {matches}")
    print(f"  Mismatches: {mismatches}")
    print(f"  Gaps: {gaps}")
    print(f"  Identity: {matches/len(seq1)*100:.1f}%")

def main():
    print("NEEDLEMAN-WUNSCH SEQUENCE ALIGNMENT")
    print("="*50)
    
    # Get user input
    seq1 = input("Enter first sequence (e.g., ACCG): ").strip().upper()
    seq2 = input("Enter second sequence (e.g., ACTG): ").strip().upper()
    
    # Validate input (optional - can be expanded)
    valid_bases = set('ACGT')
    if not all(base in valid_bases for base in seq1):
        print("Warning: Sequence 1 contains non-standard bases (only A,C,G,T recommended)")
    if not all(base in valid_bases for base in seq2):
        print("Warning: Sequence 2 contains non-standard bases (only A,C,G,T recommended)")
    
    # Ask for scoring parameters (optional)
    print("\nDefault scoring: Match=1, Mismatch=-1, Gap=-2")
    use_default = input("Use default scoring? (y/n): ").strip().lower()
    
    if use_default == 'n':
        match = int(input("Match score: "))
        mismatch = int(input("Mismatch penalty: "))
        gap = int(input("Gap penalty: "))
    else:
        match, mismatch, gap = 1, -1, -2
    
    # Run Needleman-Wunsch
    score, align1, align2, matrix, path = needleman_wunsch(seq1, seq2, match, mismatch, gap)
    
    # Display results
    display_alignment(align1, align2, score, match, mismatch, gap)

if __name__ == "__main__":
    main()