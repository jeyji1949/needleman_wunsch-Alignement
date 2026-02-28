import numpy as np

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    
    # Create score matrix
    S = np.zeros((n+1, m+1), dtype=int)
    
    # Initialize first row and column with gap penalties
    for i in range(1, n+1):
        S[i][0] = S[i-1][0] + gap
    for j in range(1, m+1):
        S[0][j] = S[0][j-1] + gap
    
    # Fill the matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = S[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up = S[i-1][j] + gap
            left = S[i][j-1] + gap
            S[i][j] = max(diag, up, left)
    
    # Print the matrix
    print("\nSCORE MATRIX:")
    print("    ", end="")
    for j in range(m):
        print(f"  {seq2[j]}", end="")
    print()
    for i in range(n+1):
        if i == 0:
            print("  ", end="")
        else:
            print(f"{seq1[i-1]} ", end="")
        for j in range(m+1):
            print(f"{S[i][j]:3}", end="")
        print()
    
    # Backtrace
    i, j = n, m
    align1, align2 = [], []
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and S[i][j] == S[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and S[i][j] == S[i-1][j] + gap:
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
    
    return S[n][m], align1, align2

# Main program
print("NEEDLEMAN-WUNSCH ALIGNMENT")
print("-" * 30)

# Get user input
seq1 = input("Enter first sequence: ").upper()
seq2 = input("Enter second sequence: ").upper()

# Run alignment
score, align1, align2 = needleman_wunsch(seq1, seq2)

# Print results
print("\n" + "=" * 40)
print("RESULTS:")
print("=" * 40)
print(f"Sequence 1: {seq1}")
print(f"Sequence 2: {seq2}")
print(f"\nAlignment score: {score}")
print(f"\nOptimal alignment:")
print(f"Seq1: {align1}")
print(f"Seq2: {align2}")

# Show match/mismatch line
match_line = ""
for i in range(len(align1)):
    if align1[i] == align2[i]:
        match_line += "|"
    elif align1[i] == "-" or align2[i] == "-":
        match_line += " "
    else:
        match_line += "*"
print(f"      {match_line}")