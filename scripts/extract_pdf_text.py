
from pdfminer.high_level import extract_text
import sys

try:
    text = extract_text(r"D:\PersonalFile\Documents\BigProject\s13321-017-0255-6.pdf")
    # Search for keywords
    keywords = ["Score", "Weighted", "Normalized", "Formula", "molecular_weight"]
    
    print("=== Extracted Text Preview ===")
    print(text[:2000])
    
    print("\n=== Keyword Search ===")
    lines = text.split('\n')
    for i, line in enumerate(lines):
        if any(k.lower() in line.lower() for k in keywords):
            # Print context
            start = max(0, i-2)
            end = min(len(lines), i+3)
            print(f"--- Line {i} ---")
            for j in range(start, end):
                 print(lines[j])
except Exception as e:
    print(f"Error: {e}")
