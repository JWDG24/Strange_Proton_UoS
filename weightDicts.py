import os
import re
import ast
import shutil

# Automatically find wcharm_analysis base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEIGHTS_DIR = os.path.join(BASE_DIR, "weights")

weight_files = [
    "WCPyPartonminus_weights.txt",
    "WCPyPartonplus_weights.txt",
    "weights.hw.particle.txt",
    "weights.py.particle.txt",
    "WCHPartonminus_weights.txt",
    "WCHPartonplus_weights.txt",
]

# Fallback mapping for Herwig parton
fallback_map = {
    "WCHPartonminus_weights.txt": "WCPyPartonminus_weights.txt",
    "WCHPartonplus_weights.txt": "WCPyPartonplus_weights.txt",
}

for filename in weight_files:

    input_path = os.path.join(WEIGHTS_DIR, filename)

    # ------------------------------------------------------------
    # If file missing and fallback defined → use fallback file
    # ------------------------------------------------------------
    if not os.path.exists(input_path):

        if filename in fallback_map:
            fallback_file = fallback_map[filename]
            fallback_path = os.path.join(WEIGHTS_DIR, fallback_file)

            if os.path.exists(fallback_path):
                print(f"{filename} not found. Using fallback {fallback_file}")
                input_path = fallback_path
            else:
                print("Fallback file also not found:", fallback_file)
                continue
        else:
            print("File not found:", filename)
            continue

    with open(input_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()

    # Extract dictionary block
    match = re.search(r"(\{.*\})", text, re.S)
    if not match:
        print("No dictionary found in:", filename)
        continue

    raw_dict = ast.literal_eval(match.group(1))

    # Clean keys
    clean_dict = {str(k).strip(): int(v) for k, v in raw_dict.items()}

    # Sort by index
    sorted_items = sorted(clean_dict.items(), key=lambda x: x[1])

    # Output clean file name (always based on original filename)
    clean_name = filename.replace(".txt", "_clean.txt")
    output_path = os.path.join(WEIGHTS_DIR, clean_name)

    with open(output_path, "w") as out:
        for name, idx in sorted_items:
            out.write(f"{idx} {name}\n")

    print("Created:", clean_name)

print("Done.")
