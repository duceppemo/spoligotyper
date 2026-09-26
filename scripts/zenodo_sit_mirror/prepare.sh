#!/usr/bin/env bash
# Prepare the files of the Zenodo mirror of the SITVIT2 pattern list published with SpolLineages (GPL-3.0).
# The list itself is not stored in the spoligotyper repository (MIT): this script downloads the version pinned in
# spoligotyper/sitdb.py, checks its SHA-256 checksum, and adds the GPL-3.0 license and a README.
# Usage: bash prepare.sh [output_folder]   Then upload the files of the folder to Zenodo (see README.md).
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
out=${1:-sit_mirror}
read -r commit sha256 < <(python3 -c "
import re, sys
text = open('$here/../../spoligotyper/sitdb.py').read()
print(re.search(r\"'commit': '(\w+)'\", text).group(1), re.search(r\"'sha256': '(\w+)'\", text).group(1))")

mkdir -p "$out"
curl -sSfL -o "$out/Spoligo_list.csv" "https://raw.githubusercontent.com/dcouvin/SpolLineages/$commit/Spoligo_list.csv"
curl -sSfL -o "$out/LICENSE" "https://raw.githubusercontent.com/dcouvin/SpolLineages/$commit/LICENSE"
echo "$sha256  $out/Spoligo_list.csv" | sha256sum -c -
sed -e "s/@COMMIT@/$commit/g" -e "s/@SHA256@/$sha256/g" "$here/README_mirror.md" > "$out/README.md"
echo "Files ready in $out/: $(ls "$out" | tr '\n' ' ')"
