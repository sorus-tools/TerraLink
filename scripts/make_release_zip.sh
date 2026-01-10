set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT_ZIP="TerraLink.zip"
STAGE="$(mktemp -d)"

cleanup() { rm -rf "$STAGE"; }
trap cleanup EXIT

rsync -a --delete \
  --exclude ".git" \
  --exclude "__MACOSX" \
  --exclude "__pycache__" \
  --exclude ".DS_Store" \
  --exclude "TerraLink.zip" \
  --exclude ".venv" \
  --exclude "venv" \
  --exclude ".idea" \
  --exclude ".vscode" \
  "$ROOT/" "$STAGE/TerraLink/"

find "$STAGE/TerraLink" -name "__pycache__" -type d -prune -exec rm -rf {} +
find "$STAGE/TerraLink" -name ".DS_Store" -type f -delete

rm -f "$ROOT/$OUT_ZIP"
(cd "$STAGE" && zip -r "$ROOT/$OUT_ZIP" "TerraLink" >/dev/null)

echo "$OUT_ZIP created at $ROOT/$OUT_ZIP"
