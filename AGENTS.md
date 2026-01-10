# TerraLink repo operating rules

## Branching
- Never commit directly to main.
- Use branches:
  - feature/<short>
  - fix/<short>
  - release/vX.Y.Z

## Quality checks (run before commit)
- Ensure metadata.txt version matches release version when releasing.
- Ensure no accidental large binaries or secrets are added.
- If a change affects UX or usage, update readme.md.

## Release policy (GitHub + QGIS)
- Release artifact is TerraLink.zip with a single top-level folder: TerraLink/
- ZIP must NOT include: .git, __MACOSX, __pycache__, .DS_Store
- Create an annotated tag: vX.Y.Z
- Create GitHub Release for the tag and upload TerraLink.zip
