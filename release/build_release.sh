#!/usr/bin/bash -e

BUILD_DIR=/data/openpilot
cd $BUILD_DIR
rm -rf .git
git init
git remote add origin https://gitcode.com/fishop/openpilot.git

# in the directory
cd $BUILD_DIR

# Cleanup
find . -name '*.a' -delete
find . -name '*.o' -delete
find . -name '*.os' -delete
find . -name '*.pyc' -delete
find . -name 'moc_*' -delete
find . -name '__pycache__' -delete
rm -rf .sconsign.dblite Jenkinsfile release/
#rm -f openpilot/selfdrive/modeld/models/*.onnx
# drop the legacy stamp inside modeld/; it is regenerated at repo root below
rm -f openpilot/selfdrive/modeld/models/.build_stamp

# ship the prebuilt release WITHOUT the .onnx model inputs (the tinygrad .pkl
# artifacts + prebuilt marker are enough). Keep the files on disk so the device
# can still rebuild if ever needed; just don't commit them.
echo 'openpilot/selfdrive/modeld/models/*.onnx' > .gitignore

find third_party/ -name '*x86*' -exec rm -r {} +
find third_party/ -name '*Darwin*' -exec rm -r {} +

# Mark as prebuilt release
touch prebuilt

# Add built files to git
git add -f .

VERSION="carrot_v$(date +%y%m%d)"
git commit -m $VERSION
git branch -m "egpucp"

# Recompute .build_stamp against the exact HEAD that will be pushed.
# launch_chffrplus.sh compares this stamp on every boot; if it doesn't match
# the pushed HEAD, FORCE_REBUILD=1 and the device tries to recompile the
# models. But the .onnx inputs are not shipped in this prebuilt release, so the
# build fails and the device hangs. Regenerate the stamp here so the release
# always matches its own HEAD.
# NOTE: the stamp value is the git tree hash of openpilot/selfdrive/modeld, so
# the stamp file itself MUST live outside that tree (repo root). If it lived
# inside modeld/, changing it would change the very hash it records, and the
# stamp could never match after commit.
STAMP="$(git rev-parse HEAD:openpilot/selfdrive/modeld HEAD:tinygrad_repo HEAD:openpilot/common/file_chunker.py | tr '\n' ':')"
echo -n "$STAMP" > .build_stamp
git add -f .build_stamp
git commit -m "${VERSION}-stamp"

git push -f origin "egpucp"
