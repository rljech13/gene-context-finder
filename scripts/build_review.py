#!/usr/bin/env python3
"""Build optional portable launchers from the same reviewed Go source.

No precompiled executable is needed for the Python package or ordinary reports.
Go is an explicit build dependency; no scientific tool or model is downloaded.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess

PLATFORMS=('darwin-arm64','darwin-amd64','linux-arm64','linux-amd64','windows-arm64','windows-amd64')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--platform',choices=('current','all',*PLATFORMS),default='current');p.add_argument('--go',default='go');args=p.parse_args()
    source=Path(__file__).resolve().parent.parent/'portable-review';out=args.out.resolve()
    if out==source or out in source.parents or source in out.parents:raise SystemExit('Choose a separate build output directory')
    choice=args.platform
    if choice=='current':
        choice={'Darwin':'darwin','Linux':'linux','Windows':'windows'}[platform.system()]+'-'+{'arm64':'arm64','aarch64':'arm64','x86_64':'amd64','AMD64':'amd64'}[platform.machine()]
    targets=PLATFORMS if choice=='all' else (choice,)
    files={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in source.iterdir() if x.suffix in {'.go','.json','.js','.mod'}}
    version=subprocess.check_output([args.go,'version'],text=True).strip();out.mkdir(parents=True,exist_ok=True);binaries={}
    for target in targets:
        osname,arch=target.split('-');name='review-'+target+('.exe' if osname=='windows' else '')
        dest=out/name
        if dest.exists():raise SystemExit('Use a fresh build directory: '+str(dest))
        subprocess.run([args.go,'build','-trimpath','-ldflags=-s -w','-o',str(dest),'.'],cwd=source,env=dict(os.environ,GOOS=osname,GOARCH=arch,CGO_ENABLED='0'),check=True)
        binaries[name]=hashlib.sha256(dest.read_bytes()).hexdigest()
    shutil.copyfile(source.parent/'genecontext/vendor/review/LICENSE-Go.txt',out/'LICENSE-Go.txt')
    (out/'build.json').write_text(json.dumps(dict(compiler=version,source_sha256=files,binaries_sha256=binaries),indent=2)+'\n')
    print(json.dumps(dict(runtime=str(out),platforms=list(targets),source_sha256=files)))

if __name__=='__main__':main()
