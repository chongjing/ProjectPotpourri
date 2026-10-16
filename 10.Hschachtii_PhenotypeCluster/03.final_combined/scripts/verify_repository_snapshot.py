#!/usr/bin/env python3
"""Verify the curated, public repository snapshot against its manifest."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import xml.etree.ElementTree as ET

TEXT_SUFFIXES={'.csv','.json','.md','.py','.txt','.tsv','.yml','.yaml','.toml','.gitignore',''}
FORBIDDEN_LITERAL=['/home/'+'cx264','C:'+'\\Users\\','github_'+'pat_','gh'+'p_']
TOKEN_PATTERNS=[re.compile(('github_'+'pat_[A-Za-z0-9_]{20,}').encode()),
                re.compile(('gh'+'p_[A-Za-z0-9]{20,}').encode())]


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def verify_repository_snapshot(base,manifest_path=None):
    base=Path(base).resolve()
    manifest_path=Path(manifest_path) if manifest_path else base/'reports/REPOSITORY_MANIFEST.json'
    if not manifest_path.is_absolute(): manifest_path=(base/manifest_path).resolve()
    spec=json.loads(manifest_path.read_text())
    rows=spec['files']
    assert rows,'Empty repository manifest'
    seen=set();pdfs=jpegs=0;total=0
    for row in rows:
        rel=Path(row['path'])
        assert not rel.is_absolute() and '..' not in rel.parts,rel
        assert rel.as_posix() not in seen,rel
        seen.add(rel.as_posix())
        path=(base/rel).resolve()
        assert path.is_relative_to(base),path
        assert path.is_file(),path
        size=path.stat().st_size
        assert size==row['size'],(rel,size,row['size'])
        assert digest(path)==row['sha256'],rel
        assert size<50*1024*1024,('Repository object >=50 MiB',rel,size)
        assert path.suffix.lower() not in {'.tiff','.svg','.log','.pyc','.tag'},rel
        data=path.read_bytes()
        if path.suffix.lower() in TEXT_SUFFIXES:
            text=data.decode('utf-8',errors='ignore')
            for literal in FORBIDDEN_LITERAL:
                assert literal not in text,(rel,literal)
        for pattern in TOKEN_PATTERNS:
            assert not pattern.search(data),('Credential-like token',rel)
        if path.suffix.lower()=='.pdf':
            import fitz
            with fitz.open(path) as doc: assert doc.page_count==1,rel
            pdfs+=1
        elif path.suffix.lower() in {'.jpeg','.jpg'}:
            from PIL import Image
            with Image.open(path) as image:
                assert image.format=='JPEG',rel
                image.verify()
            jpegs+=1
        total+=size
    if (base/'README.md').exists():
        readme=(base/'README.md').read_text()
        previews=sorted(p.relative_to(base).as_posix() for p in (base/'figures/previews').glob('*.jpeg'))
        assert len(previews)==6,previews
        assert all(path in readme for path in previews),'README does not embed every web preview'
    result={'status':'PASS','files_checked':len(rows),'bytes_checked':total,
            'pdfs_checked':pdfs,'jpeg_previews_checked':jpegs,
            'forbidden_path_scan':'PASS','credential_scan':'PASS'}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--manifest',type=Path)
    args=parser.parse_args()
    result=verify_repository_snapshot(args.base,args.manifest)
    print(json.dumps(result,indent=2))
